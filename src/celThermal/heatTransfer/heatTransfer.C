/*---------------------------------------------------------------------------*\
    Temperature in the granular material (WP T.1, decision 0136).
    See heatTransfer.H.
\*---------------------------------------------------------------------------*/

#include "heatTransfer.H"
#include "fvMesh.H"
#include "Time.H"
#include "emptyFvPatch.H"
#include "poreGas.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

Foam::IOobject thermalPropertiesIO(const Foam::fvMesh& mesh)
{
    return Foam::IOobject
    (
        "thermalProperties",
        mesh.time().constant(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
    );
}

//- Zehner and Schluender (1970) stagnant bed conductivity [W/(m K)], as
//  split between the phases by Kuipers et al. (1992); e gas fraction
//  (docs/spec/thermal-references.md, tools/thermal_reference.py zs_bed)
Foam::scalar zehnerSchlunder
(
    const Foam::scalar e,
    const Foam::scalar kGas,
    const Foam::scalar kPart
)
{
    using namespace Foam;
    const scalar omega = 7.26e-3;
    const scalar A = kPart/kGas;
    scalar B = 1.25*Foam::pow((1 - e)/e, 10.0/9.0);
    scalar C = 1;
    if (mag(A - 1) > 1e-12)
    {
        if (mag(1 - B/A) < 1e-9)
        {
            B *= 1 + 1e-6;      // removable singularity at B = A
        }
        const scalar r = B/A;
        C = 2/(1 - r)
           *(
                (A - 1)*r/sqr(1 - r)*Foam::log(A/B)
              - (B - 1)/(1 - r)
              - 0.5*(B + 1)
            );
    }
    const scalar s = Foam::sqrt(1 - e);
    return ((1 - s) + s*(omega*A + (1 - omega)*C))*kGas;
}


//- The two parts of Zehner-Schluender as split by Kuipers et al. (1992):
//  e k_g,eff (gas) and (1 - e) k_s,eff (grains) [W/(m K)]
void zehnerSchlunderParts
(
    const Foam::scalar e,
    const Foam::scalar kGas,
    const Foam::scalar kPart,
    Foam::scalar& gasPart,
    Foam::scalar& solidPart
)
{
    using namespace Foam;
    const scalar s = Foam::sqrt(1 - e);
    gasPart = (1 - s)*kGas;
    solidPart = zehnerSchlunder(e, kGas, kPart) - gasPart;
}


//- Face conductivity from the two cells' values (harmonic, series)
inline Foam::scalar faceK(const Foam::scalar kP, const Foam::scalar kN)
{
    return kP == kN ? kP : (kP + kN > 0 ? 2*kP*kN/(kP + kN) : 0);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransfer::heatTransfer
(
    eulerianSolid& solid,
    const IOobject& io,
    poreGas* gas
)
:
    solid_(solid),
    gas_(gas),
    mesh_(solid.mesh()),
    dict_(io),
    model_(dict_.get<word>("model")),
    c_(dict_.get<scalar>("c")),
    conductivity_(dict_.getOrDefault<word>("conductivity", "constant")),
    k_(dict_.getOrDefault<scalar>("k", 0)),
    kGas_(dict_.getOrDefault<scalar>("kGas", 0.02569)),
    kParticle_(dict_.getOrDefault<scalar>("kParticle", 0)),
    rhoP_(dict_.getOrDefault<scalar>("rhoP", 0)),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    heatIn_(0),
    energy0_(0),
    carriedIn_(0),
    TOld_(),
    limiter_(dict_.getOrDefault<word>("limiter", "MC")),
    voidT_(dict_.getOrDefault<scalar>("voidT", 293.15)),
    kCell_
    (
        IOobject
        (
            "kBed",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            conductivity_ == "zehnerSchlunder"
          ? IOobject::AUTO_WRITE
          : IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    wallType_(mesh_.boundary().size(), NONE),
    wallH_(mesh_.boundary().size(), Zero),
    wallTinf_(mesh_.boundary().size(), Zero),
    wallQ_(mesh_.boundary().size(), Zero),
    cg_(dict_.getOrDefault<scalar>("cGas", 1006)),
    TgVoid_(dict_.getOrDefault<scalar>("TgVoid", voidT_)),
    kGasCell_
    (
        IOobject
        (
            "kGasBed",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimPower/dimLength/dimTemperature, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    nSub_(1),
    nSubAdv_(1),
    variableGas_(dict_.getOrDefault<Switch>("variableGas", false)),
    beta_(dict_.getOrDefault<scalar>("plasticHeating", 0)),
    plasticHeat_(0),
    epsPOld_()
{
    if (rhoP_ <= 0 && gas_)
    {
        rhoP_ = gas_->rhoP();
    }
    if (model_ == "twoTemperature")
    {
        if (!gas_)
        {
            FatalIOErrorInFunction(dict_)
                << "model twoTemperature needs the pore gas "
                << "(constant/gasProperties)" << exit(FatalIOError);
        }
        Tg_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Tg",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }

    // Immersed (STL) walls and bodies (WP T.2b): heat is stored in the open
    // volume, conducted through the open face area, and mixed by mass with
    // the solid's small cut cells; the STL walls are adiabatic
    solid_.addPassenger(T_, limiter_);
    solid_.storeRemapMass();

    if (model_ != "equilibrium" && model_ != "twoTemperature")
    {
        FatalIOErrorInFunction(dict_)
            << "thermalProperties model " << model_ << " unknown; use none, "
            << "equilibrium or twoTemperature" << exit(FatalIOError);
    }
    if (conductivity_ == "constant")
    {
        if (!dict_.found("k") || k_ < 0)
        {
            FatalIOErrorInFunction(dict_)
                << "conductivity constant needs k >= 0" << exit(FatalIOError);
        }
    }
    else if (conductivity_ == "zehnerSchlunder")
    {
        if (kGas_ <= 0 || kParticle_ <= 0 || rhoP_ <= 0)
        {
            FatalIOErrorInFunction(dict_)
                << "conductivity zehnerSchlunder needs kParticle, rhoP "
                << "(and kGas, default 0.02569) positive" << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(dict_)
            << "thermalProperties conductivity " << conductivity_
            << " unknown; use constant or zehnerSchlunder"
            << exit(FatalIOError);
    }
    if (c_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "c must be positive" << exit(FatalIOError);
    }

    // Wall heat transfer (WP T.3): walls { patch { h; Tinf; } or { q; } }
    if (const dictionary* wd = dict_.findDict("walls"))
    {
        for (const entry& e : *wd)
        {
            const label patchi = mesh_.boundaryMesh().findPatchID(e.keyword());
            if (patchi < 0 || !e.isDict())
            {
                FatalIOErrorInFunction(*wd)
                    << "walls entry " << e.keyword() << " is not a patch "
                    << "sub-dictionary" << exit(FatalIOError);
            }
            const dictionary& d = e.dict();
            if (d.found("q"))
            {
                wallType_[patchi] = FLUX;
                wallQ_[patchi] = d.get<scalar>("q");
            }
            else
            {
                wallType_[patchi] = CONVECTIVE;
                wallH_[patchi] = d.get<scalar>("h");
                wallTinf_[patchi] = d.get<scalar>("Tinf");
                if (wallH_[patchi] <= 0)
                {
                    FatalIOErrorInFunction(d)
                        << "h must be positive" << exit(FatalIOError);
                }
            }
        }
    }
    if (variableGas_ && !gas_)
    {
        FatalIOErrorInFunction(dict_)
            << "variableGas needs the pore gas (constant/gasProperties)"
            << exit(FatalIOError);
    }
    if (beta_ < 0 || beta_ > 1)
    {
        FatalIOErrorInFunction(dict_)
            << "plasticHeating must be within 0 to 1" << exit(FatalIOError);
    }
    // WP T.7: the gas at its own temperature from the start (the starting
    // p is the gas pressure at that temperature: nothing is rescaled)
    if (variableGas_)
    {
        gas_->setGasTemperature(gasTemperature(), false);
    }
    correctK();

    energy0_ = energy();

    const label oldPrecision = Info().precision(16);
    Info<< "Temperature (WP T.1-T.3): model " << model_ << ", c = " << c_
        << " J/(kg K), conductivity " << conductivity_;
    if (conductivity_ == "constant")
    {
        Info<< " k = " << k_ << " W/(m K)";
    }
    else
    {
        Info<< " kGas = " << kGas_ << " kParticle = " << kParticle_
            << " W/(m K), rhoP = " << rhoP_ << " kg/m^3";
    }
    Info<< nl
        << "Thermal energy at start = " << energy0_ << " J" << nl
        << "Conduction step limit   = " << stableDeltaT() << " s" << nl;
    if (gas_)
    {
        Info<< "Gas heat (WP T.4): " << (twoT() ? "own temperature Tg, "
            "Gunn grain-gas exchange" : "in equilibrium with the grains")
            << ", cGas = " << cg_ << " J/(kg K)" << nl;
        Info<< "Gas properties (WP T.7): "
            << (
                   variableGas_
                 ? "ideal gas at the gas temperature (variableGas)"
                 : "at the case temperature T of gasProperties"
               )
            << nl;
    }
    if (beta_ > 0)
    {
        Info<< "Plastic heating (WP T.7): beta = " << beta_ << nl;
    }
    Info<< endl;
    Info().precision(oldPrecision);
}


Foam::autoPtr<Foam::heatTransfer> Foam::heatTransfer::New
(
    eulerianSolid& solid,
    poreGas* gas
)
{
    IOobject io(thermalPropertiesIO(solid.mesh()));
    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return nullptr;
    }

    const word model = IOdictionary(io).getOrDefault<word>("model", "none");
    if (model == "none")
    {
        Info<< "Temperature: model none (as without thermalProperties)" << nl
            << endl;
        return nullptr;
    }
    return autoPtr<heatTransfer>::New(solid, io, gas);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::heatTransfer::fill() const
{
    tmp<scalarField> tf
    (
        solid_.freeSurface()
      ? tmp<scalarField>::New
        (
            max(min(solid_.alpha().primitiveField(), scalar(1)), scalar(0))
        )
      : tmp<scalarField>::New(mesh_.nCells(), scalar(1))
    );
    if (solid_.immersed())
    {
        // Covered cells (inside a body) hold no material (WP T.2b)
        const tmp<scalarField> tOpen(solid_.cellOpen());
        scalarField& f = tf.ref();
        forAll(f, celli)
        {
            if (tOpen()[celli] <= 0)
            {
                f[celli] = 0;
            }
        }
    }
    return tf;
}


Foam::tmp<Foam::surfaceScalarField> Foam::heatTransfer::openArea() const
{
    // Face area open to the material: |S_f| times the open fraction with
    // immersed walls (0 where either side is covered), else |S_f| (WP T.2b)
    tmp<surfaceScalarField> tA(new surfaceScalarField("openArea", mesh_.magSf()));
    if (const surfaceScalarField* fo = solid_.faceOpenPtr())
    {
        tA.ref() *= *fo;
    }
    return tA;
}


void Foam::heatTransfer::correctK()
{
    scalarField& k = kCell_.primitiveFieldRef();
    if (conductivity_ == "constant")
    {
        k = k_;
    }
    else
    {
        const scalarField& rho = solid_.rho().primitiveField();
        const tmp<scalarField> tA(fill());
        const scalarField& a = tA();
        const tmp<scalarField> tkG(kGasField());     // WP T.7
        const scalarField& kG = tkG();
        forAll(k, celli)
        {
            if (rho[celli] <= 0 || a[celli] <= 0)
            {
                k[celli] = 0;
                continue;
            }
            // packing of the filled part, kept inside (0, 0.9]
            const scalar phiS =
                min(max(rho[celli]/(a[celli]*rhoP_), scalar(1e-3)), 0.9);
            if (twoT())
            {
                scalar kg, ks;
                zehnerSchlunderParts(1 - phiS, kG[celli], kParticle_, kg, ks);
                k[celli] = ks;
            }
            else
            {
                k[celli] = zehnerSchlunder(1 - phiS, kG[celli], kParticle_);
            }
        }
    }
    kCell_.correctBoundaryConditions();

    // Gas part e k_g,eff = (1 - sqrt(phi_s)) kGas with the cell's own solids
    // fraction (void: kGas) (WP T.4, twoTemperature)
    if (twoT())
    {
        const scalarField& phiS = gas_->phiS().primitiveField();
        const tmp<scalarField> tkG(kGasField());     // WP T.7
        const scalarField& kG = tkG();
        scalarField& kg = kGasCell_.primitiveFieldRef();
        forAll(kg, celli)
        {
            kg[celli] =
                (1 - Foam::sqrt(min(max(phiS[celli], scalar(0)), 1.0)))
               *kG[celli];
        }
        kGasCell_.correctBoundaryConditions();
    }
}


Foam::tmp<Foam::scalarField> Foam::heatTransfer::kGasField() const
{
    // WP T.7: kGas at the case temperature (gasProperties T) scaled by
    // Sutherland's law for the conductivity of air, S = 194 K (White,
    // Viscous Fluid Flow, Table 1-3), with variableGas
    if (!variableGas_ || !gas_)
    {
        return tmp<scalarField>::New(mesh_.nCells(), kGas_);
    }
    const scalar S = 194;
    auto law = [S](const scalar T) { return Foam::pow(T, 1.5)/(T + S); };
    const scalar ref = law(gas_->T());
    const scalarField& Tgas = gasTemperature().primitiveField();
    tmp<scalarField> tk(new scalarField(Tgas.size()));
    scalarField& k = tk.ref();
    forAll(k, celli)
    {
        k[celli] = kGas_*law(Tgas[celli])/ref;
    }
    return tk;
}


Foam::tmp<Foam::volScalarField> Foam::heatTransfer::fillField
(
    const bool all
) const
{
    tmp<volScalarField> tf
    (
        new volScalarField
        (
            IOobject
            (
                "thermalFill",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    if (all)
    {
        tf.ref().primitiveFieldRef() = 1;
    }
    else
    {
        tf.ref().primitiveFieldRef() = fill();
    }
    tf.ref().correctBoundaryConditions();
    return tf;
}


Foam::tmp<Foam::scalarField> Foam::heatTransfer::capacity() const
{
    return tmp<scalarField>::New
    (
        solid_.rho().primitiveField()*c_*solid_.openVolume()
    );
}


Foam::tmp<Foam::scalarField> Foam::heatTransfer::gasCapacity() const
{
    if (!gas_)
    {
        return tmp<scalarField>::New(mesh_.nCells(), Zero);
    }
    // R T per cell (WP T.7)
    const scalarField& RT = gas_->RT().primitiveField();
    return tmp<scalarField>::New
    (
        gas_->phiG().primitiveField()
       *(gas_->p0() + gas_->p().primitiveField())/RT
       *cg_*solid_.openVolume()
    );
}


Foam::scalar Foam::heatTransfer::energy() const
{
    const tmp<scalarField> tCg(gasCapacity());
    if (twoT())
    {
        return gSum(capacity()*T_.primitiveField())
             + gSum(tCg()*Tg_().primitiveField());
    }
    return gSum((capacity() + tCg())*T_.primitiveField());
}


Foam::scalar Foam::heatTransfer::stableDeltaT
(
    const volScalarField& kF,
    const volScalarField& aVol,
    const scalarField& C
) const
{
    // Explicit limit per cell: C_P / sum_f k_f |S_f| delta_f
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const tmp<surfaceScalarField> tOpenArea(openArea());
    const surfaceScalarField& magSf = tOpenArea();
    const surfaceScalarField& ndc = mesh_.nonOrthDeltaCoeffs();
    const scalarField& a = aVol.primitiveField();
    const scalarField& kc = kF.primitiveField();
    scalarField sumG(mesh_.nCells(), Zero);
    forAll(nei, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar g =
            faceK(kc[P], kc[N])*min(a[P], a[N])*magSf[facei]*ndc[facei];
        sumG[own[facei]] += g;
        sumG[nei[facei]] += g;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (isA<emptyFvPatch>(p))
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const scalarField& pMag = magSf.boundaryField()[patchi];
        const scalarField& pNdc = ndc.boundaryField()[patchi];
        const scalarField aN
        (
            p.coupled()
          ? aVol.boundaryField()[patchi].patchNeighbourField()
          : tmp<scalarField>::New(fc.size(), scalar(1))
        );
        const scalarField kN
        (
            p.coupled()
          ? kF.boundaryField()[patchi].patchNeighbourField()
          : tmp<scalarField>::New(scalarField(kc, fc))
        );
        forAll(fc, i)
        {
            sumG[fc[i]] +=
                faceK(kc[fc[i]], kN[i])*min(a[fc[i]], aN[i])*pMag[i]*pNdc[i];
        }
    }

    scalar dtMin = GREAT;
    forAll(C, celli)
    {
        if (sumG[celli] > VSMALL && C[celli] > VSMALL)
        {
            dtMin = min(dtMin, C[celli]/sumG[celli]);
        }
    }
    return returnReduce(dtMin, minOp<scalar>());
}


Foam::scalar Foam::heatTransfer::conduct
(
    volScalarField& Tfield,
    const volScalarField& kF,
    const volScalarField& aVol,
    const scalarField& C,
    const bool walls,
    const scalar dt
)
{
    Tfield.correctBoundaryConditions();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const tmp<surfaceScalarField> tOpenArea(openArea());
    const surfaceScalarField& magSf = tOpenArea();
    const scalarField& w = mesh_.weights();
    const surfaceScalarField& ndc = mesh_.nonOrthDeltaCoeffs();
    const surfaceVectorField& corr = mesh_.nonOrthCorrectionVectors();
    const scalarField& V = mesh_.V();
    const scalarField& a = aVol.primitiveField();
    const scalarField& kc = kF.primitiveField();
    scalarField& T = Tfield.primitiveFieldRef();

    // Gauss gradient with linear face values
    volVectorField gradT
    (
        IOobject
        (
            "gradT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimTemperature/dimLength, Zero),
        fvPatchFieldBase::calculatedType()
    );
    vectorField& g = gradT.primitiveFieldRef();
    forAll(nei, facei)
    {
        // a face to an empty cell takes the material side's value
        const label P = own[facei];
        const label N = nei[facei];
        const scalar wf =
            a[P] > 0 && a[N] > 0 ? w[facei] : (a[P] > 0 ? 1 : 0);
        const scalar Tf = wf*T[P] + (1 - wf)*T[N];
        g[P] += Sf[facei]*Tf;
        g[N] -= Sf[facei]*Tf;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const fvPatchScalarField& pT = Tfield.boundaryField()[patchi];
        if (pT.coupled())
        {
            // a processor patch holds the neighbour cell values
            const scalarField TN(pT.patchNeighbourField());
            const scalarField aN
            (
                aVol.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pw = p.weights();
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar wf =
                    a[P] > 0 && aN[i] > 0 ? pw[i] : (a[P] > 0 ? 1 : 0);
                g[P] += pSf[i]*(wf*T[P] + (1 - wf)*TN[i]);
            }
        }
        else
        {
            forAll(fc, i)
            {
                g[fc[i]] += pSf[i]*pT[i];
            }
        }
    }
    g /= V;
    gradT.correctBoundaryConditions();

    // Heat flow into each cell over dt
    scalarField dE(mesh_.nCells(), Zero);
    const scalarField& ndcI = ndc.primitiveField();
    const vectorField& corrI = corr.primitiveField();
    forAll(nei, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar kf = faceK(kc[P], kc[N])*min(a[P], a[N]);
        if (kf <= 0)
        {
            continue;
        }
        const vector gf = w[facei]*g[P] + (1 - w[facei])*g[N];
        const scalar snGrad =
            (T[N] - T[P])*ndcI[facei] + (corrI[facei] & gf);
        const scalar q = dt*kf*magSf[facei]*snGrad;    // into P
        dE[P] += q;
        dE[N] -= q;
    }

    scalar inflow = 0;
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (isA<emptyFvPatch>(p))
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const fvPatchScalarField& pT = Tfield.boundaryField()[patchi];
        const scalarField& pMag = magSf.boundaryField()[patchi];
        const scalarField& pNdc = ndc.boundaryField()[patchi];

        if (pT.coupled())
        {
            const scalarField TN(pT.patchNeighbourField());
            const scalarField aN(aVol.boundaryField()[patchi].patchNeighbourField());
            const scalarField kN
            (
                kF.boundaryField()[patchi].patchNeighbourField()
            );
            const vectorField gN
            (
                gradT.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pw = p.weights();
            const vectorField& pCorr = corr.boundaryField()[patchi];
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar kf = faceK(kc[P], kN[i])*min(a[P], aN[i]);
                if (kf <= 0)
                {
                    continue;
                }
                const vector gf = pw[i]*g[P] + (1 - pw[i])*gN[i];
                const scalar snGrad =
                    (TN[i] - T[P])*pNdc[i] + (pCorr[i] & gf);
                dE[P] += dt*kf*pMag[i]*snGrad;
            }
        }
        else if (walls && wallType_[patchi] == FLUX)
        {
            // walls { patch { q; } }: heat flux into the material [W/m^2]
            forAll(fc, i)
            {
                const scalar q = dt*wallQ_[patchi]*a[fc[i]]*pMag[i];
                dE[fc[i]] += q;
                inflow += q;
            }
        }
        else if (walls && wallType_[patchi] == CONVECTIVE)
        {
            // walls { patch { h; Tinf; } }: the film and the half cell in
            // series, q = (Tinf - T_P)/(1/h + 1/(k delta)) (WP T.3)
            const scalar h = wallH_[patchi];
            const scalar Tinf = wallTinf_[patchi];
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar kd = kc[P]*pNdc[i];
                if (kd <= 0 || a[P] <= 0)
                {
                    continue;
                }
                const scalar U = 1/(1/h + 1/kd);
                const scalar q = dt*U*a[P]*pMag[i]*(Tinf - T[P]);
                dE[P] += q;
                inflow += q;
            }
        }
        else if (pT.fixesValue())
        {
            // Dirichlet: (T_b - T_P) delta + (n - delta d) . grad T_P
            // Cf - C (fvPatch::delta() is the normal part only for
            // non-coupled patches)
            const vectorField d(p.Cf() - p.Cn());
            const vectorField nf(p.nf());
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar kf = kc[P]*a[P];
                const scalar dl = 1/max(nf[i] & d[i], 0.05*mag(d[i]));
                const scalar snGrad =
                    (pT[i] - T[P])*dl + ((nf[i] - dl*d[i]) & g[P]);
                const scalar q = dt*kf*pMag[i]*snGrad;
                dE[P] += q;
                inflow += q;
            }
        }
        else
        {
            // zeroGradient, fixedGradient (heat flux), symmetry, wedge
            const scalarField sn(pT.snGrad());
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar q = dt*kc[P]*a[P]*pMag[i]*sn[i];
                dE[P] += q;
                inflow += q;
            }
        }
    }

    forAll(T, celli)
    {
        if (C[celli] > VSMALL)
        {
            T[celli] += dE[celli]/C[celli];
        }
    }
    Tfield.correctBoundaryConditions();

    return returnReduce(inflow, sumOp<scalar>());
}


void Foam::heatTransfer::beforeLagrange()
{
    if (beta_ > 0)
    {
        epsPOld_ = solid_.epsP().primitiveField();
        // The heat capacity and material volume at the start of the step:
        // the Lagrange step keeps each cell's mass but sets rho = m/V_L,
        // so capacity() before the remap would not be m c (WP T.7)
        CPre_ = capacity();
        if (gas_ && !twoT())
        {
            CPre_ += gasCapacity();
        }
        VMat_ = fill()*solid_.openVolume();
    }
}


void Foam::heatTransfer::beforeRemap()
{
    // WP T.7: plastic heating. The Lagrange step's plastic work per unit
    // material volume sigma : d(epsP) (stress after the update), times
    // beta, heats the material where it is (fill times the open volume at
    // the start of the step is the material's volume, the capacity is that
    // of the step's start, m c); with equilibrium the gas shares it
    if (beta_ > 0 && epsPOld_.size() == mesh_.nCells())
    {
        const symmTensorField& sig = solid_.sigma().primitiveField();
        const symmTensorField& epsP = solid_.epsP().primitiveField();
        const scalarField& C = CPre_;
        scalarField& T = T_.primitiveFieldRef();
        scalar added = 0;
        forAll(T, celli)
        {
            if (C[celli] <= VSMALL || VMat_[celli] <= 0)
            {
                continue;
            }
            const scalar dE =
                beta_*(sig[celli] && (epsP[celli] - epsPOld_[celli]))
               *VMat_[celli];
            T[celli] += dE/C[celli];
            added += dE;
        }
        plasticHeat_ += returnReduce(added, sumOp<scalar>());
        T_.correctBoundaryConditions();
    }
    TOld_ = T_.primitiveField();
}


void Foam::heatTransfer::afterRemap()
{
    const surfaceScalarField& dm = solid_.remapMass();
    scalar in = 0;
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (p.coupled() || isA<emptyFvPatch>(p))
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const scalarField& pdm = dm.boundaryField()[patchi];
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        forAll(fc, i)
        {
            if (pdm[i] > 0)
            {
                in -= c_*pdm[i]*TOld_[fc[i]];       // out with the cell's T
            }
            else if (pdm[i] < 0)
            {
                in -= c_*pdm[i]*pT[i];              // in at the inlet's T
            }
        }
    }
    carriedIn_ += returnReduce(in, sumOp<scalar>());

    // Empty cells hold no heat; write voidT there
    if (solid_.freeSurface())
    {
        const scalarField& rho = solid_.rho().primitiveField();
        scalarField& T = T_.primitiveFieldRef();
        forAll(T, celli)
        {
            if (rho[celli] <= 0)
            {
                T[celli] = voidT_;
            }
        }
    }
    T_.correctBoundaryConditions();

    // WP T.7: the gas in each cell at its moved temperature (equilibrium)
    if (variableGas_ && !twoT())
    {
        gas_->setGasTemperature(T_, true);
    }
}


Foam::label Foam::heatTransfer::nSubCycles
(
    const scalar dt,
    const scalar dtc
)
{
    // Sub-cycle when the step exceeds half the explicit limit
    return dt > 0.5*dtc ? label(std::ceil(dt/(0.5*dtc))) : label(1);
}


Foam::scalar Foam::heatTransfer::stableDeltaT() const
{
    const tmp<volScalarField> tAVol(fillField());
    const tmp<scalarField> tCs(capacity());
    if (!gas_)
    {
        return stableDeltaT(kCell_, tAVol(), tCs());
    }
    const tmp<scalarField> tCg(gasCapacity());
    if (!twoT())
    {
        return stableDeltaT(kCell_, tAVol(), tCs() + tCg());
    }
    const tmp<volScalarField> tOne(fillField(true));
    return min
    (
        stableDeltaT(kCell_, tAVol(), tCs()),
        stableDeltaT(kGasCell_, tOne(), tCg())
    );
}


void Foam::heatTransfer::totalGasFlux(surfaceScalarField& F, const scalar dt)
    const
{
    // Gas mass flux out of the owner [kg/s]: relative to the grains (the
    // gas solve) plus the gas the grains carry, (dt phiL - dm/rhoP) at the
    // donor's density (decision 0128), not across faces open to the void
    // (decision 0129)
    F = gas_->massFlux();
    if (!gas_->coupled() || dt <= 0)
    {
        return;
    }
    // R T per cell, the donor's (WP T.7)
    const scalarField& RT = gas_->RT().primitiveField();
    const scalar p0 = gas_->p0();
    const scalar rhoP = gas_->rhoP();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& p = gas_->p().primitiveField();
    const scalarField& phiL = solid_.phiL().primitiveField();
    const scalarField& dm = solid_.remapMass().primitiveField();
    const bool fs = solid_.freeSurface();
    const scalarField& vd = gas_->voidCells().primitiveField();
    scalarField& Fi = F.primitiveFieldRef();
    forAll(own, facei)
    {
        const bool vP = fs && vd[own[facei]] > 0.5;
        const bool vN = fs && vd[nei[facei]] > 0.5;
        const scalar G =
            (gas_->flow() ? (vP || vN) : vP != vN)
          ? -dm[facei]/rhoP
          : dt*phiL[facei] - dm[facei]/rhoP;
        const scalar pd = G > 0 ? p[own[facei]] : p[nei[facei]];
        const scalar RTd = G > 0 ? RT[own[facei]] : RT[nei[facei]];
        Fi[facei] += G*(p0 + pd)/RTd/dt;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& pp = mesh_.boundary()[patchi];
        if (!pp.coupled())
        {
            continue;
        }
        const labelUList& fc = pp.faceCells();
        const scalarField& pPhi = solid_.phiL().boundaryField()[patchi];
        const scalarField& pdm = solid_.remapMass().boundaryField()[patchi];
        const scalarField pN(gas_->p().boundaryField()[patchi].patchNeighbourField());
        const scalarField RTN(gas_->RT().boundaryField()[patchi].patchNeighbourField());
        const scalarField vN
        (
            gas_->voidCells().boundaryField()[patchi].patchNeighbourField()
        );
        scalarField& Fp = F.boundaryFieldRef()[patchi];
        forAll(fc, i)
        {
            const bool vP = fs && vd[fc[i]] > 0.5;
            const bool vB = fs && vN[i] > 0.5;
            const scalar G =
                (gas_->flow() ? (vP || vB) : vP != vB)
              ? -pdm[i]/rhoP
              : dt*pPhi[i] - pdm[i]/rhoP;
            const scalar pd = G > 0 ? p[fc[i]] : pN[i];
            const scalar RTd = G > 0 ? RT[fc[i]] : RTN[i];
            Fp[i] += G*(p0 + pd)/RTd/dt;
        }
    }
}


Foam::scalar Foam::heatTransfer::advect
(
    volScalarField& Tfield,
    const surfaceScalarField& F,
    const scalarField& C,
    const scalar dt
)
{
    // Upwind, advective form: C_P dT_P = dt c_g sum_in |F| (T_up - T_P),
    // so a uniform temperature stays uniform under any flow; the heat that
    // crosses non-coupled boundaries is returned (in at the patch value,
    // out at the cell's). Sub-cycled to inflow Courant 0.5.
    // C is the capacity at the end of the step. The gas mass, and so C,
    // changed over the step by c_g dt (net inflow); a sub-step k of n ends
    // at C - c_g h (n - 1 - k) (net inflow), which keeps the energy exactly
    // when it is sub-cycled (WP T.7; one sub-step: C itself)
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& Fi = F.primitiveField();

    scalarField sumIn(mesh_.nCells(), Zero);
    scalarField netIn(mesh_.nCells(), Zero);
    forAll(own, facei)
    {
        if (Fi[facei] > 0) sumIn[nei[facei]] += Fi[facei];
        else sumIn[own[facei]] -= Fi[facei];
        netIn[nei[facei]] += Fi[facei];
        netIn[own[facei]] -= Fi[facei];
    }
    forAll(mesh_.boundary(), patchi)
    {
        const scalarField& Fp = F.boundaryField()[patchi];
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        if (isA<emptyFvPatch>(mesh_.boundary()[patchi])) continue;
        forAll(fc, i)
        {
            if (Fp[i] < 0) sumIn[fc[i]] -= Fp[i];
            netIn[fc[i]] -= Fp[i];
        }
    }
    scalar rate = 0;
    forAll(C, celli)
    {
        if (C[celli] > VSMALL)
        {
            rate = max(rate, cg_*sumIn[celli]/C[celli]);
        }
    }
    reduce(rate, maxOp<scalar>());
    nSubAdv_ = max(label(1), label(std::ceil(dt*rate/0.5)));
    const scalar h = dt/nSubAdv_;

    scalar in = 0;
    for (label n = 0; n < nSubAdv_; ++n)
    {
        Tfield.correctBoundaryConditions();
        const scalarField& T = Tfield.primitiveField();
        scalarField dE(mesh_.nCells(), Zero);
        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            if (Fi[facei] > 0)
            {
                dE[N] += h*cg_*Fi[facei]*(T[P] - T[N]);
            }
            else if (Fi[facei] < 0)
            {
                dE[P] -= h*cg_*Fi[facei]*(T[N] - T[P]);
            }
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& pp = mesh_.boundary()[patchi];
            if (isA<emptyFvPatch>(pp)) continue;
            const labelUList& fc = pp.faceCells();
            const scalarField& Fp = F.boundaryField()[patchi];
            const fvPatchScalarField& pT = Tfield.boundaryField()[patchi];
            const scalarField TN
            (
                pp.coupled() ? pT.patchNeighbourField() : tmp<scalarField>(pT)
            );
            forAll(fc, i)
            {
                const label P = fc[i];
                if (Fp[i] < 0)
                {
                    dE[P] -= h*cg_*Fp[i]*(TN[i] - T[P]);
                    if (!pp.coupled())
                    {
                        in -= h*cg_*Fp[i]*TN[i];
                    }
                }
                else if (Fp[i] > 0 && !pp.coupled())
                {
                    in -= h*cg_*Fp[i]*T[P];
                }
            }
        }
        scalarField& Tw = Tfield.primitiveFieldRef();
        const label left = nSubAdv_ - 1 - n;
        forAll(Tw, celli)
        {
            if (C[celli] > VSMALL)
            {
                const scalar Ck =
                    left > 0 ? C[celli] - cg_*h*left*netIn[celli] : C[celli];
                Tw[celli] += dE[celli]/(Ck > VSMALL ? Ck : C[celli]);
            }
        }
    }
    Tfield.correctBoundaryConditions();
    return returnReduce(in, sumOp<scalar>());
}


void Foam::heatTransfer::exchange
(
    const scalar dt,
    const scalarField& Cs,
    const scalarField& Cg
)
{
    // Grain-gas heat transfer, Gunn (1978), implicit and exact per cell:
    // the difference decays as exp(-H (1/Cs + 1/Cg) dt), Cs Ts + Cg Tg kept
    const scalarField& rho = solid_.rho().primitiveField();
    const tmp<scalarField> tA(fill());
    const scalarField& a = tA();
    const vectorField& q = gas_->qGas().primitiveField();
    const tmp<scalarField> tV(solid_.openVolume());
    const scalarField& V = tV();
    const scalar d = gas_->d();
    // Viscosity, R T and conductivity per cell (WP T.7)
    const scalarField& mu = gas_->muField().primitiveField();
    const scalarField& RT = gas_->RT().primitiveField();
    const tmp<scalarField> tkG(kGasField());
    const scalarField& kG = tkG();
    const scalarField& p = gas_->p().primitiveField();
    scalarField& Ts = T_.primitiveFieldRef();
    scalarField& Tg = Tg_->primitiveFieldRef();
    forAll(Ts, celli)
    {
        if (Cs[celli] <= VSMALL || Cg[celli] <= VSMALL || a[celli] <= 0)
        {
            continue;
        }
        const scalar phiSm =
            min(max(rho[celli]/(a[celli]*rhoP_), scalar(1e-6)), 0.9);
        const scalar e = 1 - phiSm;
        const scalar rhoG = (gas_->p0() + p[celli])/RT[celli];
        const scalar Re = rhoG*mag(q[celli])*d/mu[celli];
        const scalar Pr = mu[celli]*cg_/kG[celli];
        const scalar Nu =
            (7 - 10*e + 5*e*e)*(1 + 0.7*Foam::pow(Re, 0.2)*Foam::cbrt(Pr))
          + (1.33 - 2.4*e + 1.2*e*e)*Foam::pow(Re, 0.7)*Foam::cbrt(Pr);
        const scalar hAs = Nu*kG[celli]/d*6*(rho[celli]/rhoP_)/d;  // W/(m^3 K)
        const scalar H = hAs*V[celli];
        const scalar mean = (Cs[celli]*Ts[celli] + Cg[celli]*Tg[celli])
                           /(Cs[celli] + Cg[celli]);
        const scalar f =
            Foam::exp(-H*(1/Cs[celli] + 1/Cg[celli])*dt);
        const scalar diff = (Ts[celli] - Tg[celli])*f;
        Ts[celli] = mean + Cg[celli]/(Cs[celli] + Cg[celli])*diff;
        Tg[celli] = mean - Cs[celli]/(Cs[celli] + Cg[celli])*diff;
    }
    T_.correctBoundaryConditions();
    Tg_->correctBoundaryConditions();
}


void Foam::heatTransfer::solve(const scalar dt)
{
    // Capacities, fill fractions and the limit are computed once per step
    correctK();
    const tmp<volScalarField> tAVol(fillField());
    const tmp<scalarField> tCs(capacity());
    if (!gas_)
    {
        nSub_ = nSubCycles(dt, stableDeltaT(kCell_, tAVol(), tCs()));
        const scalar h = dt/nSub_;
        for (label n = 0; n < nSub_; ++n)
        {
            heatIn_ += conduct(T_, kCell_, tAVol(), tCs(), true, h);
        }
        return;
    }

    // With gas (WP T.4)
    const tmp<scalarField> tCg(gasCapacity());
    surfaceScalarField F
    (
        IOobject
        (
            "gasHeatFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    );
    totalGasFlux(F, dt);

    if (!twoT())
    {
        const scalarField C(tCs() + tCg());
        // Carried first (WP T.7): the gas moved in the gas solve at the
        // temperatures of the start of the step, which conduction would
        // change before the advection (an energy error of the order of the
        // flux times the conduction's change, decision 0136)
        carriedIn_ += advect(T_, F, C, dt);
        nSub_ = nSubCycles(dt, stableDeltaT(kCell_, tAVol(), C));
        const scalar h = dt/nSub_;
        for (label n = 0; n < nSub_; ++n)
        {
            heatIn_ += conduct(T_, kCell_, tAVol(), C, true, h);
        }
        if (variableGas_)
        {
            gas_->setGasTemperature(T_, true);      // WP T.7
        }
        return;
    }

    const tmp<volScalarField> tOne(fillField(true));
    nSub_ = nSubCycles
    (
        dt,
        min
        (
            stableDeltaT(kCell_, tAVol(), tCs()),
            stableDeltaT(kGasCell_, tOne(), tCg())
        )
    );
    const scalar h = dt/nSub_;
    // Carried first, as with equilibrium (WP T.7)
    carriedIn_ += advect(Tg_(), F, tCg(), dt);
    for (label n = 0; n < nSub_; ++n)
    {
        heatIn_ += conduct(T_, kCell_, tAVol(), tCs(), true, h);
        heatIn_ += conduct(Tg_(), kGasCell_, tOne(), tCg(), false, h);
    }
    exchange(dt, tCs(), tCg());

    // Void cells of a poreGas (not gasFlow) case are the atmosphere
    if (solid_.freeSurface() && !gas_->flow())
    {
        const scalarField& vd = gas_->voidCells().primitiveField();
        scalarField& Tg = Tg_->primitiveFieldRef();
        forAll(Tg, celli)
        {
            if (vd[celli] > 0.5)
            {
                Tg[celli] = TgVoid_;
            }
        }
        Tg_->correctBoundaryConditions();
    }

    // WP T.7: each cell's gas pressure follows its new temperature
    if (variableGas_)
    {
        gas_->setGasTemperature(Tg_(), true);
    }
}


void Foam::heatTransfer::report() const
{
    const scalar E = energy();
    const label oldPrecision = Info().precision(16);
    Info<< "Thermal energy  = " << E << nl
        << "Heat in         = " << heatIn_ << nl
        << "Heat carried in = " << carriedIn_ << nl;
    if (beta_ > 0)
    {
        Info<< "Plastic heat    = " << plasticHeat_ << nl
            << "Energy balance  = "
            << E - energy0_ - heatIn_ - carriedIn_ - plasticHeat_
            << " (energy - start - heat in - carried in - plastic heat)"
            << endl;
    }
    else
    {
        Info<< "Energy balance  = " << E - energy0_ - heatIn_ - carriedIn_
            << " (energy - start - heat in - carried in)" << endl;
    }
    Info().precision(oldPrecision);
}


// ************************************************************************* //
