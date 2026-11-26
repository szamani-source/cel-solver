/*---------------------------------------------------------------------------*\
    Temperature in the granular material (WP T.1, decision 0135).
    See heatTransfer.H.
\*---------------------------------------------------------------------------*/

#include "heatTransfer.H"
#include "fvMesh.H"
#include "Time.H"
#include "emptyFvPatch.H"

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
    const IOobject& io
)
:
    solid_(solid),
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
    limiter_(dict_.getOrDefault<word>("limiter", "superbee")),
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
    nSub_(1)
{
    if (solid_.immersed())
    {
        FatalIOErrorInFunction(dict_)
            << "temperature with immersed (STL) walls or bodies is not yet "
            << "supported (WP T.2b)" << exit(FatalIOError);
    }
    solid_.addPassenger(T_, limiter_);
    solid_.storeRemapMass();

    if (model_ != "equilibrium")
    {
        FatalIOErrorInFunction(dict_)
            << "thermalProperties model " << model_ << " unknown; use none or "
            << "equilibrium (twoTemperature comes with WP T.4)"
            << exit(FatalIOError);
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
        << "Conduction step limit   = " << stableDeltaT() << " s" << nl << endl;
    Info().precision(oldPrecision);
}


Foam::autoPtr<Foam::heatTransfer> Foam::heatTransfer::New
(
    eulerianSolid& solid
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
    return autoPtr<heatTransfer>::New(solid, io);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::heatTransfer::fill() const
{
    if (solid_.freeSurface())
    {
        return tmp<scalarField>::New
        (
            max(min(solid_.alpha().primitiveField(), scalar(1)), scalar(0))
        );
    }
    return tmp<scalarField>::New(mesh_.nCells(), scalar(1));
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
            k[celli] = zehnerSchlunder(1 - phiS, kGas_, kParticle_);
        }
    }
    kCell_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField> Foam::heatTransfer::fillField() const
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
    tf.ref().primitiveFieldRef() = fill();
    tf.ref().correctBoundaryConditions();
    return tf;
}


Foam::tmp<Foam::scalarField> Foam::heatTransfer::capacity() const
{
    return tmp<scalarField>::New
    (
        solid_.rho().primitiveField()*c_*mesh_.V().field()
    );
}


Foam::scalar Foam::heatTransfer::energy() const
{
    return gSum(capacity()*T_.primitiveField());
}


Foam::scalar Foam::heatTransfer::stableDeltaT() const
{
    // Explicit limit per cell: C_P / sum_f k_f |S_f| delta_f
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& ndc = mesh_.nonOrthDeltaCoeffs();
    const tmp<volScalarField> tAVol(fillField());
    const volScalarField& aVol = tAVol();
    const scalarField& a = aVol.primitiveField();

    const scalarField& kc = kCell_.primitiveField();
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
          ? kCell_.boundaryField()[patchi].patchNeighbourField()
          : tmp<scalarField>::New(scalarField(kc, fc))
        );
        forAll(fc, i)
        {
            sumG[fc[i]] +=
                faceK(kc[fc[i]], kN[i])*min(a[fc[i]], aN[i])*pMag[i]*pNdc[i];
        }
    }

    const tmp<scalarField> tC(capacity());
    const scalarField& C = tC();
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


void Foam::heatTransfer::conduct(const scalar dt)
{
    T_.correctBoundaryConditions();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const scalarField& w = mesh_.weights();
    const surfaceScalarField& ndc = mesh_.nonOrthDeltaCoeffs();
    const surfaceVectorField& corr = mesh_.nonOrthCorrectionVectors();
    const scalarField& V = mesh_.V();
    const tmp<volScalarField> tAVol(fillField());
    const volScalarField& aVol = tAVol();
    const scalarField& a = aVol.primitiveField();
    const scalarField& kc = kCell_.primitiveField();
    scalarField& T = T_.primitiveFieldRef();

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
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
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
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        const scalarField& pMag = magSf.boundaryField()[patchi];
        const scalarField& pNdc = ndc.boundaryField()[patchi];

        if (pT.coupled())
        {
            const scalarField TN(pT.patchNeighbourField());
            const scalarField aN(aVol.boundaryField()[patchi].patchNeighbourField());
            const scalarField kN
            (
                kCell_.boundaryField()[patchi].patchNeighbourField()
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
        else if (wallType_[patchi] == FLUX)
        {
            // walls { patch { q; } }: heat flux into the material [W/m^2]
            forAll(fc, i)
            {
                const scalar q = dt*wallQ_[patchi]*a[fc[i]]*pMag[i];
                dE[fc[i]] += q;
                inflow += q;
            }
        }
        else if (wallType_[patchi] == CONVECTIVE)
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

    const tmp<scalarField> tC(capacity());
    const scalarField& C = tC();
    forAll(T, celli)
    {
        if (C[celli] > VSMALL)
        {
            T[celli] += dE[celli]/C[celli];
        }
    }
    T_.correctBoundaryConditions();

    heatIn_ += returnReduce(inflow, sumOp<scalar>());
}


void Foam::heatTransfer::beforeRemap()
{
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
}


void Foam::heatTransfer::solve(const scalar dt)
{
    correctK();
    const scalar dtc = stableDeltaT();
    nSub_ = 1;
    if (dt > 0.5*dtc)
    {
        nSub_ = label(std::ceil(dt/(0.5*dtc)));
    }
    const scalar h = dt/nSub_;
    for (label n = 0; n < nSub_; ++n)
    {
        conduct(h);
    }
}


void Foam::heatTransfer::report() const
{
    const scalar E = energy();
    const label oldPrecision = Info().precision(16);
    Info<< "Thermal energy  = " << E << nl
        << "Heat in         = " << heatIn_ << nl
        << "Heat carried in = " << carriedIn_ << nl
        << "Energy balance  = " << E - energy0_ - heatIn_ - carriedIn_
        << " (energy - start - heat in - carried in)" << endl;
    Info().precision(oldPrecision);
}


// ************************************************************************* //
