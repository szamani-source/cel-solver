/*---------------------------------------------------------------------------*\
    poreGas: gas fractions and permeability (WP F.1, decisions 0125, 0126)
\*---------------------------------------------------------------------------*/

#include "poreGas.H"
#include "Time.H"
#include "fvm.H"
#include "fvc.H"
#include "gaussGrad.H"
#include "laplacianScheme.H"
#include "fixedGradientFvPatchFields.H"
#include "Function1.H"
#include "IStringStream.H"
#include "scalarMatrices.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * Local helpers * * * * * * * * * * * * * * //

namespace
{

Foam::IOobject gasPropertiesIO(const Foam::fvMesh& mesh)
{
    return Foam::IOobject
    (
        "gasProperties",
        mesh.time().constant(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
    );
}


Foam::IOobject fractionIO(const Foam::word& name, const Foam::fvMesh& mesh)
{
    // Computed from rho and alpha, never read (a restart recomputes them)
    return Foam::IOobject
    (
        name,
        mesh.time().timeName(),
        mesh,
        Foam::IOobject::NO_READ,
        Foam::IOobject::AUTO_WRITE
    );
}

Foam::autoPtr<Foam::volScalarField> makePressure
(
    const Foam::eulerianSolid& solid
)
{
    using namespace Foam;
    const fvMesh& mesh = solid.mesh();

    IOobject io
    (
        "p",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if (io.typeHeaderOk<volScalarField>(true))
    {
        Info<< "Pore gas: reading p from " << mesh.time().timeName() << endl;
        return autoPtr<volScalarField>::New(io, mesh);
    }

    // No 0/p: p0 everywhere, every patch sealed but the outlets of the
    // solid, at the ambient pressure (WP F.4)
    io.readOpt(IOobject::NO_READ);
    Info<< "Pore gas: no p file, p = 0 (gauge) and every patch sealed"
        << " but outlets (p = 0)" << endl;
    wordList types
    (
        mesh.boundary().size(), fvPatchFieldBase::zeroGradientType()
    );
    forAll(types, patchi)
    {
        if (solid.isOutlet(mesh.boundary()[patchi].name()))
        {
            types[patchi] = "fixedValue";
        }
    }
    return autoPtr<volScalarField>::New
    (
        io,
        mesh,
        dimensionedScalar(dimPressure, Zero),
        types
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Static functions  * * * * * * * * * * * * * //

Foam::scalar Foam::poreGas::sutherland(const scalar T)
{
    // tools/gas_reference.py: mu0 1.716e-5 Pa s, T0 273.15 K, S 110.4 K
    const scalar mu0 = 1.716e-5, T0 = 273.15, S = 110.4;
    return mu0*Foam::pow(T/T0, 1.5)*(T0 + S)/(T + S);
}


Foam::scalar Foam::poreGas::permeability(const scalar e, const scalar d)
{
    if (e <= 0)
    {
        return 0;
    }
    if (e < eSwitch)
    {
        // Ergun's viscous term
        return e*e*e*d*d/(150*(1 - e)*(1 - e));
    }
    if (e >= 1)
    {
        return GREAT;
    }
    // Wen and Yu with Schiller-Naumann Cd at Re -> 0:
    // K = 18 mu (1-e) e^-2.65 / d^2, k = mu e^2 / K
    return d*d*Foam::pow(e, 4.65)/(18*(1 - e));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreGas::poreGas(const eulerianSolid& solid, const IOobject& io)
:
    solid_(solid),
    dict_(io),
    rhoP_(dict_.get<scalar>("rhoP")),
    d_(dict_.get<scalar>("d")),
    T_(dict_.getOrDefault<scalar>("T", 293.15)),
    p0_(dict_.getOrDefault<scalar>("p0", 101325)),
    drag_(dict_.getOrDefault<word>("drag", "ergunWenYu")),
    dragModel_
    (
        gasDragModel::New(drag_, dict_.optionalSubDict(drag_ + "Coeffs"))
    ),
    rhoG_(p0_/(Rair*T_)),
    mu_(sutherland(T_)),
    phiS_
    (
        fractionIO("phiS", solid.mesh()),
        solid.mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    phiG_
    (
        fractionIO("phiG", solid.mesh()),
        solid.mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    k_
    (
        fractionIO("k", solid.mesh()),
        solid.mesh(),
        dimensionedScalar(dimArea, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    eM_
    (
        IOobject
        (
            "gasFractionFilled",
            solid.mesh().time().timeName(),
            solid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solid.mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    pPtr_(makePressure(solid)),
    qGas_
    (
        fractionIO("qGas", solid.mesh()),
        solid.mesh(),
        dimensionedVector(dimVelocity, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    massFlux_
    (
        IOobject
        (
            "gasMassFlux",
            solid.mesh().time().timeName(),
            solid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solid.mesh(),
        dimensionedScalar(dimMass/dimTime, Zero)
    ),
    inflow_(),
    nOuter_(dict_.getOrDefault<label>("nOuter", 50)),
    nAnderson_(dict_.getOrDefault<label>("nAnderson", 5)),
    pTol_(dict_.getOrDefault<scalar>("pTolerance", 1e-14*p0_)),
    laplacian_
    (
        dict_.getOrDefault<string>("laplacian", "Gauss linear corrected")
    ),
    solverDict_(),
    netIn_(0),
    coupled_(dict_.getOrDefault<Switch>("coupled", true)),
    phiG0_(),
    pOldOld_(),
    dtOld_(-1),
    freeSurface_(solid.freeSurface()),
    alphaVoid_(dict_.getOrDefault<scalar>("alphaVoid", 1e-4)),
    void_
    (
        IOobject
        (
            "gasVoid",
            solid.mesh().time().timeName(),
            solid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solid.mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    pf_
    (
        IOobject
        (
            "pFace",
            solid.mesh().time().timeName(),
            solid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solid.mesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    voidIn_(0),
    flow_(dict_.getOrDefault<word>("model", "poreGas") == "gasFlow"),
    nCorr_(dict_.getOrDefault<label>("nCorr", 1)),
    divUg_(dict_.getOrDefault<string>("divUg", "Gauss upwind")),
    USolverDict_(),
    poreInertia_(flow_ && dict_.getOrDefault<Switch>("poreInertia", true)),
    surfaceCarry_(dict_.getOrDefault<Switch>("surfaceCarry", false)),
    eMin_(dict_.getOrDefault<scalar>("eMin", 0.2)),
    kappa_(dict_.getOrDefault<scalar>("diffusivity", 0))
{
    if (rhoP_ <= 0 || d_ <= 0 || T_ <= 0 || p0_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "rhoP, d, T and p0 must be positive" << exit(FatalIOError);
    }
    if (solid.rho0() > rhoP_)
    {
        WarningInFunction
            << "materialProperties rho " << solid.rho0()
            << " is above the grain density rhoP " << rhoP_
            << ": the packing would exceed 1" << endl;
    }

    const label oldPrecision = Info().precision(8);
    Info<< "Pore gas (WP F.1, F.2): air at T = " << T_ << " K, p0 = " << p0_
        << " Pa: rho_g = " << rhoG_ << " kg/m^3, mu = " << mu_ << " Pa s;"
        << nl << "    grains rhoP = " << rhoP_ << " kg/m^3, d = " << d_
        << " m, drag " << drag_ << "; set packing rho/rhoP = "
        << solid.rho0()/rhoP_ << nl << endl;
    Info().precision(oldPrecision);

    // Linear solver controls for each Newton correction (defaults: PCG/DIC,
    // residual down by 1e-3: with 0.1 (F.2) each outer iteration gained
    // one digit and a step took 8 to 10 of them, each with the face
    // coefficients; 1e-3 needs 3 per step at a few more PCG sweeps, WP F.3)
    solverDict_.add("solver", "PCG");
    solverDict_.add("preconditioner", "DIC");
    solverDict_.add("tolerance", 1e-30);
    solverDict_.add("relTol", 1e-3);
    solverDict_.add("maxIter", 1000);
    solverDict_.add("log", 0);
    if (const dictionary* sd = dict_.findDict("solver"))
    {
        solverDict_.merge(*sd);
    }

    // The gas mass flux comes from the matrix
    const_cast<fvMesh&>(solid.mesh()).setFluxRequired(pPtr_->name());

    // Inflow patches: must be fixedGradient in p
    if (const dictionary* in = dict_.findDict("inflow"))
    {
        const volScalarField& p = *pPtr_;
        for (const entry& e : *in)
        {
            const label patchi =
                solid.mesh().boundaryMesh().findPatchID(e.keyword());
            if
            (
                patchi < 0
             || !isA<fixedGradientFvPatchScalarField>
                (
                    p.boundaryField()[patchi]
                )
            )
            {
                FatalIOErrorInFunction(dict_)
                    << "inflow patch " << e.keyword()
                    << " must exist and be fixedGradient in p"
                    << exit(FatalIOError);
            }
            // A number or, WP F.7, any Function1 of time (a ramp: table)
            inflow_.set(e.keyword(), Function1<scalar>::New(e.keyword(), *in));
            const scalar U0 = inflow_[e.keyword()]->value(0);
            Info<< "Pore gas inflow on " << e.keyword() << ": "
                << U0 << " m/s at p0 at t = 0 ("
                << U0*rhoG_ << " kg/(m^2 s))" << endl;
        }
    }

    // Outlets (WP F.4): at a fixed pressure once open (sealed until then)
    forAll(solid.mesh().boundary(), patchi)
    {
        const word& name = solid.mesh().boundary()[patchi].name();
        if
        (
            solid.isOutlet(name)
         && !pPtr_->boundaryField()[patchi].fixesValue()
        )
        {
            FatalIOErrorInFunction(dict_)
                << "outlet " << name << " must be fixedValue in 0/p "
                << "(the pressure outside, usually 0 gauge)"
                << exit(FatalIOError);
        }
    }
    if (alphaVoid_ <= 0 || alphaVoid_ >= 1)
    {
        FatalIOErrorInFunction(dict_)
            << "alphaVoid must be between 0 and 1" << exit(FatalIOError);
    }

    correctFractions();
    phiG0_ = phiG_.primitiveField();
    if (flow_)
    {
        initFlow();
    }
    if (freeSurface_ && !flow_)
    {
        // Void cells start at the void pressure
        scalarField& pc = pPtr_->primitiveFieldRef();
        forAll(pc, celli)
        {
            if (void_[celli] > 0.5)
            {
                pc[celli] = 0;
            }
        }
        pPtr_->correctBoundaryConditions();
        Info<< "Pore gas with a free surface (WP F.4): void below alpha "
            << alphaVoid_ << " at p0, " << gSum(void_.primitiveField())
            << " void cells" << endl;
    }
    correctFacePressure();
    if (freeSurface_ && surfaceCarry_)
    {
        Info<< "Pore gas (WP F.8): the moving surface carries gas to and "
            << "from the void (surfaceCarry)" << endl;
    }
    if (kappa_ < 0 || (kappa_ > 0 && flow_))
    {
        FatalIOErrorInFunction(dict_)
            << "diffusivity must be positive and needs model poreGas"
            << exit(FatalIOError);
    }
    if (kappa_ > 0)
    {
        Info<< "Pore gas (WP F.8e): pore pressure advected with the grains "
            << "and diffused with diffusivity " << kappa_ << " m^2/s "
            << "(no compression or suction by the grains; gas mass is not "
            << "conserved)" << endl;
    }
    if (eMin_ != 0.2)
    {
        Info<< "Pore gas (WP F.8c): gas fraction of the filled part at "
            << "least eMin = " << eMin_ << endl;
    }

    Info<< "Pore gas coupling (WP F.3): "
        << (coupled_ ? "two-way (effective stress, gas carried by the grains)"
                     : "none (solid frozen for the gas, F.2)")
        << nl << endl;
}


Foam::autoPtr<Foam::poreGas> Foam::poreGas::New(const eulerianSolid& solid)
{
    IOobject io(gasPropertiesIO(solid.mesh()));
    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return nullptr;
    }

    const word model = IOdictionary(io).getOrDefault<word>("model", "poreGas");
    if (model == "none")
    {
        Info<< "Pore gas: model none (as without gasProperties)" << nl << endl;
        return nullptr;
    }
    if (model != "poreGas" && model != "gasFlow")
    {
        FatalErrorInFunction
            << "gasProperties model " << model
            << " unknown; use none, poreGas or gasFlow" << exit(FatalError);
    }
    return autoPtr<poreGas>::New(solid, io);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::poreGas::correctFractions()
{
    const scalarField& rho = solid_.rho().primitiveField();
    const scalarField& alpha = solid_.alpha().primitiveField();
    const tmp<scalarField> tOpen(solid_.cellOpen());
    const scalarField& open = tOpen();

    // Gas fraction of the last solve (WP F.3): the old gas mass
    phiG0_ = phiG_.primitiveField();

    scalarField& phiS = phiS_.primitiveFieldRef();
    scalarField& phiG = phiG_.primitiveFieldRef();
    scalarField& k = k_.primitiveFieldRef();

    forAll(rho, celli)
    {
        if (open[celli] < eulerianSolid::coveredOpen)
        {
            // Inside a body: no grains and no gas
            phiS[celli] = 0;
            phiG[celli] = 0;
            k[celli] = 0;
            eM_[celli] = 0;
            continue;
        }

        // Fractions of the open volume: alpha rho_m / rhoP = rho / rhoP
        phiS[celli] = rho[celli]/rhoP_;
        phiG[celli] = 1 - phiS[celli];

        // Permeability of the filled part (decision 0125 item 4)
        if (alpha[celli] < alphaMin)
        {
            eM_[celli] = 0;
            k[celli] = 0;
        }
        else
        {
            // WP F.8 (decision 0133): at least eMin (default 0.2). The remap can
            // leave a nearly empty cell with more mass than alpha holds
            // (rho/alpha > rhoP); its gas fraction then went negative,
            // its permeability 0, and the gas in it was shut in (the
            // pressure diverged at a thin spreading front)
            eM_[celli] =
                max(1 - rho[celli]/(alpha[celli]*rhoP_), eMin_);
            k[celli] =
                dragModel_->permeability(eM_[celli], d_, mu_, rhoG_);
        }
    }

    // Void cells for the gas (WP F.4): free surface only
    if (freeSurface_)
    {
        scalarField& v = void_.primitiveFieldRef();
        forAll(v, celli)
        {
            v[celli] =
                (open[celli] >= eulerianSolid::coveredOpen
              && alpha[celli] < alphaVoid_) ? 1 : 0;
        }
        void_.correctBoundaryConditions();
    }

    eM_.correctBoundaryConditions();
    phiS_.correctBoundaryConditions();
    phiG_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}


Foam::scalar Foam::poreGas::gasMass() const
{
    const scalarField& p = pPtr_->primitiveField();
    const scalarField& phiG = phiG_.primitiveField();
    const tmp<scalarField> tVol(solid_.openVolume());
    // Reference and gauge parts summed apart: adding p0 cell by cell would
    // round each change in p to the last digit of p0
    const scalarField gasVol(tVol()*phiG);
    return (p0_*gSum(gasVol) + gSum(gasVol*p))/(Rair*T_);
}


void Foam::poreGas::faceCoefficients
(
    const volVectorField& gradP,
    surfaceScalarField& DfS,
    surfaceScalarField& DfT
) const
{
    // rho_g at the face times the harmonic mean of the two cells'
    // mobilities (0 if either side has none), each at the face pressure
    // gradient: the interpolated cell gradient with its normal part from
    // the two cell values (corrected snGrad); times the open face fraction
    const fvMesh& mesh = solid_.mesh();
    const volScalarField& p = *pPtr_;
    const scalar RT = Rair*T_;
    // WP F.8 (decision 0133): the face gas density from at least 1e-3 p0,
    // so a diverging iterate (gauge p below -p0) cannot make it negative
    const scalar pFloor = 1e-3*p0_;
    const surfaceScalarField& w = mesh.weights();
    const surfaceVectorField nf(mesh.Sf()/mesh.magSf());
    const surfaceScalarField& dc = mesh.nonOrthDeltaCoeffs();
    const surfaceVectorField& corr = mesh.nonOrthCorrectionVectors();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& eM = eM_.primitiveField();
    const scalarField& k = k_.primitiveField();
    // WP F.7: with the pore gas inertia, the drag law K + e rhoG/dt and the
    // gradient less the old flux's part, rho_g q_o/(e dt), on the faces
    // between gas cells (and gas and void); plain on non-coupled boundaries
    const gasDragModel* lawPtr =
        inertialDrag_ ? inertialDrag_.get() : dragModel_.get();
    const bool inertia = bool(inertialDrag_);

    // Harmonic mean of the two sides at gradient magnitude G
    auto face =
        [&]
        (
            const scalar eP, const scalar kP,
            const scalar eN, const scalar kN,
            const scalar wf, const scalar rho, const scalar G,
            scalar& dS, scalar& dT
        )
    {
        if (kP <= 0 || kN <= 0)
        {
            dS = 0;
            dT = 0;
            return;
        }
        const gasDragModel& drag = *lawPtr;
        const scalar sP = drag.mobility(eP, G, d_, mu_, rho);
        const scalar sN = drag.mobility(eN, G, d_, mu_, rho);
        const scalar tP = drag.tangentMobility(eP, G, d_, mu_, rho);
        const scalar tN = drag.tangentMobility(eN, G, d_, mu_, rho);
        dS = rho/(wf/sP + (1 - wf)/sN);
        dT = rho/(wf/tP + (1 - wf)/tN);
    };

    // Face between a gas cell and a void cell (WP F.4): the gas cell's
    // mobility over the distance from its centre to the face, fP of the
    // centre distance (a pressure face at the void pressure, 0 gauge)
    // With gas flow in the void (WP F.6) the void side is the void cell's
    // pressure, and the gas density the face's (pV the void cell's p)
    auto faceVoid =
        [&]
        (
            const scalar eP, const scalar kP, const scalar fP,
            const scalar G, scalar& dS, scalar& dT, const scalar pV = 0
        )
    {
        if (kP <= 0 || fP <= 0)
        {
            dS = 0;
            dT = 0;
            return;
        }
        const scalar rho = max(flow_ ? p0_ + pV : p0_, pFloor)/RT;
        const gasDragModel& drag = *lawPtr;
        dS = rho*drag.mobility(eP, G, d_, mu_, rho)/fP;
        dT = rho*drag.tangentMobility(eP, G, d_, mu_, rho)/fP;
    };
    const scalarField& vd = void_.primitiveField();
    const bool fs = freeSurface_;

    // Face gradient magnitude from the interpolated gradient g and the
    // normal derivative gn
    auto mag3 =
        [](const vector& g, const vector& n, const scalar gn)
    {
        return mag(g + n*(gn - (n & g)));
    };

    const vectorField gI_
    (
        inertia ? vectorField(gInertia_->primitiveField()) : vectorField()
    );
    {
        scalarField& DS = DfS.primitiveFieldRef();
        scalarField& DT = DfT.primitiveFieldRef();
        forAll(DS, facei)
        {
            const label P = own[facei], N = nei[facei];
            const scalar wf = w[facei];
            if (fs && (vd[P] > 0.5 || vd[N] > 0.5))
            {
                if (vd[P] > 0.5 && vd[N] > 0.5)
                {
                    // Gas flow inside the void (WP F.6): rho_g 1/A of the
                    // momentum equation; none without (F.4)
                    DS[facei] =
                        flow_
                      ? max(p0_ + wf*p[P] + (1 - wf)*p[N], pFloor)/RT
                       *rAUfPtr_()[facei]
                      : 0;
                    DT[facei] = DS[facei];
                    continue;
                }
                const bool ownGas = vd[P] < 0.5;
                const label C = ownGas ? P : N;
                const scalar fC = ownGas ? 1 - wf : wf;
                // The void side is at 0 without gas flow (F.4)
                const scalar gn = (p[N] - p[P])*dc[facei]/max(fC, SMALL);
                const vector gI = inertia ? gI_[facei] : Zero;
                faceVoid
                (
                    eM[C], k[C], fC,
                    mag3(gradP[C] - gI, nf[facei], gn - (nf[facei] & gI)),
                    DS[facei], DT[facei], ownGas ? p[N] : p[P]
                );
                continue;
            }
            const vector gI = inertia ? gI_[facei] : Zero;
            const vector g = wf*gradP[P] + (1 - wf)*gradP[N] - gI;
            const scalar gn =
                dc[facei]*(p[N] - p[P]) + (corr[facei] & (g + gI))
              - (nf[facei] & gI);
            face
            (
                eM[P], k[P], eM[N], k[N], wf,
                max(p0_ + wf*p[P] + (1 - wf)*p[N], pFloor)/RT,
                mag3(g, nf[facei], gn),
                DS[facei], DT[facei]
            );
        }
    }

    forAll(DfS.boundaryField(), patchi)
    {
        lawPtr = inertia ? inertialDrag_.get() : dragModel_.get();
        const fvPatchScalarField& pB = p.boundaryField()[patchi];
        const fvPatchVectorField& gB = gradP.boundaryField()[patchi];
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField& n = nf.boundaryField()[patchi];
        scalarField& DS = DfS.boundaryFieldRef()[patchi];
        scalarField& DT = DfT.boundaryFieldRef()[patchi];

        if (pB.coupled())
        {
            const scalarField pN(pB.patchNeighbourField());
            const scalarField eN(eM_.boundaryField()[patchi].patchNeighbourField());
            const scalarField kN(k_.boundaryField()[patchi].patchNeighbourField());
            const vectorField gN(gB.patchNeighbourField());
            const scalarField& wf = w.boundaryField()[patchi];
            const scalarField& dcB = dc.boundaryField()[patchi];
            const vectorField& cB = corr.boundaryField()[patchi];
            const scalarField vN
            (
                fs ? void_.boundaryField()[patchi].patchNeighbourField()
                   : tmp<scalarField>::New(DS.size(), Zero)
            );
            const vectorField gIB
            (
                inertia
              ? vectorField(gInertia_->boundaryField()[patchi])
              : vectorField()
            );
            forAll(DS, i)
            {
                const label P = fc[i];
                if (fs && (vd[P] > 0.5 || vN[i] > 0.5))
                {
                    if (vd[P] > 0.5 && vN[i] > 0.5)
                    {
                        DS[i] =
                            flow_
                          ? max(p0_ + wf[i]*p[P] + (1 - wf[i])*pN[i], pFloor)/RT
                           *rAUfPtr_().boundaryField()[patchi][i]
                          : 0;
                        DT[i] = DS[i];
                        continue;
                    }
                    const vector gI = inertia ? gIB[i] : Zero;
                    if (vd[P] > 0.5)
                    {
                        // The gas side is across the cut: the neighbour
                        // processor's face takes its mobility; this side
                        // must give the same coefficient
                        const scalar fC = wf[i];
                        const scalar gn = (pN[i] - p[P])*dcB[i]/max(fC, SMALL);
                        faceVoid
                        (
                            eN[i], kN[i], fC,
                            mag3(gN[i] - gI, n[i], gn - (n[i] & gI)),
                            DS[i], DT[i], p[P]
                        );
                    }
                    else
                    {
                        const scalar fC = 1 - wf[i];
                        const scalar gn = (pN[i] - p[P])*dcB[i]/max(fC, SMALL);
                        faceVoid
                        (
                            eM[P], k[P], fC,
                            mag3(gradP[P] - gI, n[i], gn - (n[i] & gI)),
                            DS[i], DT[i], pN[i]
                        );
                    }
                    continue;
                }
                const vector gI = inertia ? gIB[i] : Zero;
                const vector g = wf[i]*gradP[P] + (1 - wf[i])*gN[i] - gI;
                const scalar gn =
                    dcB[i]*(pN[i] - p[P]) + (cB[i] & (g + gI)) - (n[i] & gI);
                face
                (
                    eM[P], k[P], eN[i], kN[i], wf[i],
                    max(p0_ + wf[i]*p[P] + (1 - wf[i])*pN[i], pFloor)/RT,
                    mag3(g, n[i], gn),
                    DS[i], DT[i]
                );
            }
        }
        else
        {
            // Pore inertia (WP F.7) on fixed-pressure faces; an inflow
            // (fixed flux) or a wall keeps the Darcy-Ergun law
            const bool inertiaB = inertia && pB.fixesValue();
            if (!inertiaB)
            {
                lawPtr = dragModel_.get();
            }
            const vectorField gIB
            (
                inertiaB
              ? vectorField(gInertia_->boundaryField()[patchi])
              : vectorField(DS.size(), Zero)
            );
            const scalarField gn(pB.snGrad());
            // Closed outlets are walls for the gas too (WP F.4)
            const word& pname = mesh.boundary()[patchi].name();
            const bool shut =
                solid_.isOutlet(pname) && !solid_.outletOpen(pname);
            forAll(DS, i)
            {
                const label P = fc[i];
                if (shut || (fs && vd[P] > 0.5))
                {
                    // Gas flow in the void (WP F.6): rho_g 1/A
                    DS[i] =
                        (!shut && flow_)
                      ? max(p0_ + pB[i], pFloor)/RT
                       *rAUfPtr_().boundaryField()[patchi][i]
                      : 0;
                    DT[i] = DS[i];
                    continue;
                }
                face
                (
                    eM[P], k[P], eM[P], k[P], 1,
                    max(p0_ + pB[i], pFloor)/RT,
                    mag3(gradP[P] - gIB[i], n[i], gn[i] - (n[i] & gIB[i])),
                    DS[i], DT[i]
                );
            }
        }
    }

    if (const surfaceScalarField* faceOpen = solid_.faceOpenPtr())
    {
        DfS *= *faceOpen;
        DfT *= *faceOpen;
    }
}


void Foam::poreGas::solve(const scalar dt)
{
    const fvMesh& mesh = solid_.mesh();
    volScalarField& p = *pPtr_;
    const scalar RT = Rair*T_;

    const scalarField& phiG = phiG_.primitiveField();
    const tmp<scalarField> tOpen(solid_.cellOpen());
    const scalarField& open = tOpen();

    // Gas storage per cell volume and pressure [kg/(m^3 Pa)]. Cells without
    // gas (inside a body) get 1/RT so the matrix stays regular; they have
    // no mobility, so their p never changes and they hold no mass.
    scalarField S(open*phiG/RT);
    if (kappa_ > 0)
    {
        // WP F.8e: storage that gives the Darcy term the diffusivity kappa
        // in the filled part of each cell
        const scalarField& alpha = solid_.alpha().primitiveField();
        const scalarField& kc = k_.primitiveField();
        forAll(S, celli)
        {
            S[celli] = open[celli]*alpha[celli]*rhoG_*kc[celli]/(mu_*kappa_);
        }
    }
    forAll(S, celli)
    {
        if (S[celli] <= 0)
        {
            S[celli] = 1/RT;
        }
    }
    const scalarField pOld(p.primitiveField());

    const bool fs = freeSurface_;
    const scalarField& vd = void_.primitiveField();
    // Void cells held at the void pressure (WP F.4); with gas flow in the
    // void (WP F.6) they solve for their own
    const bool fixVoid = fs && !flow_;

    // Gas carried and squeezed by the grains (WP F.3, decision 0128): per
    // cell the gas mass X [kg] that the solid's step would add to the gas
    // stored at the old pressure,
    //   X = (phiG - phiG0) Vo rho_g,P + sum_f rho_g,donor G_f,
    //   G_f = dt phiL_f - dm_f/rhoP (gas volume out across the face),
    // with rho_g at the start of the step. The p0 part and the gauge part
    // are summed apart (round-off of p0). The balance below is then
    //   S Vo/dt (p - pOld) + X/dt + sum(Darcy flux out) = 0.
    scalarField X(p.size(), Zero);
    scalar carriedOut = 0;      // through non-coupled faces (WP F.4) [kg]
    if (kappa_ > 0 && coupled_ && solid_.mesh().time().timeIndex() > 0)
    {
        // WP F.8e: the material crossing a face (its bulk volume, mass over
        // the donor's bulk density) brings the donor's pressure into the
        // receiving cell; the cell it leaves keeps its own (upwind
        // advection, non-conservative as dp/dt + u.grad p)
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& dm = solid_.remapMass().primitiveField();
        const scalarField& rhoC = solid_.rho().primitiveField();
        const scalarField& alpha = solid_.alpha().primitiveField();
        auto bulk = [&](const label c)
        {
            return alpha[c] > alphaMin ? rhoC[c]/alpha[c] : rhoP_;
        };
        forAll(own, facei)
        {
            const scalar m = dm[facei];
            if (m == 0)
            {
                continue;
            }
            const label d = m > 0 ? own[facei] : nei[facei];
            const label r = m > 0 ? nei[facei] : own[facei];
            const scalar F = mag(m)/bulk(d);
            X[r] -= S[r]/max(alpha[r], alphaMin)*F*(pOld[d] - pOld[r]);
        }
        forAll(mesh.boundary(), patchi)
        {
            if (!mesh.boundary()[patchi].coupled())
            {
                continue;
            }
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const scalarField& pdm =
                solid_.remapMass().boundaryField()[patchi];
            const scalarField pN
            (
                p.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField rN
            (
                solid_.rho().boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField aN
            (
                solid_.alpha().boundaryField()[patchi].patchNeighbourField()
            );
            forAll(fc, i)
            {
                // Material coming in from the neighbour processor
                if (pdm[i] < 0)
                {
                    const label r = fc[i];
                    const scalar F =
                        -pdm[i]/(aN[i] > alphaMin ? rN[i]/aN[i] : rhoP_);
                    X[r] -= S[r]/max(alpha[r], alphaMin)*F*(pN[i] - pOld[r]);
                }
            }
        }
    }
    else if (coupled_ && solid_.mesh().time().timeIndex() > 0)
    {
        const scalarField Vo(open*mesh.V().field());
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& phiL = solid_.phiL().primitiveField();
        const scalarField& dm = solid_.remapMass().primitiveField();

        scalarField vol((phiG - phiG0_)*Vo);    // p0 part, volumes
        scalarField gauge(vol*pOld);            // gauge part, volume x Pa
        // A face between a gas cell and a void cell (WP F.4) is open to
        // the void: the surface moving does not squeeze the gas, only the
        // grains crossing it displace their volume of gas (out of the void
        // when they leave the cell, into it when they come in)
        // With gas flow in the void (WP F.6) the gas there moves by its
        // own velocity: on faces of void cells only the grains crossing
        // displace gas
        // With surfaceCarry (WP F.8, decision 0133) the surface's own
        // motion carries gas too: the part of a cell that the material
        // leaves fills with gas from the void across the surface face, and
        // material rising into a cell pushes the gas there into the void,
        // without a pressure drop through the grains. Only faces between
        // two void cells (gas flow there moves by its own velocity) keep
        // the grains-only rule.
        auto surfaceFace = [&](const bool vA, const bool vB)
        {
            if (surfaceCarry_)
            {
                return fs && flow_ && vA && vB;
            }
            return fs && (flow_ ? (vA || vB) : vA != vB);
        };
        forAll(own, facei)
        {
            const scalar G =
                surfaceFace(vd[own[facei]] > 0.5, vd[nei[facei]] > 0.5)
              ? -dm[facei]/rhoP_
              : dt*phiL[facei] - dm[facei]/rhoP_;
            const scalar pd = G > 0 ? pOld[own[facei]] : pOld[nei[facei]];
            vol[own[facei]] += G;
            vol[nei[facei]] -= G;
            gauge[own[facei]] += G*pd;
            gauge[nei[facei]] -= G*pd;
        }
        forAll(mesh.boundary(), patchi)
        {
            if (!mesh.boundary()[patchi].coupled())
            {
                // Non-coupled boundaries: no mixture volume crosses them in
                // the remap of full cells. With a free surface (WP F.4) the
                // Lagrangian volume counts every face: gas leaves with the
                // material through outlets at the cell's pressure and comes
                // in through inlets (or open faces) at the boundary's
                if (!fs)
                {
                    continue;
                }
                // Only faces open to the gas outside carry it: open outlets
                // and inlets, and fixedValue patches of p. On the others
                // (walls) the material that pulls away opens void in the
                // cell, which the cell's gas fills: no gas crosses them.
                const word& pname = mesh.boundary()[patchi].name();
                const bool open =
                    solid_.isOutlet(pname)
                  ? solid_.outletOpen(pname)
                  : (
                        solid_.inletOpen(pname)
                     || p.boundaryField()[patchi].fixesValue()
                    );
                if (!open)
                {
                    continue;
                }
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                const scalarField& pPhi = solid_.phiL().boundaryField()[patchi];
                const scalarField& pdm =
                    solid_.remapMass().boundaryField()[patchi];
                const scalarField& pB = p.boundaryField()[patchi];
                forAll(fc, i)
                {
                    const scalar G =
                        (flow_ && vd[fc[i]] > 0.5)
                      ? -pdm[i]/rhoP_
                      : dt*pPhi[i] - pdm[i]/rhoP_;
                    if (G == 0)
                    {
                        continue;
                    }
                    const scalar pd = G > 0 ? pOld[fc[i]] : pB[i];
                    vol[fc[i]] += G;
                    gauge[fc[i]] += G*pd;
                    carriedOut += (p0_*G + G*pd)/RT;
                }
                continue;
            }
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const scalarField& pPhi = solid_.phiL().boundaryField()[patchi];
            const scalarField& pdm =
                solid_.remapMass().boundaryField()[patchi];
            const scalarField pN
            (
                p.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField vN
            (
                fs ? void_.boundaryField()[patchi].patchNeighbourField()
                   : tmp<scalarField>::New(fc.size(), Zero)
            );
            forAll(fc, i)
            {
                const scalar G =
                    surfaceFace(vd[fc[i]] > 0.5, vN[i] > 0.5)
                  ? -pdm[i]/rhoP_
                  : dt*pPhi[i] - pdm[i]/rhoP_;
                const scalar pd = G > 0 ? pOld[fc[i]] : pN[i];
                vol[fc[i]] += G;
                gauge[fc[i]] += G*pd;
            }
        }
        forAll(X, celli)
        {
            X[celli] = (p0_*vol[celli] + gauge[celli])/RT;
        }
    }

    // Void cells hold the void pressure (WP F.4); a cell that has just
    // become void gives its gas to the atmosphere (counted below). After
    // the carried gas, which takes the old pressure on both sides of every
    // face (processor faces too)
    if (fixVoid)
    {
        scalarField& pc = p.primitiveFieldRef();
        forAll(pc, celli)
        {
            if (vd[celli] > 0.5)
            {
                pc[celli] = 0;
            }
        }
        p.correctBoundaryConditions();
    }

    // First guess (WP F.3): the pressure extrapolated linearly from the
    // last two steps, so a smooth evolution needs fewer iterations. The
    // converged result does not depend on it beyond pTolerance.
    if (dtOld_ > 0 && pOldOld_.size() == pOld.size())
    {
        p.primitiveFieldRef() += (pOld - pOldOld_)*(dt/dtOld_);
        if (fixVoid)
        {
            scalarField& pc = p.primitiveFieldRef();
            forAll(pc, celli)
            {
                if (vd[celli] > 0.5)
                {
                    pc[celli] = 0;
                }
            }
        }
        p.correctBoundaryConditions();
    }
    pOldOld_ = pOld;
    dtOld_ = dt;

    volScalarField::Internal SdtV
    (
        IOobject("gasStorage", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimDensity/dimPressure/dimTime, Zero)
    );
    SdtV.field() = S/dt;

    const dimensionSet dimD(dimDensity*dimArea/dimPressure/dimTime);
    surfaceScalarField DfS
    (
        IOobject("gasDiffusivity", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimD, Zero)
    );
    surfaceScalarField DfT
    (
        IOobject("gasDiffusivityTangent", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimD, Zero)
    );

    tmp<fv::laplacianScheme<scalar, scalar>> tLap;
    {
        IStringStream is(laplacian_);
        tLap = fv::laplacianScheme<scalar, scalar>::New(mesh, is);
    }

    const scalar rhoG0 = rhoG_;
    label outer = 0;
    scalar change = GREAT, prevChange = GREAT, errorEst = GREAT;
    label nNewton = 0;
    label nInner = 0;
    label nOuterAll = 0;

    // Gas flow in the void (WP F.6): momentum predictor, then nCorr
    // pressure solves (PISO), each with H/A of the latest velocity
    if (flow_)
    {
        predictMomentum(dt, pOld);
    }
    const label nCorr = flow_ ? max(nCorr_, label(1)) : 1;

    // Pore gas momentum (WP F.7): the old relative fluxes of this step
    startInertia(dt);

    for (label corr = 0; corr < nCorr; ++corr)
    {
    if (flow_)
    {
        flowFluxes();
    }
    outer = 0;
    change = GREAT;
    prevChange = GREAT;
    errorEst = GREAT;

    // Anderson history
    DynamicList<scalarField> dF, dG;
    scalarField fPrev, gPrev;
    bool lastNewton = false;
    scalarField delta(p.size(), Zero);

    for (outer = 0; outer < nOuter_; ++outer)
    {
        // Face coefficients at the current pressure gradient (Ergun's
        // inertial term through the slip) and gas density
        faceCoefficients
        (
            fv::gaussGrad<scalar>(mesh).calcGrad(p, "grad(p)"),
            DfS,
            DfT
        );

        // Inflow: the gradient that carries rho_g0 U in through the face
        // (a fixed mass flux, so it adds nothing to the Newton matrix)
        forAllConstIters(inflow_, iter)
        {
            const label patchi = mesh.boundaryMesh().findPatchID(iter.key());
            auto& pB = refCast<fixedGradientFvPatchScalarField>
            (
                p.boundaryFieldRef()[patchi]
            );
            const scalarField& D = DfS.boundaryField()[patchi];
            scalarField& g = pB.gradient();
            const scalar Uin = iter.val()->value(mesh.time().value());
            forAll(g, i)
            {
                g[i] = D[i] > 0 ? rhoG0*Uin/D[i] : 0;
            }
            pB.evaluate();
        }

        // Mass balance with the flux law at the current p: its residual
        fvScalarMatrix eqnS
        (
            fvm::Sp(SdtV, p)
          - tLap.ref().fvmLaplacian(DfS, p)
        );
        eqnS.source() += SdtV.field()*mesh.V().field()*pOld - X/dt;
        if (flow_)
        {
            // Gas carried by H/A through faces of void cells (WP F.6)
            eqnS.source() -=
                fvc::surfaceIntegrate(mHbyAPtr_())().primitiveField()
               *mesh.V().field();
        }
        if (poreInertia_)
        {
            // Gas carried by the old pore-gas momentum (WP F.7)
            inertiaFlux(DfS);
            eqnS.source() -=
                fvc::surfaceIntegrate(mInertia_())().primitiveField()
               *mesh.V().field();
        }

        if (outer > 0 && errorEst <= pTol_)
        {
            massFlux_ = eqnS.flux();
            break;
        }
        // Residual b - A p. (fvMatrix::residual() adds the processor-face
        // terms a second time on top of the interface update, so it is
        // assembled here as fvMatrix::solve does: boundary diagonal and
        // the sources of the non-coupled patches only.)
        scalarField r;
        {
            scalarField src(eqnS.source());
            scalarField bDiag(p.size(), Zero);
            forAll(eqnS.internalCoeffs(), patchi)
            {
                const labelUList& fc = eqnS.lduAddr().patchAddr(patchi);
                const scalarField& ic = eqnS.internalCoeffs()[patchi];
                forAll(fc, i)
                {
                    bDiag[fc[i]] += ic[i];
                }
                if (!p.boundaryField()[patchi].coupled())
                {
                    const scalarField& bc = eqnS.boundaryCoeffs()[patchi];
                    forAll(fc, i)
                    {
                        src[fc[i]] += bc[i];
                    }
                }
            }
            src -= bDiag*p.primitiveField();
            r = eqnS.lduMatrix::residual
            (
                p.primitiveField(),
                src,
                eqnS.boundaryCoeffs(),
                p.boundaryField().scalarInterfaces(),
                0
            );
        }
        if (fixVoid)
        {
            // Void cells are fixed (WP F.4)
            forAll(r, celli)
            {
                if (vd[celli] > 0.5)
                {
                    r[celli] = 0;
                }
            }
        }

        // Matrix for the correction: Newton (tangent mobility) once the
        // corrections shrink; Picard (secant mobility, the flux law
        // itself) on the first iteration and whenever the last
        // correction grew. Newton alone can overshoot from a far start
        // (Ergun's q ~ sqrt(|grad p|)), Picard alone converges slowly.
        const bool newton = outer > 0 && change < prevChange;
        if (newton)
        {
            ++nNewton;
        }
        fvScalarMatrix eqnT
        (
            fvm::Sp(SdtV, p)
          - tLap.ref().fvmLaplacian(newton ? DfT : DfS, p)
        );
        if (fixVoid)
        {
            // Void cells are fixed (WP F.4): no coupling to them, so their
            // correction is zero and a gas cell beside one sees it as a
            // fixed pressure (its face coefficient stays on its diagonal);
            // the matrix stays symmetric, across processors too
            const labelUList& own = mesh.owner();
            const labelUList& nei = mesh.neighbour();
            scalarField& up = eqnT.upper();
            forAll(up, facei)
            {
                if (vd[own[facei]] > 0.5 || vd[nei[facei]] > 0.5)
                {
                    up[facei] = 0;
                }
            }
            forAll(mesh.boundary(), patchi)
            {
                if (!p.boundaryField()[patchi].coupled())
                {
                    continue;
                }
                const scalarField vN
                (
                    void_.boundaryField()[patchi].patchNeighbourField()
                );
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                scalarField& bc = eqnT.boundaryCoeffs()[patchi];
                forAll(fc, i)
                {
                    if (vd[fc[i]] > 0.5 || vN[i] > 0.5)
                    {
                        bc[i] = 0;
                    }
                }
            }
        }
        {
            scalarField saveDiag(eqnT.diag());
            // Boundary contributions to the diagonal, as fvMatrix::solve
            scalarField& diag = eqnT.diag();
            forAll(eqnT.internalCoeffs(), patchi)
            {
                const labelUList& fc = eqnT.lduAddr().patchAddr(patchi);
                const scalarField& ic = eqnT.internalCoeffs()[patchi];
                forAll(fc, i)
                {
                    diag[fc[i]] += ic[i];
                }
            }

            delta = Zero;
            const SolverPerformance<scalar> perf =
                lduMatrix::solver::New
                (
                    p.name() + "Correction",
                    eqnT,
                    eqnT.boundaryCoeffs(),
                    eqnT.internalCoeffs(),
                    p.boundaryField().scalarInterfaces(),
                    solverDict_
                )->solve(delta, r);
            nInner += perf.nIterations();

            eqnT.diag() = saveDiag;
        }

        prevChange = change;
        change = gMax(mag(delta)());

        // Error left after this correction (WP F.3): the next correction,
        // estimated from the contraction of the last two when they shrink
        // fast (linear convergence of the lagged gas density); without
        // that estimate, the correction itself
        {
            const scalar rate = outer > 0 ? change/prevChange : 1;
            errorEst = rate < 0.1 ? rate*change : change;
        }

        // Anderson acceleration (depth nAnderson) of the fixed-point map
        // p -> p + delta: p_new = g_k - dG gamma, gamma minimising
        // |f_k - dF gamma| with f = delta, g = p + delta. A combination
        // with weights summing to one keeps the gas mass of the iterates.
        // History restarts when the matrix changes kind (Picard/Newton).
        scalarField g(p.primitiveField() + delta);
        if (nAnderson_ > 0)
        {
            if (newton != lastNewton || outer == 0)
            {
                dF.clear();
                dG.clear();
            }
            else
            {
                dF.append(scalarField(delta - fPrev));
                dG.append(scalarField(g - gPrev));
                if (dF.size() > nAnderson_)
                {
                    dF.remove(0);
                    dG.remove(0);
                }
            }
            fPrev = delta;
            gPrev = g;
            lastNewton = newton;

            const label m = dF.size();
            if (m > 0)
            {
                // Normal equations H gamma = b, all processors
                scalarSquareMatrix H(m, Zero);
                scalarField b(m, Zero);
                for (label i = 0; i < m; ++i)
                {
                    for (label j = 0; j <= i; ++j)
                    {
                        H(i, j) = H(j, i) = sum(dF[i]*dF[j]);
                    }
                    b[i] = sum(dF[i]*delta);
                }
                {
                    scalarList flat(m*m);
                    for (label i = 0; i < m; ++i)
                    {
                        for (label j = 0; j < m; ++j) flat[i*m + j] = H(i, j);
                    }
                    Pstream::listCombineReduce(flat, plusEqOp<scalar>());
                    for (label i = 0; i < m; ++i)
                    {
                        for (label j = 0; j < m; ++j) H(i, j) = flat[i*m + j];
                    }
                }
                Pstream::listCombineReduce(b, plusEqOp<scalar>());
                for (label i = 0; i < m; ++i)
                {
                    H(i, i) *= 1 + 1e-10;
                }
                // Gaussian elimination with partial pivoting (m <= 10)
                scalarField gamma(b);
                bool ok = true;
                for (label c = 0; c < m && ok; ++c)
                {
                    label piv = c;
                    for (label r = c + 1; r < m; ++r)
                    {
                        if (mag(H(r, c)) > mag(H(piv, c))) piv = r;
                    }
                    if (mag(H(piv, c)) < VSMALL) { ok = false; break; }
                    if (piv != c)
                    {
                        for (label j = 0; j < m; ++j) std::swap(H(c, j), H(piv, j));
                        std::swap(gamma[c], gamma[piv]);
                    }
                    for (label r = c + 1; r < m; ++r)
                    {
                        const scalar f = H(r, c)/H(c, c);
                        for (label j = c; j < m; ++j) H(r, j) -= f*H(c, j);
                        gamma[r] -= f*gamma[c];
                    }
                }
                if (ok)
                {
                    for (label c = m - 1; c >= 0; --c)
                    {
                        for (label j = c + 1; j < m; ++j) gamma[c] -= H(c, j)*gamma[j];
                        gamma[c] /= H(c, c);
                    }
                    for (label i = 0; i < m; ++i)
                    {
                        g -= gamma[i]*dG[i];
                    }
                }
            }
        }

        p.primitiveFieldRef() = g;
        p.correctBoundaryConditions();
    }

    if (outer == nOuter_)
    {
        // Not converged: flux of the last state for the bookkeeping
        fvScalarMatrix eqnS
        (
            fvm::Sp(SdtV, p)
          - tLap.ref().fvmLaplacian(DfS, p)
        );
        massFlux_ = eqnS.flux();
    }
    nOuterAll += outer;

    if (poreInertia_)
    {
        // The old pore-gas momentum's part (WP F.7)
        massFlux_ += mInertia_();
    }
    if (flow_)
    {
        // Total gas mass flux: the pressure part and H/A (WP F.6)
        massFlux_ += mHbyAPtr_();
        if (corr < nCorr - 1)
        {
            correctVelocity();
        }
    }
    }   // pressure corrections (WP F.6)

    // Conservative finish: the new pressure from the old one and the face
    // fluxes of the converged state, so the gas mass changes by exactly
    // the boundary fluxes (to round-off) whatever is left of the
    // iteration (at most the residual, a fraction of pTolerance)
    {
        const scalarField netOut
        (
            fvc::surfaceIntegrate(massFlux_)().primitiveField()
        );
        scalarField& pc = p.primitiveFieldRef();
        scalar fromVoid = 0;
        forAll(pc, celli)
        {
            if (fixVoid && vd[celli] > 0.5)
            {
                // Void at the void pressure (WP F.4): the atmosphere gives
                // what its balance lacks
                pc[celli] = 0;
                fromVoid +=
                    dt*mesh.V()[celli]
                   *(SdtV[celli]*(0 - pOld[celli]) + netOut[celli])
                  + X[celli];
                continue;
            }
            pc[celli] =
                pOld[celli]
              - (netOut[celli] + X[celli]/(dt*mesh.V()[celli]))/SdtV[celli];
        }
        p.correctBoundaryConditions();
        if (fixVoid)
        {
            fromVoid = returnReduce(fromVoid, sumOp<scalar>());
            voidIn_ += fromVoid;
            netIn_ += fromVoid;
        }
    }

    // Gas that came in through the (non-coupled) boundaries
    scalar in = 0;
    forAll(massFlux_.boundaryField(), patchi)
    {
        if (!massFlux_.boundaryField()[patchi].coupled())
        {
            in -= sum(massFlux_.boundaryField()[patchi]);
        }
    }
    netIn_ += returnReduce(in, sumOp<scalar>())*dt;
    if (fs)
    {
        netIn_ -= returnReduce(carriedOut, sumOp<scalar>());
    }
    correctFacePressure();

    // Superficial flux q = massFlux / (rho_g |Sf|), reconstructed to cells
    {
        const dimensionedScalar dRT(dimPressure/dimDensity, RT);
        const dimensionedScalar dp0(dimPressure, p0_);
        surfaceScalarField volFlux
        (
            massFlux_/((fvc::interpolate(p) + dp0)/dRT)
        );
        qGas_ = fvc::reconstruct(volFlux);
        qGas_.correctBoundaryConditions();
    }
    endInertia(dt);     // WP F.7
    if (flow_)
    {
        correctVelocity();
    }

    const label oldPrecision = Info().precision(16);
    if (fixVoid)
    {
        Info<< "Pore gas void: " << gSum(vd) << " cells, from the void "
            << voidIn_ << " kg" << nl;
    }
    if (flow_)
    {
        const volVectorField& Ug = *UgPtr_;
        scalar uMax = 0;
        forAll(vd, celli)
        {
            if (vd[celli] > 0.5)
            {
                uMax = max(uMax, mag(Ug[celli]));
            }
        }
        Info<< "Gas flow: " << nCorr << " corrections, largest |Ug| in the "
            << "void " << returnReduce(uMax, maxOp<scalar>()) << " m/s" << nl;
    }
    Info<< "Pore gas: " << (flow_ ? nOuterAll : outer) << " solves (" << nNewton << " Newton, "
        << nInner << " iterations), change " << change
        << " Pa, mass " << gasMass() << " kg, net in " << netIn_ << " kg"
        << nl;
    Info().precision(oldPrecision);

    if (outer == nOuter_)
    {
        WarningInFunction
            << "pore pressure not converged in " << nOuter_
            << " iterations (last change " << change << " Pa)" << endl;
    }
}


void Foam::poreGas::correctFacePressure()
{
    // Linear between cells (F.3); the void pressure, 0 gauge, on faces
    // between material and void (WP F.4)
    pf_ = fvc::interpolate(*pPtr_);
    if (!freeSurface_)
    {
        return;
    }
    const fvMesh& mesh = solid_.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& vd = void_.primitiveField();
    scalarField& pfi = pf_.primitiveFieldRef();
    const scalarField& pc = pPtr_->primitiveField();
    forAll(pfi, facei)
    {
        const bool vO = vd[own[facei]] > 0.5, vN = vd[nei[facei]] > 0.5;
        if (vO || vN)
        {
            // With gas flow (WP F.6) the void cell's pressure
            pfi[facei] =
                !flow_ ? 0
              : (vO && vN) ? pfi[facei]
              : vO ? pc[own[facei]] : pc[nei[facei]];
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        scalarField& pp = pf_.boundaryFieldRef()[patchi];
        if (mesh.boundary()[patchi].coupled())
        {
            const scalarField vN
            (
                void_.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField pN
            (
                pPtr_->boundaryField()[patchi].patchNeighbourField()
            );
            forAll(fc, i)
            {
                const bool vO = vd[fc[i]] > 0.5, vB = vN[i] > 0.5;
                if (vO || vB)
                {
                    pp[i] =
                        !flow_ ? 0
                      : (vO && vB) ? pp[i]
                      : vO ? pc[fc[i]] : pN[i];
                }
            }
        }
    }
}


// ************************************************************************* //
