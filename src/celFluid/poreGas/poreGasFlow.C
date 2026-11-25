/*---------------------------------------------------------------------------*\
    poreGas: gas momentum in the void (WP F.6, decisions 0125, 0131)

    With `model gasFlow;` in constant/gasProperties the void cells (alpha
    below alphaVoid of a freeSurface solid) carry a gas velocity Ug and
    their own pressure. Isothermal ideal gas, low Mach: per step

        rho_o (Ug - Ug_o)/dt + div(m Ug) - Ug div(m) - div(mu grad Ug)
            = -grad p

    (m the gas mass flux of the last step), solved once as a predictor;
    then nCorr pressure solves (PISO): the gas mass balance of poreGas with
    the face flux of the void faces m_f = rho_f (HbyA_f . S_f - rAU_f
    |S_f| snGrad p), so the pressure solve is the same Newton solve as the
    pore gas, with rho_f rAU_f as the face coefficient of the void faces and
    the H/A flux as a source. The time step stays the solid's (the
    pressure is implicit, so the gas sound speed sets no limit). Faces
    between the bed and the void keep the F.4 rule (the pore side's
    mobility over its half distance), now against the void cell's own
    pressure. Cells with grains are not solved for Ug: they hold the grain
    velocity plus the pore flux over phi_g (the void sees them as moving
    walls). The weight of the gas is left out, as in the pore gas (12 Pa/m
    against about 16 kPa/m for the grains). No drag on grains (F.7).
\*---------------------------------------------------------------------------*/

#include "poreGas.H"
#include "Time.H"
#include "fvm.H"
#include "fvc.H"
#include "gaussGrad.H"
#include "laplacianScheme.H"
#include "convectionScheme.H"
#include "constrainHbyA.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "IStringStream.H"
#include "linear.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::poreGas::initFlow()
{
    const fvMesh& mesh = solid_.mesh();

    if (!freeSurface_)
    {
        FatalIOErrorInFunction(dict_)
            << "model gasFlow needs a freeSurface solid: the gas moves in "
            << "the void cells (alpha below alphaVoid)" << exit(FatalIOError);
    }

    IOobject io
    (
        "Ug",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    if (io.typeHeaderOk<volVectorField>(true))
    {
        UgPtr_.reset(new volVectorField(io, mesh));
    }
    else
    {
        // Default: no slip on walls, zeroGradient where p is fixed (open
        // to the outside), constraint patches as they are
        const volScalarField& p = *pPtr_;
        wordList types
        (
            mesh.boundary().size(),
            fixedValueFvPatchVectorField::typeName
        );
        forAll(types, patchi)
        {
            const word& pt = mesh.boundary()[patchi].type();
            if (polyPatch::constraintType(pt))
            {
                types[patchi] = pt;
            }
            else if (p.boundaryField()[patchi].fixesValue())
            {
                types[patchi] = zeroGradientFvPatchVectorField::typeName;
            }
        }
        io.readOpt(IOobject::NO_READ);
        UgPtr_.reset
        (
            new volVectorField
            (
                io,
                mesh,
                dimensionedVector(dimVelocity, Zero),
                types
            )
        );
    }

    USolverDict_.add("solver", "PBiCGStab");
    USolverDict_.add("preconditioner", "DILU");
    USolverDict_.add("tolerance", 1e-30);
    USolverDict_.add("relTol", 1e-10);
    USolverDict_.add("maxIter", 1000);
    USolverDict_.add("log", 0);
    if (const dictionary* sd = dict_.findDict("USolver"))
    {
        USolverDict_.merge(*sd);
    }

    // Mass flux of the starting velocity through the void faces, for the
    // convection of the first step (a restart starts from Ug)
    {
        const scalar RT = Rair*T_;
        const volScalarField& p = *pPtr_;
        const surfaceScalarField rhof
        (
            (linearInterpolate(p) + dimensionedScalar(dimPressure, p0_))
           /dimensionedScalar(dimPressure/dimDensity, RT)
        );
        massFlux_ = voidFaces()*rhof*(linearInterpolate(*UgPtr_) & mesh.Sf());
    }

    Info<< "Gas flow in the void (WP F.6): " << nCorr_
        << " pressure corrections per step, convection " << divUg_
        << ", Ug " << (io.readOpt() == IOobject::MUST_READ ? "read" : "0")
        << nl << endl;
}


Foam::tmp<Foam::surfaceScalarField> Foam::poreGas::voidFaces() const
{
    // 1 on faces of the void flow, else 0. Internal face: both cells void.
    // Coupled patch face: both sides void. Other patch face: the cell void
    // and, for an outlet, open
    const fvMesh& mesh = solid_.mesh();
    const scalarField& vd = void_.primitiveField();
    tmp<surfaceScalarField> tm
    (
        new surfaceScalarField
        (
            IOobject("gasVoidFaces", mesh.time().timeName(), mesh),
            mesh,
            dimensionedScalar(dimless, Zero)
        )
    );
    surfaceScalarField& m = tm.ref();
    scalarField& mi = m.primitiveFieldRef();
    forAll(mi, facei)
    {
        mi[facei] =
            (vd[mesh.owner()[facei]] > 0.5 && vd[mesh.neighbour()[facei]] > 0.5)
          ? 1 : 0;
    }
    forAll(m.boundaryField(), patchi)
    {
        const fvPatch& pp = mesh.boundary()[patchi];
        const labelUList& fc = pp.faceCells();
        scalarField& mb = m.boundaryFieldRef()[patchi];
        if (pp.coupled())
        {
            const scalarField vN
            (
                void_.boundaryField()[patchi].patchNeighbourField()
            );
            forAll(mb, i)
            {
                mb[i] = (vd[fc[i]] > 0.5 && vN[i] > 0.5) ? 1 : 0;
            }
        }
        else
        {
            const word& name = pp.name();
            const bool shut =
                solid_.isOutlet(name) && !solid_.outletOpen(name);
            forAll(mb, i)
            {
                mb[i] = (!shut && vd[fc[i]] > 0.5) ? 1 : 0;
            }
        }
    }
    return tm;
}


void Foam::poreGas::predictMomentum(const scalar dt, const scalarField& pOld)
{
    const fvMesh& mesh = solid_.mesh();
    volVectorField& U = *UgPtr_;
    const volScalarField& p = *pPtr_;
    const scalar RT = Rair*T_;
    const scalarField& vd = void_.primitiveField();

    // Velocities of the cells that are not solved: grain velocity plus
    // the pore flux over phi_g; none inside bodies
    {
        const tmp<scalarField> tOpen(solid_.cellOpen());
        const scalarField& open = tOpen();
        const vectorField& Us = solid_.U().primitiveField();
        const scalarField& phiG = phiG_.primitiveField();
        vectorField& Uc = U.primitiveFieldRef();
        forAll(Uc, celli)
        {
            if (vd[celli] < 0.5)
            {
                Uc[celli] =
                    open[celli] < eulerianSolid::coveredOpen
                  ? vector::zero
                  : Us[celli] + qGas_[celli]/max(phiG[celli], SMALL);
            }
        }
        U.correctBoundaryConditions();
    }
    const vectorField U0(U.primitiveField());

    // Old velocity and volume flux of the void faces, for the time
    // derivative's part of the face flux (ddtCorr, Rhie-Chow): without it
    // a pressure checkerboard in the void survives (seen in the venting
    // column, 3 Pa odd-even)
    Uold_.reset(new volVectorField("UgOld", U));
    {
        const surfaceScalarField rhof
        (
            (linearInterpolate(p) + dimensionedScalar(dimPressure, p0_))
           /dimensionedScalar(dimPressure/dimDensity, RT)
        );
        phiOld_.reset(new surfaceScalarField("phiUgOld", massFlux_/rhof));
    }

    volScalarField::Internal rhoDt
    (
        IOobject("gasRhoByDt", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
    rhoDt.field() = (p0_ + pOld)/(RT*dt);

    volScalarField::Internal divM
    (
        IOobject("gasDivMassFlux", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimDensity/dimTime, Zero)
    );
    divM.field() = fvc::surfaceIntegrate(massFlux_)().primitiveField();

    tmp<fv::convectionScheme<vector>> tConv;
    {
        IStringStream is(divUg_);
        tConv = fv::convectionScheme<vector>::New(mesh, massFlux_, is);
    }
    tmp<fv::laplacianScheme<vector, scalar>> tLap;
    {
        IStringStream is(laplacian_);
        tLap = fv::laplacianScheme<vector, scalar>::New(mesh, is);
    }
    const surfaceScalarField muf
    (
        IOobject("gasMu", mesh.time().timeName(), mesh),
        mesh,
        dimensionedScalar(dimDynamicViscosity, mu_)
    );

    UEqnPtr_.reset
    (
        new fvVectorMatrix
        (
            fvm::Sp(rhoDt, U)
          + tConv.ref().fvmDiv(massFlux_, U)
          - fvm::Sp(divM, U)
          - tLap.ref().fvmLaplacian(muf, U)
        )
    );
    fvVectorMatrix& UEqn = UEqnPtr_();
    UEqn.source() += rhoDt.field()*mesh.V().field()*U0;

    // Cells with grains (or inside bodies) hold their velocity
    {
        DynamicList<label> cells(vd.size());
        forAll(vd, celli)
        {
            if (vd[celli] < 0.5)
            {
                cells.append(celli);
            }
        }
        UEqn.setValues(cells, vectorField(U0, cells));
    }

    // Predictor with the pressure gradient of the first guess
    fvVectorMatrix UPred(UEqn);
    {
        const tmp<vectorField> tgradP(flowGradP());
        const vectorField& gradP = tgradP();
        vectorField& src = UPred.source();
        forAll(vd, celli)
        {
            if (vd[celli] > 0.5)
            {
                src[celli] -= mesh.V()[celli]*gradP[celli];
            }
        }
    }
    UPred.solve(USolverDict_);
}


void Foam::poreGas::flowFluxes()
{
    const fvMesh& mesh = solid_.mesh();
    const volVectorField& U = *UgPtr_;
    const volScalarField& p = *pPtr_;
    const fvVectorMatrix& UEqn = UEqnPtr_();
    const scalar RT = Rair*T_;

    rAUPtr_.reset(new volScalarField("rAUg", 1/UEqn.A()));
    rAUPtr_->correctBoundaryConditions();
    HbyAPtr_.reset
    (
        new volVectorField(constrainHbyA(rAUPtr_()*UEqn.H(), U, p))
    );
    rAUfPtr_.reset(new surfaceScalarField("rAUgf", linearInterpolate(rAUPtr_())));

    const surfaceScalarField rhof
    (
        (linearInterpolate(p) + dimensionedScalar(dimPressure, p0_))
       /dimensionedScalar(dimPressure/dimDensity, RT)
    );
    // Face flux of H/A with the time derivative's part from the old face
    // flux instead of the interpolated old velocity (ddtCorr): rho_f/dt
    // 1/A_f (phi_o - U_o,f . S_f)
    const scalar dt = mesh.time().deltaTValue();
    mHbyAPtr_.reset
    (
        new surfaceScalarField
        (
            "gasMassFluxHbyA",
            rhof
           *(
                (linearInterpolate(HbyAPtr_()) & mesh.Sf())
              + rhof*rAUfPtr_()/dimensionedScalar(dimTime, dt)
               *(phiOld_() - (linearInterpolate(Uold_()) & mesh.Sf()))
            )
        )
    );
    surfaceScalarField& m = mHbyAPtr_();
    m *= voidFaces();
    if (const surfaceScalarField* faceOpen = solid_.faceOpenPtr())
    {
        m *= *faceOpen;
    }
}


Foam::tmp<Foam::vectorField> Foam::poreGas::flowGradP()
{
    // Gauss gradient of the face pressures as the void sees them: linear
    // between cells, the void cell's own pressure on faces between the bed
    // and the void (pFace). With the linear face value there, the void
    // cell above a bed felt half the jump at the bed's surface as a
    // gradient no face flux balances, and its velocity grew without bound
    correctFacePressure();
    const fvMesh& mesh = solid_.mesh();
    return tmp<vectorField>::New
    (
        fvc::surfaceIntegrate(mesh.Sf()*pf_)().primitiveField()
    );
}


void Foam::poreGas::correctVelocity()
{
    volVectorField& U = *UgPtr_;
    const scalarField& vd = void_.primitiveField();
    const tmp<vectorField> tgradP(flowGradP());
    const vectorField& gradP = tgradP();
    const scalarField& rAU = rAUPtr_();
    const vectorField& HbyA = HbyAPtr_();
    const tmp<scalarField> tOpen(solid_.cellOpen());
    const scalarField& open = tOpen();
    const vectorField& Us = solid_.U().primitiveField();
    const scalarField& phiG = phiG_.primitiveField();

    vectorField& Uc = U.primitiveFieldRef();
    forAll(Uc, celli)
    {
        if (vd[celli] > 0.5)
        {
            Uc[celli] = HbyA[celli] - rAU[celli]*gradP[celli];
        }
        else
        {
            Uc[celli] =
                open[celli] < eulerianSolid::coveredOpen
              ? vector::zero
              : Us[celli] + qGas_[celli]/max(phiG[celli], SMALL);
        }
    }
    U.correctBoundaryConditions();
}


// ************************************************************************* //
