/*---------------------------------------------------------------------------*\
    eulerianSolid: see eulerianSolid.H
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

//- IOobject for a field read if present and written at output times
static IOobject fieldIO(const word& name, const fvMesh& mesh)
{
    return IOobject
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );
}

//- Gravity from constant/g if present, else zero
static vector readGravity(const fvMesh& mesh)
{
    IOobject io
    (
        "g",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );
    if (!io.typeHeaderOk<uniformDimensionedVectorField>(true))
    {
        return Zero;
    }
    return uniformDimensionedVectorField(io).value();
}


//- Wall patches and tan(friction) from materialProperties walls
static HashTable<scalar> readWalls(const fvMesh& mesh, const dictionary& dict)
{
    HashTable<scalar> walls;
    const dictionary* wd = dict.findDict("walls");
    if (!wd)
    {
        return walls;
    }
    for (const entry& e : *wd)
    {
        if (mesh.boundaryMesh().findPatchID(e.keyword()) < 0)
        {
            FatalIOErrorInFunction(*wd)
                << "walls: no patch " << e.keyword() << exit(FatalIOError);
        }
        const scalar deg = e.dict().get<scalar>("friction");
        if (deg < 0 || deg >= 90)
        {
            FatalIOErrorInFunction(*wd)
                << "walls: friction of " << e.keyword()
                << " must be in [0, 90) degrees, got " << deg
                << exit(FatalIOError);
        }
        walls.set(e.keyword(), Foam::tan(deg*constant::mathematical::pi/180));
    }
    return walls;
}

} // End namespace Foam


Foam::eulerianSolid::eulerianSolid(const fvMesh& mesh)
:
    mesh_(mesh),
    materialProperties_
    (
        IOobject
        (
            "materialProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    material_(constitutiveModel::New(materialProperties_)),
    maxCo_(mesh.time().controlDict().getOrDefault<scalar>("maxCo", 0.4)),
    solveMomentum_
    (
        mesh.time().controlDict().getOrDefault<Switch>("solveMomentum", true)
    ),
    secondOrder_
    (
        mesh.time().controlDict().getOrDefault<Switch>("secondOrder", true)
    ),
    rho_
    (
        fieldIO("rho", mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, materialProperties_),
        fvPatchFieldBase::zeroGradientType()
    ),
    U_
    (
        fieldIO("U", mesh),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    sigma_
    (
        fieldIO("sigma", mesh),
        mesh,
        dimensionedSymmTensor(dimPressure, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    epsP_
    (
        fieldIO("epsP", mesh),
        mesh,
        dimensionedSymmTensor(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    epsPEq_
    (
        fieldIO("epsPEq", mesh),
        mesh,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    phiL_
    (
        IOobject
        (
            "phiL",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVolume/dimTime, Zero)
    ),
    remapLimiter_(readRemapLimiters(materialProperties_)),
    g_(readGravity(mesh)),
    gravityRamp_
    (
        mesh.time().controlDict().getOrDefault<scalar>("gravityRamp", 0)
    ),
    damping_(mesh.time().controlDict().getOrDefault<scalar>("damping", 0)),
    wallFriction_(readWalls(mesh, materialProperties_)),
    vStar_
    (
        IOobject
        (
            "vStar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    ),
    Fstar_
    (
        IOobject
        (
            "Fstar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimForce, Zero)
    ),
    freeSurface_(materialProperties_.getOrDefault<Switch>("freeSurface", false)),
    alpha_
    (
        fieldIO("alpha", mesh),
        mesh,
        dimensionedScalar(dimless, 1),
        fvPatchFieldBase::zeroGradientType()
    ),
    rho0_(materialProperties_.get<scalar>("rho")),
    noTension_
    (
        (materialProperties_.found("c") && materialProperties_.get<scalar>("c") == 0)
     || (
            materialProperties_.found("sigmaT")
         && materialProperties_.get<scalar>("sigmaT") == 0
        )
    )
{
    if (maxCo_ <= 0)
    {
        FatalIOErrorInFunction(mesh.time().controlDict())
            << "maxCo must be positive, got " << maxCo_
            << exit(FatalIOError);
    }

    if (freeSurface_)
    {
        setCartesianGeometry();
        lean_.setSize(mesh.nCells(), Zero);
        alphaL_ = alpha_.primitiveField();
        smooth_.setSize(mesh.nCells(), Zero);
        Info<< "Free surface: material plus void, alpha in ["
            << gMin(alpha_.primitiveField()) << ", "
            << gMax(alpha_.primitiveField()) << "], "
            << (noTension_ ? "no tensile strength" : "tensile strength")
            << nl;
    }
    else
    {
        alpha_.writeOpt(IOobject::NO_WRITE);
    }

    const scalar minRho = gMin(rho_.primitiveField());
    if (minRho < 0 || (minRho == 0 && !freeSurface_))
    {
        FatalErrorInFunction
            << "rho must be positive in every cell, smallest is " << minRho
            << exit(FatalError);
    }

    Info<< "Eulerian solid: rho in [" << minRho << ", "
        << gMax(rho_.primitiveField()) << "] kg/m^3, maxCo = " << maxCo_
        << ", solveMomentum " << solveMomentum_
        << ", secondOrder " << secondOrder_ << nl;
    if (mag(g_) > 0)
    {
        Info<< "Gravity " << g_ << " m/s^2, ramp " << gravityRamp_ << " s"
            << nl;
    }
    if (damping_ > 0)
    {
        Info<< "Velocity damping " << damping_ << " 1/s" << nl;
    }
    forAllConstIters(wallFriction_, iter)
    {
        Info<< "Wall " << iter.key() << ": friction coefficient "
            << iter.val() << nl;
    }
    Info<< endl;
}


Foam::vector Foam::eulerianSolid::gravityAt(const scalar t) const
{
    if (gravityRamp_ <= 0 || t >= gravityRamp_)
    {
        return g_;
    }
    return
        0.5*(1 - Foam::cos(constant::mathematical::pi*max(t, 0)/gravityRamp_))
       *g_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eulerianSolid::setCartesianGeometry()
{
    // h_a = V / A_a with A_a half the summed |S_a| of the cell's faces;
    // every face must be normal to an axis (decision 0031)
    const vectorField& Sf = mesh_.Sf().primitiveField();
    const scalarField& magSf = mesh_.magSf().primitiveField();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    vectorField A(mesh_.nCells(), Zero);

    auto check = [&](const vector& S, const scalar mS)
    {
        if (cmptMax(cmptMag(S)) < (1 - 1e-6)*mS)
        {
            FatalErrorInFunction
                << "freeSurface needs a Cartesian hex mesh: face normal "
                << S/mS << " is not along an axis" << exit(FatalError);
        }
    };
    forAll(own, facei)
    {
        check(Sf[facei], magSf[facei]);
        A[own[facei]] += 0.5*cmptMag(Sf[facei]);
        A[nei[facei]] += 0.5*cmptMag(Sf[facei]);
    }
    // Poly patches: empty patches (2D front and back) count here too
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& p = mesh_.boundaryMesh()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField pS(p.faceAreas());
        forAll(fc, i)
        {
            check(pS[i], mag(pS[i]));
            A[fc[i]] += 0.5*cmptMag(pS[i]);
        }
    }
    const scalarField& V = mesh_.V().field();
    hCell_.setSize(mesh_.nCells());
    forAll(hCell_, celli)
    {
        for (direction a = 0; a < vector::nComponents; ++a)
        {
            hCell_[celli][a] = V[celli]/A[celli][a];
        }
    }
}


Foam::scalar Foam::eulerianSolid::stableDeltaT() const
{
    if (freeSurface_)
    {
        // Wave speed of the material (density rho/alpha, bounded as in the
        // Lagrange step); empty cells take the reference density
        const scalarField& r = rho_.primitiveField();
        const scalarField& a = alpha_.primitiveField();
        const scalarField delta(cbrt(mesh_.V().field()));
        scalar dt = GREAT;
        forAll(r, celli)
        {
            const scalar rm =
                a[celli] > 1e-12
              ? min(max(r[celli]/a[celli], 0.5*rho0_), 2*rho0_)
              : rho0_;
            dt = min(dt, delta[celli]/sqrt(material_->waveModulus()/rm));
        }
        return maxCo_*returnReduce(dt, minOp<scalar>());
    }

    // Cell size cbrt(V) over P-wave speed sqrt(M/rho), smallest over the mesh
    const scalarField delta(cbrt(mesh_.V().field()));
    const scalarField c(sqrt(material_->waveModulus()/rho_.primitiveField()));

    return maxCo_*gMin(delta/c);
}


// remap(dt): eulerianSolidRemap.C (WP 2.4)


Foam::scalar Foam::eulerianSolid::totalMass() const
{
    return gSum(rho_.primitiveField()*mesh_.V().field());
}


Foam::vector Foam::eulerianSolid::totalMomentum() const
{
    return gSum
    (
        rho_.primitiveField()*mesh_.V().field()*U_.primitiveField()
    );
}


// ************************************************************************* //
