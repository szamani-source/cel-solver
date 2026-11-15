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
    )
{
    if (maxCo_ <= 0)
    {
        FatalIOErrorInFunction(mesh.time().controlDict())
            << "maxCo must be positive, got " << maxCo_
            << exit(FatalIOError);
    }

    const scalar minRho = gMin(rho_.primitiveField());
    if (minRho <= 0)
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

Foam::scalar Foam::eulerianSolid::stableDeltaT() const
{
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
