/*---------------------------------------------------------------------------*\
    eulerianSolid: see eulerianSolid.H
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "zeroGradientFvPatchFields.H"

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
    maxCo_(mesh.time().controlDict().getOrDefault<scalar>("maxCo", 0.5)),
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
        << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::eulerianSolid::stableDeltaT() const
{
    // Cell size cbrt(V) over P-wave speed sqrt(M/rho), smallest over the mesh
    const scalarField delta(cbrt(mesh_.V().field()));
    const scalarField c(sqrt(material_->waveModulus()/rho_.primitiveField()));

    return maxCo_*gMin(delta/c);
}


void Foam::eulerianSolid::lagrangeStep(const scalar)
{
    // WP 2.3: stress update with the constitutive model, momentum update.
    // The discretisation is chosen in WP 1.3b. Nothing happens yet.
}


void Foam::eulerianSolid::remap(const scalar)
{
    // WP 2.4: remap mass, momentum, stress and history onto the fixed mesh.
    // Nothing happens yet.
}


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
