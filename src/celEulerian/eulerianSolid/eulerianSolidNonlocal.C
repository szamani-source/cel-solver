/*---------------------------------------------------------------------------*\
    Nonlocal granular fluidity (WP N.1, decision 0111; reference
    prototype/constitutive/ngf.py). Each Lagrange step, from the stress
    at its start:

        k g - div(a_f grad g) = k g_loc,   k = |mu - mu_s|/(A d)^2

    mu = sqrt(J2)/P of the material stress, a_f = min of the material
    fractions on either side of a face (no flux into void; walls zero
    gradient). Empty cells: k = 1/(A d)^2 and g_loc = 0, so g = 0 there.
    A floor of 1e-8/(A d)^2 on k keeps a layer at mu = mu_s solvable.
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "nonlocalFluidity.H"
#include "celNonlocal.H"
#include "fvm.H"
#include "laplacianScheme.H"
#include "IStringStream.H"
#include "SolverPerformance.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eulerianSolid::solveFluidity()
{
    nonlocalFluidity
    (
        mesh_, *material_->nonlocal(), sigma_,
        freeSurface_
      ? scalarField(alpha_.primitiveField())
      : scalarField(mesh_.nCells(), 1),
        fluidity_.ref()
    );
}


void Foam::nonlocalFluidity
(
    const fvMesh& mesh,
    const cel::constitutive::Nonlocal& law,
    const volSymmTensorField& sigma,
    const scalarField& a,
    volScalarField& g
)
{

    const scalar Ad2 = sqr(law.A()*law.d());
    const scalar kMin = 1e-8/Ad2;

    volScalarField::Internal kf
    (
        IOobject("gK", mesh.time().timeName(), mesh, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar(dimless/dimArea, Zero)
    );
    volScalarField::Internal rhs
    (
        IOobject("gRhs", mesh.time().timeName(), mesh, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar(dimless/dimArea/dimTime, Zero)
    );
    forAll(kf, celli)
    {
        if (a[celli] <= 1e-6)
        {
            kf[celli] = 1/Ad2;
            continue;
        }
        double s[6];
        constitutiveModel::toVoigt(sigma[celli], s);
        double mu = 0, P = 0;
        law.muP(s, mu, P);
        kf[celli] = max(law.k(mu), kMin);
        rhs[celli] = kf[celli]*law.gLoc(mu, P);
    }

    // Face diffusivity: the smaller material fraction of the two sides
    surfaceScalarField af
    (
        IOobject("gAf", mesh.time().timeName(), mesh, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar(dimless, Zero)
    );
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        forAll(own, facei)
        {
            af[facei] = min(a[own[facei]], a[nei[facei]]);
        }
        volScalarField aw
        (
            IOobject("gA", mesh.time().timeName(), mesh, IOobject::NO_REGISTER),
            mesh,
            dimensionedScalar(dimless, Zero),
            g.boundaryField().types()
        );
        aw.primitiveFieldRef() = a;
        aw.correctBoundaryConditions();
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            scalarField& pf = af.boundaryFieldRef()[patchi];
            const scalarField ai(aw.boundaryField()[patchi].patchInternalField());
            if (p.coupled())
            {
                pf = min(ai, aw.boundaryField()[patchi].patchNeighbourField());
            }
            else
            {
                pf = ai;
            }
        }
    }

    IStringStream schemeIs("Gauss linear corrected");
    tmp<fv::laplacianScheme<scalar, scalar>> lap
    (
        fv::laplacianScheme<scalar, scalar>::New(mesh, schemeIs)
    );
    fvScalarMatrix gEqn
    (
        fvm::Sp(kf, g) - lap.ref().fvmLaplacian(af, g) == rhs
    );
    IStringStream solverIs
    (
        "solver PCG; preconditioner DIC; tolerance 1e-14; relTol 0;"
        " maxIter 10000;"
    );
    const dictionary solverDict(solverIs);
    const int dbg = SolverPerformance<scalar>::debug;
    SolverPerformance<scalar>::debug = 0;
    gEqn.solve(solverDict);
    SolverPerformance<scalar>::debug = dbg;
    g.primitiveFieldRef() = max(g.primitiveField(), scalar(0));
    g.correctBoundaryConditions();
}


// ************************************************************************* //
