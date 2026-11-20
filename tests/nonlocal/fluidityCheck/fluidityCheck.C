/*---------------------------------------------------------------------------*\
Application
    celFluidityCheck

Description
    WP N.1 test driver (decision 0111): the nonlocal fluidity of a prescribed
    stress field. system/fluidityCheckDict gives the law (phi, mu2, I0,
    nonlocalA, d, rhoS) and linear profiles in y of the pressure and the
    friction, P = P0 + P1 y, mu = mu0 + mu1 y (sigma = -P I + mu P (xy));
    optional emptyAbove: cells above that height are void. Writes
    fluidity.dat: y, g per cell.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nonlocalFluidity.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "fluidityCheckDict", runTime.system(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );
    cel::constitutive::Parameters p;
    for (const word k : {"phi", "mu2", "I0", "nonlocalA", "d", "rhoS"})
    {
        p[k] = dict.get<scalar>(k);
    }
    auto law = cel::constitutive::Nonlocal::New(p, p["phi"]);

    const scalar P0 = dict.get<scalar>("P0"), P1 = dict.get<scalar>("P1");
    const scalar m0 = dict.get<scalar>("mu0"), m1 = dict.get<scalar>("mu1");
    const scalar yEmpty = dict.getOrDefault<scalar>("emptyAbove", GREAT);

    volSymmTensorField sigma
    (
        IOobject("sigma", runTime.timeName(), mesh),
        mesh,
        dimensionedSymmTensor(dimPressure, Zero)
    );
    scalarField a(mesh.nCells(), 1);
    const volVectorField& C = mesh.C();
    forAll(sigma, i)
    {
        const scalar y = C[i].y();
        const scalar P = P0 + P1*y;
        sigma[i] = symmTensor(-P, (m0 + m1*y)*P, 0, -P, 0, -P);
        if (y > yEmpty)
        {
            a[i] = 0;
            sigma[i] = Zero;
        }
    }
    wordList types(mesh.boundary().size(), "zeroGradient");
    forAll(mesh.boundary(), patchi)
    {
        if (polyPatch::constraintType(mesh.boundary()[patchi].type()))
        {
            types[patchi] = mesh.boundary()[patchi].type();
        }
    }
    volScalarField g
    (
        IOobject("g", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimless/dimTime, Zero),
        types
    );
    nonlocalFluidity(mesh, *law, sigma, a, g);

    OFstream os("fluidity.dat");
    forAll(g, i)
    {
        os << C[i].y() << ' ' << g[i] << nl;
    }
    Info<< "fluidity.dat written, max g " << gMax(g) << nl << "End" << nl;
    return 0;
}
