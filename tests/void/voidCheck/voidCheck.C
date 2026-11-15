/*---------------------------------------------------------------------------*\
Application
    celVoidCheck

Description
    WP 3.2a test driver: one material plus void (materialProperties
    freeSurface yes) run with the celFoam cycle (lagrangeStep then remap)
    at a fixed time step. Run by tests/test_void.sh through
    tests/void/void_checks.py.

    system/voidCheckDict:

    \verbatim
        fill        (x0 y0 x1 y1);  // material rectangle (cell fractions)
        U0          (0 0 0);        // initial velocity of the material
        xSplit      0.3;            // optional: material left of xSplit
        Uleft       (-0.5 0 0);     //   starts with Uleft instead
        initK0      0.4286;         // optional: K0 stress of the weight
                                    //   above (g along -y, top at y1)
        tEnd        0.1;            // run to tEnd with dt = tEnd/ceil(...)
        nSteps      1000;           // or this many steps of the stable dt
    \endverbatim

    Output (master unless noted):
        history.dat  step t mass px py cx cy max|s| max|s|_empty max|U|
                     Fy_bottom max(s_xx, s_yy)
                     (s: material stress; Fy_bottom: force of the patch
                     bottom on the material)
        cells.dat    (cells.dat.procN in parallel) at the end:
                     x y alpha rho Ux Uy s_xx s_yy s_xy
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "OFstream.H"
#include "boundBox.H"
#include <cmath>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("WP 3.2a: material plus void checks");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    eulerianSolid solid(mesh);
    if (!solid.freeSurface())
    {
        FatalErrorInFunction << "needs freeSurface yes" << exit(FatalError);
    }

    IOdictionary dict
    (
        IOobject
        (
            "voidCheckDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalarField& V = mesh.V().field();
    const vectorField& C = mesh.C().primitiveField();
    const scalar rho0 = solid.rho0();
    const scalar gMag = mag(solid.g());

    // Initial state
    const scalarList box(dict.get<scalarList>("fill"));
    const vector U0(dict.get<vector>("U0"));
    scalarField& alpha = solid.alphaRef().primitiveFieldRef();
    scalarField& rho = solid.rhoRef().primitiveFieldRef();
    vectorField& U = solid.URef().primitiveFieldRef();
    symmTensorField& s = solid.sigmaRef().primitiveFieldRef();
    const labelListList& cp = mesh.cellPoints();
    forAll(alpha, celli)
    {
        const boundBox bb(mesh.points(), cp[celli], false);
        const scalar hx = bb.max().x() - bb.min().x();
        const scalar hy = bb.max().y() - bb.min().y();
        const scalar fx =
            min(max((min(bb.max().x(), box[2]) - max(bb.min().x(), box[0]))/hx, 0), 1);
        const scalar fy =
            min(max((min(bb.max().y(), box[3]) - max(bb.min().y(), box[1]))/hy, 0), 1);
        scalar a = fx*fy;
        a = a >= 1 - 1e-12 ? 1 : (a > 1e-12 ? a : 0);
        alpha[celli] = a;
        rho[celli] = rho0*a;
        U[celli] = a > 0 ? U0 : Zero;
        if (a > 0 && dict.found("xSplit") && C[celli].x() < dict.get<scalar>("xSplit"))
        {
            U[celli] = dict.get<vector>("Uleft");
        }
        s[celli] = Zero;
        if (a > 0 && dict.found("initK0"))
        {
            // weight above the centre per area over alpha (material stress)
            const scalar K0 = dict.get<scalar>("initK0");
            const scalar above =
                rho0*max(box[3] - bb.max().y(), 0) + 0.5*rho[celli]*hy;
            const scalar syy = -gMag*above/a;
            s[celli] = symmTensor(K0*syy, 0, 0, syy, 0, K0*syy);
        }
    }
    solid.alphaRef().correctBoundaryConditions();
    solid.rhoRef().correctBoundaryConditions();
    solid.URef().correctBoundaryConditions();
    solid.sigmaRef().correctBoundaryConditions();

    const scalar dt0 = solid.stableDeltaT();
    label nSteps;
    scalar dt;
    if (dict.found("tEnd"))
    {
        const scalar tEnd = dict.get<scalar>("tEnd");
        nSteps = label(std::ceil(tEnd/dt0 - 1e-9));
        dt = tEnd/nSteps;
    }
    else
    {
        nSteps = dict.get<label>("nSteps");
        dt = dt0;
    }
    Info<< "dt " << dt << " s, " << nSteps << " steps" << nl << endl;

    const label bottom = mesh.boundary().findPatchID("bottom");

    autoPtr<OFstream> hist;
    if (Pstream::master())
    {
        hist.reset(new OFstream("history.dat"));
        hist().precision(17);
        hist() << "# step t mass px py cx cy smax smaxEmpty Umax FyBottom snMax"
            << nl;
    }

    auto record = [&](const label step)
    {
        const scalarField& r = solid.rho().primitiveField();
        const scalarField& a = solid.alpha().primitiveField();
        const vectorField& u = solid.U().primitiveField();
        const symmTensorField& sg = solid.sigma().primitiveField();
        const scalarField m(r*V);
        const scalar mass = gSum(m);
        const vector p = gSum(m*u);
        const vector mc = gSum(m*C);
        scalar smax = 0, sEmpty = 0, umax = 0, snMax = -GREAT;
        forAll(m, celli)
        {
            const scalar sm = cmptMax(cmptMag(sg[celli]));
            smax = max(smax, sm);
            if (a[celli] <= 1e-12)
            {
                sEmpty = max(sEmpty, sm);
            }
            umax = max(umax, mag(u[celli]));
            snMax = max(snMax, max(sg[celli].xx(), sg[celli].yy()));
        }
        reduce(smax, maxOp<scalar>());
        reduce(sEmpty, maxOp<scalar>());
        reduce(umax, maxOp<scalar>());
        reduce(snMax, maxOp<scalar>());
        scalar fy = 0;
        if (bottom >= 0)
        {
            fy = sum(solid.Fstar().boundaryField()[bottom].component(1));
        }
        reduce(fy, sumOp<scalar>());
        if (Pstream::master())
        {
            hist() << step << ' ' << runTime.value() << ' ' << mass << ' '
                << p.x() << ' ' << p.y() << ' ' << mc.x()/mass << ' '
                << mc.y()/mass << ' ' << smax << ' ' << sEmpty << ' ' << umax
                << ' ' << fy << ' ' << snMax << nl;
        }
    };

    record(0);
    for (label step = 1; step <= nSteps; ++step)
    {
        runTime.setDeltaT(dt);
        ++runTime;
        solid.lagrangeStep(dt);
        solid.remap(dt);
        record(step);
    }

    {
        const fileName name
        (
            Pstream::parRun()
          ? fileName("cells.dat.proc" + Foam::name(Pstream::myProcNo()))
          : fileName("cells.dat")
        );
        OFstream os(runTime.globalPath()/name);
        os.precision(17);
        os << "# x y alpha rho Ux Uy sxx syy sxy" << nl;
        const scalarField& r = solid.rho().primitiveField();
        const scalarField& a = solid.alpha().primitiveField();
        const vectorField& u = solid.U().primitiveField();
        const symmTensorField& sg = solid.sigma().primitiveField();
        forAll(r, celli)
        {
            os << C[celli].x() << ' ' << C[celli].y() << ' ' << a[celli] << ' '
                << r[celli] << ' ' << u[celli].x() << ' ' << u[celli].y()
                << ' ' << sg[celli].xx() << ' ' << sg[celli].yy() << ' '
                << sg[celli].xy() << nl;
        }
    }

    Info<< "End" << nl << endl;
    return 0;
}
