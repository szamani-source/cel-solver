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
        maxCo       0.4;            // optional (WP 3.2b): run to tEnd at
                                    //   the stable dt of each step instead
        writeTimes  (0.5 1);        // optional (WP 3.2b): also write
                                    //   cells_<t>.dat at these times
        layer { theta 31; kzz 0.6; Us 1; }
                                    // optional (WP 3.6b): steady inclined
                                    //   layer stress and Bagnold velocity
                                    //   (fill top = surface, floor y = 0)
    \endverbatim

    STL walls (WP 4.3c): fixed bodies in constant/rigidBodies, if present,
    are cut-cell walls (decision 0047, 0049). alpha is then the material
    fraction of each cell's open volume; mass and momentum are summed over
    open volumes, and history.dat gets three more columns: the force of the
    bodies on the material (x, y) and the mass in cells whose centre lies
    below y = belowY (voidCheckDict, default -GREAT: none).

    Output (master unless noted):
        history.dat  step t mass px py cx cy max|s| max|s|_empty max|U|
                     Fy_bottom max(s_xx, s_yy) KE nCapped discharged
                     (s: material stress; Fy_bottom: force of the patch
                     bottom on the material; discharged: mass that
                     left through outlets, WP 3.3)
        cells.dat    (cells.dat.procN in parallel) at the end:
                     x y alpha rho Ux Uy s_xx s_yy s_xy
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "rigidBodies.H"
#include "bodyCells.H"
#include "immersedWalls.H"
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
        if (a > 0 && dict.found("layer"))
        {
            // Steady inclined layer (WP 3.6b): x down the slope, y normal
            // to it, surface at y1, gravity (g sin theta, -g cos theta).
            // Normal stresses equal in plane (s_xx = s_yy = sn), s_zz =
            // kzz sn, shear rho g sin theta (y1 - y); velocity
            // Us (1 - (1 - y/y1)^1.5) along x (Bagnold, cell centre).
            const dictionary& ld = dict.subDict("layer");
            const scalar th =
                ld.get<scalar>("theta")*constant::mathematical::pi/180;
            const scalar depth = max(box[3] - C[celli].y(), 0);
            const scalar sn = -rho0*gMag*Foam::cos(th)*depth;
            const scalar tau = rho0*gMag*Foam::sin(th)*depth;
            s[celli] =
                symmTensor(sn, tau, 0, sn, 0, ld.get<scalar>("kzz")*sn);
            U[celli] = vector
            (
                ld.get<scalar>("Us")
               *(1 - Foam::pow(max(1 - C[celli].y()/box[3], 0), 1.5)),
                0,
                0
            );
        }
    }
    solid.alphaRef().correctBoundaryConditions();
    solid.rhoRef().correctBoundaryConditions();
    solid.URef().correctBoundaryConditions();
    solid.sigmaRef().correctBoundaryConditions();

    // STL walls (WP 4.3c): fixed bodies as cut-cell walls
    rigidBodies bodies(runTime);
    autoPtr<bodyCells> cut;
    if (bodies.size())
    {
        cut.reset(new bodyCells(mesh, bodies));
        setImmersedWalls(solid, cut.ref(), bodies);
    }
    const scalarField Vopen(solid.openVolume());
    const scalar belowY = dict.getOrDefault<scalar>("belowY", -GREAT);

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
            << " KE nCapped discharged FxBodies FyBodies massBelow" << nl;
    }

    auto record = [&](const label step)
    {
        const scalarField& r = solid.rho().primitiveField();
        const scalarField& a = solid.alpha().primitiveField();
        const vectorField& u = solid.U().primitiveField();
        const symmTensorField& sg = solid.sigma().primitiveField();
        const scalarField m(r*Vopen);
        const scalar mass = gSum(m);
        const vector p = gSum(m*u);
        const vector mc = gSum(m*C);
        const scalar ke = 0.5*gSum(m*magSqr(u));
        const label nCap = returnReduce(solid.nCapped(), sumOp<label>());
        const scalar out = solid.discharged();
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
        vector fb(Zero);
        forAll(solid.bodyNames(), b)
        {
            fb += solid.bodyForce(b);
        }
        scalar mBelow = 0;
        forAll(m, celli)
        {
            if (C[celli].y() < belowY) mBelow += m[celli];
        }
        reduce(mBelow, sumOp<scalar>());
        if (Pstream::master())
        {
            hist() << step << ' ' << runTime.value() << ' ' << mass << ' '
                << p.x() << ' ' << p.y() << ' ' << mc.x()/mass << ' '
                << mc.y()/mass << ' ' << smax << ' ' << sEmpty << ' ' << umax
                << ' ' << fy << ' ' << snMax << ' ' << ke << ' ' << nCap << ' ' << out
                << ' ' << fb.x() << ' ' << fb.y() << ' ' << mBelow << nl;
        }
    };

    auto writeCells = [&](const word& stem)
    {
        const fileName name
        (
            Pstream::parRun()
          ? fileName(stem + ".dat.proc" + Foam::name(Pstream::myProcNo()))
          : fileName(stem + ".dat")
        );
        OFstream os(runTime.globalPath()/name);
        os.precision(17);
        os << "# x y alpha rho Ux Uy sxx syy sxy szz open" << nl;
        const scalarField open(solid.cellOpen());
        const scalarField& r = solid.rho().primitiveField();
        const scalarField& a = solid.alpha().primitiveField();
        const vectorField& u = solid.U().primitiveField();
        const symmTensorField& sg = solid.sigma().primitiveField();
        forAll(r, celli)
        {
            os << C[celli].x() << ' ' << C[celli].y() << ' ' << a[celli] << ' '
                << r[celli] << ' ' << u[celli].x() << ' ' << u[celli].y()
                << ' ' << sg[celli].xx() << ' ' << sg[celli].yy() << ' '
                << sg[celli].xy() << ' ' << sg[celli].zz() << ' '
                << open[celli] << nl;
        }
    };

    record(0);
    if (dict.found("maxCo"))
    {
        // Variable step (WP 3.2b): the stable dt of each step at maxCo,
        // shortened to land on tEnd and on each of writeTimes
        const scalar tEnd = dict.get<scalar>("tEnd");
        const scalar co = dict.get<scalar>("maxCo");
        const scalarList wt(dict.getOrDefault<scalarList>("writeTimes", scalarList()));
        label next = 0;
        label step = 0;
        while (runTime.value() < tEnd - 1e-12)
        {
            scalar dtk = co/solid.maxCo()*solid.stableDeltaT();
            scalar target = tEnd;
            if (next < wt.size())
            {
                target = min(target, wt[next]);
            }
            dtk = min(dtk, target - runTime.value());
            runTime.setDeltaT(dtk);
            ++runTime;
            solid.lagrangeStep(dtk);
            solid.remap(dtk);
            record(++step);
            if (next < wt.size() && runTime.value() >= wt[next] - 1e-12)
            {
                writeCells("cells_" + Foam::name(wt[next]));
                ++next;
            }
        }
        Info<< step << " steps to t = " << runTime.value() << " s" << nl;
    }
    else
    {
        for (label step = 1; step <= nSteps; ++step)
        {
            runTime.setDeltaT(dt);
            ++runTime;
            solid.lagrangeStep(dt);
            solid.remap(dt);
            record(step);
        }
    }
    writeCells("cells");
    Info<< "End" << nl << endl;
    return 0;
}
