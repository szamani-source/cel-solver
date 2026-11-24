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
        fillAngle   -30;            // optional (WP G.3): fill and the
                                    //   cell bounds in a frame turned by
                                    //   this angle (deg) about z, for a
                                    //   mesh of blocks turned the same way
        fillRight   ((y0 x0) (y1 x1) ...);
                                    // optional: material only left of
                                    //   x(y) (linear in y, cell centre
                                    //   height), e.g. inside a silo wall
        fillSlope   30;             // optional (WP G.3): fill and cell
                                    //   bounds in (x, y + x tan(slope)),
                                    //   for a strip sheared along a slope
        U0          (0 0 0);        // initial velocity of the material
        xSplit      0.3;            // optional: material left of xSplit
        Uleft       (-0.5 0 0);     //   starts with Uleft instead
        sigma0      (-1e5 0 0 -1e5 0 -1e5);
                                    // optional (WP 4.4b): uniform initial
                                    //   stress (xx xy xz yy yz zz)
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

    STL walls (WP 4.3c): bodies in constant/rigidBodies, if present,
    are cut-cell walls (decision 0047, 0051); moving bodies push the
    material (WP 4.3b, decision 0053). alpha is then the material
    fraction of each cell's open volume; mass and momentum are summed over
    open volumes, and history.dat gets three more columns: the force of the
    bodies on the material (x, y) and the mass in cells whose centre lies
    below y = belowY (voidCheckDict, default -GREAT: none); then (WP 4.3b)
    the force of the mesh walls on the material (x, y: all non-coupled
    patches), the mass covered cells could not hand over, and the largest
    material density rho/alpha over rho0 (cells with alpha > 1e-3); then
    (WP 3.7a) the mass fed in through inlets.

    Output (master unless noted):
        history.dat  step t mass px py cx cy max|s| max|s|_empty max|U|
                     Fy_bottom max(s_xx, s_yy) KE nCapped discharged
                     (s: material stress; Fy_bottom: force of the patch
                     bottom on the material; discharged: mass that
                     left through outlets, WP 3.3)
        cells.dat    (cells.dat.procN in parallel) at the end:
                     x y alpha rho Ux Uy s_xx s_yy s_xy s_zz open
                     [pc, with the compaction cap, WP K.4]
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "rigidBodies.H"
#include "bodyCells.H"
#include "immersedWalls.H"
#include "OFstream.H"
#include "boundBox.H"
#include "IFstream.H"
#include "IStringStream.H"
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

    // Restart (decision 0117): started from a written time (controlDict
    // startFrom latestTime), the fields come from that time; no fill
    const bool restart = runTime.value() > 0;
    if (restart)
    {
        Info<< "Restart from t = " << runTime.timeName() << nl;
    }
    else
    {
        // Initial state
        const scalarList box(dict.get<scalarList>("fill"));
        const vector U0(dict.get<vector>("U0"));
        scalarField& alpha = solid.alphaRef().primitiveFieldRef();
        scalarField& rho = solid.rhoRef().primitiveFieldRef();
        vectorField& U = solid.URef().primitiveFieldRef();
        symmTensorField& s = solid.sigmaRef().primitiveFieldRef();
        const labelListList& cp = mesh.cellPoints();
        // Points in the fill frame (WP G.3): x' = x cos a + y sin a,
        // y' = -x sin a + y cos a, a = fillAngle
        pointField fillPts(mesh.points());
        if (dict.found("fillAngle"))
        {
            const scalar fa =
                dict.get<scalar>("fillAngle")*constant::mathematical::pi/180;
            const scalar ca = Foam::cos(fa), sa = Foam::sin(fa);
            forAll(fillPts, pi)
            {
                const point& q = mesh.points()[pi];
                fillPts[pi] =
                    point(q.x()*ca + q.y()*sa, -q.x()*sa + q.y()*ca, q.z());
            }
        }
        if (dict.found("fillSlope"))
        {
            const scalar ts = Foam::tan
            (
                dict.get<scalar>("fillSlope")*constant::mathematical::pi/180
            );
            forAll(fillPts, pi)
            {
                fillPts[pi].y() += fillPts[pi].x()*ts;
            }
        }
        forAll(alpha, celli)
        {
            const boundBox bb(fillPts, cp[celli], false);
            const scalar hx = bb.max().x() - bb.min().x();
            const scalar hy = bb.max().y() - bb.min().y();
            scalar fx =
                min(max((min(bb.max().x(), box[2]) - max(bb.min().x(), box[0]))/hx, 0), 1);
            const scalar fy =
                min(max((min(bb.max().y(), box[3]) - max(bb.min().y(), box[1]))/hy, 0), 1);
            if (dict.found("fillRight"))
            {
                const List<vector2D> xr(dict.get<List<vector2D>>("fillRight"));
                const scalar yc = 0.5*(bb.min().y() + bb.max().y());
                scalar xl = xr.last().y();
                if (yc <= xr.first().x())
                {
                    xl = xr.first().y();
                }
                for (label k = 1; k < xr.size(); ++k)
                {
                    if (yc >= xr[k-1].x() && yc <= xr[k].x())
                    {
                        const scalar w =
                            (yc - xr[k-1].x())/max(xr[k].x() - xr[k-1].x(), SMALL);
                        xl = xr[k-1].y() + w*(xr[k].y() - xr[k-1].y());
                        break;
                    }
                }
                fx *= min(max((xl - bb.min().x())/hx, 0), 1);
            }
            scalar a = fx*fy;
            // round-off of the fill frame (fillAngle, fillSlope: turned points)
            // leaves slivers of 1e-9 next to the fill; treat them as empty/full
            const scalar tol = (dict.found("fillAngle") || dict.found("fillSlope")) ? 1e-6 : 1e-12;
            a = a >= 1 - tol ? 1 : (a > tol ? a : 0);
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
            if (a > 0 && dict.found("sigma0"))
            {
                // uniform initial material stress (WP 4.4b): a confining
                // pressure in weightless runs
                s[celli] = dict.get<symmTensor>("sigma0");
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
    }

    // STL walls (WP 4.3c): fixed bodies as cut-cell walls
    rigidBodies bodies(runTime);
    autoPtr<bodyCells> cut;
    if (bodies.size())
    {
        cut.reset(new bodyCells(mesh, bodies));
        setImmersedWalls(solid, cut.ref(), bodies);
    }
    // Initial heights of a deformable first body's nodes (WP 5.3)
    scalarField y0Body;
    scalar y0Top = 0, y0Bot = 0;
    if (bodies.size() && bodies[0].deformable())
    {
        y0Body = bodies[0].hex().points().component(1);
        y0Top = max(y0Body);
        y0Bot = min(y0Body);
    }

    auto cycle = [&](const scalar dtk)
    {
        solid.lagrangeStep(dtk);
        moveBodies(solid, bodies);
        solid.remap(dtk);
        if (cut)
        {
            moveImmersedWalls(solid, cut.ref(), bodies);
        }
    };
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
    label restartStep = 0;
    if (Pstream::master())
    {
        // on a restart keep the lines up to the restart time and append
        DynamicList<string> keep;
        label step0 = 0;
        if (restart && isFile("history.dat"))
        {
            IFstream is("history.dat");
            string line;
            while (is.getLine(line) && is.good())
            {
                if (line.empty()) continue;
                if (line[0] != '#')
                {
                    IStringStream ls(line);
                    label st;
                    scalar tl;
                    ls >> st >> tl;
                    if (tl > runTime.value() + 1e-12) break;
                    step0 = st;
                }
                keep.append(line);
            }
        }
        restartStep = step0;
        hist.reset(new OFstream("history.dat"));
        hist().precision(17);
        for (const string& l : keep)
        {
            hist() << l.c_str() << nl;
        }
        if (keep.empty())
        hist() << "# step t mass px py cx cy smax smaxEmpty Umax FyBottom snMax"
            << " KE nCapped discharged FxBodies FyBodies massBelow"
            << " FxWalls FyWalls lostMass maxRhoRatio charged"
            << " body0: cx cy vx vy omegaz (deformable: height)" << nl;
    }

    auto record = [&](const label step)
    {
        const scalarField& r = solid.rho().primitiveField();
        const scalarField& a = solid.alpha().primitiveField();
        const vectorField& u = solid.U().primitiveField();
        const symmTensorField& sg = solid.sigma().primitiveField();
        const scalarField Vopen(solid.openVolume());
        const scalarField m(r*Vopen);
        const scalar mass = gSum(m);
        const vector p = gSum(m*u);
        const vector mc = gSum(m*C);
        const scalar ke = 0.5*gSum(m*magSqr(u));
        const label nCap = returnReduce(solid.nCapped(), sumOp<label>());
        const scalar out = solid.discharged();
        const scalar fed = solid.charged();
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
        vector fw(Zero);
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled()) continue;
            fw += sum(solid.Fstar().boundaryField()[patchi]);
        }
        reduce(fw, sumOp<vector>());
        const scalar lost = solid.lostMass();
        scalar rhoMax = 0;
        forAll(m, celli)
        {
            if (a[celli] > 1e-3)
            {
                rhoMax = max(rhoMax, r[celli]/a[celli]);
            }
        }
        reduce(rhoMax, maxOp<scalar>());
        if (Pstream::master())
        {
            hist() << step << ' ' << runTime.value() << ' ' << mass << ' '
                << p.x() << ' ' << p.y() << ' '
                << (mass > 0 ? mc.x()/mass : 0) << ' '
                << (mass > 0 ? mc.y()/mass : 0) << ' ' << smax << ' ' << sEmpty << ' ' << umax
                << ' ' << fy << ' ' << snMax << ' ' << ke << ' ' << nCap << ' ' << out
                << ' ' << fb.x() << ' ' << fb.y() << ' ' << mBelow
                << ' ' << fw.x() << ' ' << fw.y() << ' ' << lost << ' '
                << rhoMax/solid.rho0() << ' ' << fed;
            if (bodies.size())
            {
                // First body's pose and velocity (WP 4.3d)
                const rigidBody& b0 = bodies[0];
                if (b0.deformable())
                {
                    // Deformable (WP 5.3): mean node position and velocity,
                    // and the height (mean y of the nodes that started on
                    // top less that of the nodes that started at the
                    // bottom) in place of omega
                    const pointField& xp = b0.hex().points();
                    const vector c(average(xp));
                    const vector u(average(b0.hex().velocity()));
                    scalar yt = 0, yb = 0;
                    label nt = 0, nb = 0;
                    forAll(xp, i)
                    {
                        if (y0Body[i] > y0Top - 1e-9) { yt += xp[i].y(); ++nt; }
                        if (y0Body[i] < y0Bot + 1e-9) { yb += xp[i].y(); ++nb; }
                    }
                    hist() << ' ' << c.x() << ' ' << c.y() << ' ' << u.x()
                        << ' ' << u.y() << ' '
                        << yt/max(nt, label(1)) - yb/max(nb, label(1));
                }
                else
                {
                const vector v0(b0.velocity(runTime.value()));
                hist() << ' ' << b0.centre().x() << ' ' << b0.centre().y()
                    << ' ' << v0.x() << ' ' << v0.y() << ' '
                    << b0.omega(runTime.value()).z();
                }
            }
            hist() << nl;
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
        // compaction cap (WP K.4): consolidation pressure pc last, if any
        const volScalarField* pc =
            mesh.findObject<volScalarField>("pc");
        os << "# x y alpha rho Ux Uy sxx syy sxy szz open"
            << (pc ? " pc" : "") << nl;
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
                << open[celli];
            if (pc)
            {
                os << ' ' << pc->primitiveField()[celli];
            }
            os << nl;
        }
    };

    if (!dict.found("maxCo"))
    {
        record(0);
    }
    if (dict.found("maxCo"))
    {
        // Variable step (WP 3.2b): the stable dt of each step at maxCo,
        // shortened to land on tEnd and on each of writeTimes
        const scalar tEnd = dict.get<scalar>("tEnd");
        const scalar co = dict.get<scalar>("maxCo");
        const scalarList wt(dict.getOrDefault<scalarList>("writeTimes", scalarList()));
        // Checkpoints (decision 0117): every checkpointEvery seconds the
        // fields and the solid's restart state are written as a time
        // directory, from which the run restarts (startFrom latestTime)
        const scalar ckEvery = dict.getOrDefault<scalar>("checkpointEvery", 0);
        label next = 0;
        while (next < wt.size() && wt[next] <= runTime.value() + 1e-12)
        {
            ++next;
        }
        label nextCk =
            ckEvery > 0 ? label(std::floor(runTime.value()/ckEvery + 1e-9)) + 1 : 0;
        Pstream::broadcast(restartStep);
        label step = restartStep;
        if (!restart)
        {
            record(0);
        }
        while (runTime.value() < tEnd - 1e-12)
        {
            scalar dtk = co/solid.maxCo()*solid.stableDeltaT();
            scalar target = tEnd;
            if (next < wt.size())
            {
                target = min(target, wt[next]);
            }
            if (ckEvery > 0)
            {
                target = min(target, nextCk*ckEvery);
            }
            dtk = min(dtk, target - runTime.value());
            runTime.setDeltaT(dtk);
            ++runTime;
            cycle(dtk);
            record(++step);
            if (next < wt.size() && runTime.value() >= wt[next] - 1e-12)
            {
                writeCells("cells_" + Foam::name(wt[next]));
                ++next;
            }
            if (ckEvery > 0 && runTime.value() >= nextCk*ckEvery - 1e-12)
            {
                runTime.writeNow();
                solid.writeState();
                if (Pstream::master())
                {
                    hist().flush();
                }
                Info<< "Checkpoint t = " << runTime.timeName() << endl;
                ++nextCk;
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
            cycle(dt);
            record(step);
        }
    }
    writeCells("cells");
    Info<< "End" << nl << endl;
    return 0;
}
