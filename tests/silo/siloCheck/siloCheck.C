/*---------------------------------------------------------------------------*\
Application
    celSiloCheck

Description
    WP 3.1 test driver: a tank under gravity with wall patches, run with the
    celFoam cycle (stable time step, lagrangeStep then remap). Run by
    tests/test_silo_settle.sh through tests/silo/silo_checks.py.

    system/siloCheckDict:

    \verbatim
        mode        settle;     // settle | fixed
        nSteps      1000;       // fixed: number of cycles
        keTol       1e-6;       // settle: KE <= keTol * strain energy,
        checkEvery  20;         //   checked every checkEvery steps after
        tMax        1.0;        //   the gravity ramp, up to tMax
        continueFraction 0.1;   // settle: then this fraction of the steps
                                //   again with damping off
        initK0      0.4286;     // optional: start in discrete equilibrium
                                //   under full gravity (uniform rho, g
                                //   along -y, free top at the highest y):
                                //   s_yy from the weight above, s_xx = s_zz
                                //   = initK0 s_yy
    \endverbatim

    Output (master unless noted):
        history.dat   step t KE strainEnergy max|U| mass
        summary.dat   key value lines (steps, settled, forces, ...)
        cells.dat     (cells.dat.procN in parallel) at the end of settling or
                      of the fixed run: x y Ux Uy s_xx s_yy s_zz s_xy rho epsPEq
        walls.dat     (walls.dat.procN) wall faces at the same time:
                      patch x y tx ty vx vy slip
                      (traction = face force / area on the material)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "rigidBodies.H"
#include "bodyCells.H"
#include "immersedWalls.H"
#include "OFstream.H"
#include <cmath>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("WP 3.1: tank settling under gravity with walls");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    eulerianSolid solid(mesh);

    // STL walls (WP 4.3a): fixed rigid bodies from constant/rigidBodies,
    // if present, as cut-cell walls
    rigidBodies bodies(runTime);
    autoPtr<bodyCells> cut;
    if (bodies.size())
    {
        cut.reset(new bodyCells(mesh, bodies));
        setImmersedWalls(solid, cut.ref(), bodies);
    }
    const scalarField Vopen(solid.openVolume());

    IOdictionary dict
    (
        IOobject
        (
            "siloCheckDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word mode(dict.get<word>("mode"));
    const vectorField& C = mesh.C().primitiveField();
    const scalar G = solid.material().shearModulus();
    const scalar K = solid.material().waveModulus() - 4.0*G/3.0;
    const scalar gMag = mag(solid.g());

    if (dict.found("initK0"))
    {
        const scalar K0 = dict.get<scalar>("initK0");
        const scalar yTop = gMax(mesh.points().component(vector::Y)());
        const scalarField& rho = solid.rho().primitiveField();
        symmTensorField& s = solid.sigmaRef().primitiveFieldRef();
        forAll(s, celli)
        {
            const scalar syy = -rho[celli]*gMag*(yTop - C[celli].y());
            s[celli] = symmTensor(K0*syy, 0, 0, syy, 0, K0*syy);
        }
        solid.sigmaRef().correctBoundaryConditions();
    }

    labelList wallPatches;
    forAll(mesh.boundary(), patchi)
    {
        if (solid.wallFriction().found(mesh.boundary()[patchi].name()))
        {
            wallPatches.append(patchi);
        }
    }
    List<scalarField> slip(mesh.boundary().size());
    forAll(slip, patchi)
    {
        slip[patchi].setSize(mesh.boundary()[patchi].size(), 0);
    }

    const scalar mass0 = solid.totalMass();

    auto kinetic = [&]()
    {
        return gSum
        (
            0.5*solid.rho().primitiveField()*Vopen
           *magSqr(solid.U().primitiveField())
        );
    };
    auto strainEnergy = [&]()
    {
        const symmTensorField& s = solid.sigma().primitiveField();
        scalar u = 0;
        forAll(s, celli)
        {
            const symmTensor d(dev(s[celli]));
            u += Vopen[celli]
               *((d && d)/(4*G) + sqr(tr(s[celli]))/(18*K));
        }
        return returnReduce(u, sumOp<scalar>());
    };
    auto maxSpeed = [&]()
    {
        const vectorField& U = solid.U().primitiveField();
        scalar v = 0;
        forAll(U, celli)
        {
            if (Vopen[celli] > 0) v = max(v, mag(U[celli]));
        }
        return returnReduce(v, maxOp<scalar>());
    };

    autoPtr<OFstream> hist;
    if (Pstream::master())
    {
        hist.reset(new OFstream("history.dat"));
        hist().precision(12);
        hist() << "# step t KE strainEnergy max|U| mass" << nl;
    }

    label step = 0;
    scalar dtLast = 0;

    // One celFoam cycle
    auto cycle = [&]()
    {
        const scalar dt = solid.stableDeltaT();
        runTime.setDeltaT(dt);
        ++runTime;
        solid.lagrangeStep(dt);
        solid.remap(dt);
        ++step;
        dtLast = dt;
        for (const label patchi : wallPatches)
        {
            const vectorField& v = solid.vStar().boundaryField()[patchi];
            const vectorField nf(mesh.boundary()[patchi].nf());
            forAll(v, i)
            {
                slip[patchi][i] += dt*mag(v[i] - (v[i] & nf[i])*nf[i]);
            }
        }
    };

    auto record = [&](scalar& ke, scalar& se)
    {
        ke = kinetic();
        se = strainEnergy();
        const scalar vmax = maxSpeed();
        const scalar mass = solid.totalMass();
        if (Pstream::master())
        {
            hist() << step << ' ' << runTime.value() << ' ' << ke << ' '
                << se << ' ' << vmax << ' ' << mass << nl;
        }
    };

    auto procName = [&](const word& base)
    {
        return runTime.globalPath()/
        (
            Pstream::parRun()
          ? fileName(base + ".proc" + Foam::name(Pstream::myProcNo()))
          : fileName(base)
        );
    };

    auto writeState = [&]()
    {
        {
            OFstream os(procName("cells.dat"));
            os.precision(17);
            os  << "# x y Ux Uy s_xx s_yy s_zz s_xy rho epsPEq open" << nl;
            const scalarField open(solid.cellOpen());
            forAll(C, celli)
            {
                const vector& u = solid.U()[celli];
                const symmTensor& s = solid.sigma()[celli];
                os  << C[celli].x() << ' ' << C[celli].y() << ' '
                    << u.x() << ' ' << u.y() << ' '
                    << s.xx() << ' ' << s.yy() << ' ' << s.zz() << ' '
                    << s.xy() << ' ' << solid.rho()[celli] << ' '
                    << solid.epsPEq()[celli] << ' ' << open[celli] << nl;
            }
        }
        OFstream os(procName("walls.dat"));
        os.precision(17);
        os  << "# patch x y tx ty vx vy slip" << nl;
        for (const label patchi : wallPatches)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const vectorField& Cf = p.Cf();
            const scalarField& mS = p.magSf();
            const vectorField& F = solid.Fstar().boundaryField()[patchi];
            const vectorField& v = solid.vStar().boundaryField()[patchi];
            forAll(Cf, i)
            {
                const vector t(F[i]/mS[i]);
                os  << p.name() << ' ' << Cf[i].x() << ' ' << Cf[i].y() << ' '
                    << t.x() << ' ' << t.y() << ' ' << v[i].x() << ' '
                    << v[i].y() << ' ' << slip[patchi][i] << nl;
            }
        }
    };

    // Upward force on the material from each wall patch, per unit depth
    auto patchForces = [&](autoPtr<OFstream>& os, const word& prefix)
    {
        const scalar depth =
            gMax(mesh.points().component(vector::Z)())
          - gMin(mesh.points().component(vector::Z)());
        for (const label patchi : wallPatches)
        {
            const scalar fy =
                gSum(solid.Fstar().boundaryField()[patchi].component(1)());
            if (Pstream::master())
            {
                os() << prefix << "force_" << mesh.boundary()[patchi].name()
                    << ' ' << -fy*sign(solid.g().y())/depth << nl;
            }
        }
        const scalar weight = solid.totalMass()*gMag/depth;
        if (Pstream::master())
        {
            os() << prefix << "weight " << weight << nl;
        }
        // STL walls: force of each body on the material per unit depth,
        // and the mass that entered through its wall faces (leak report)
        forAll(solid.bodyNames(), b)
        {
            const vector f(solid.bodyForce(b)/depth);
            const scalar mIn = solid.bodyMassIn(b);
            if (Pstream::master())
            {
                const word& n = solid.bodyNames()[b];
                os() << prefix << "bforce_" << n << "_x " << f.x() << nl
                    << prefix << "bforce_" << n << "_y " << f.y() << nl
                    << prefix << "massIn_" << n << ' ' << mIn << nl;
            }
        }
        const scalar mass = solid.totalMass();
        if (Pstream::master())
        {
            os() << prefix << "mass " << mass << nl
                << prefix << "gx " << solid.g().x() << nl
                << prefix << "gy " << solid.g().y() << nl;
        }
    };

    autoPtr<OFstream> sum;
    if (Pstream::master())
    {
        sum.reset(new OFstream("summary.dat"));
        sum().precision(12);
    }

    scalar ke = 0, se = 0;

    if (mode == "fixed")
    {
        const label nSteps = dict.get<label>("nSteps");
        scalar vmax = 0;
        for (label k = 0; k < nSteps; ++k)
        {
            cycle();
            vmax = max(vmax, maxSpeed());
        }
        record(ke, se);
        writeState();
        if (Pstream::master())
        {
            sum() << "steps " << step << nl << "t " << runTime.value() << nl
                << "dt " << dtLast << nl << "vmax " << vmax << nl;
        }
        patchForces(sum, "");
    }
    else if (mode == "settle")
    {
        const scalar keTol = dict.get<scalar>("keTol");
        const label every = dict.getOrDefault<label>("checkEvery", 20);
        const scalar tMax = dict.get<scalar>("tMax");
        const scalar frac = dict.getOrDefault<scalar>("continueFraction", 0.1);
        const scalar tRamp =
            runTime.controlDict().getOrDefault<scalar>("gravityRamp", 0);

        bool settled = false;
        while (runTime.value() < tMax)
        {
            cycle();
            if (step % every == 0)
            {
                record(ke, se);
                if (runTime.value() >= tRamp && se > 0 && ke <= keTol*se)
                {
                    settled = true;
                    break;
                }
            }
        }
        writeState();
        const scalar massNow = solid.totalMass();
        if (Pstream::master())
        {
            sum() << "settled " << label(settled) << nl
                << "steps " << step << nl << "t " << runTime.value() << nl
                << "dt " << dtLast << nl
                << "massChange " << (massNow - mass0)/mass0 << nl;
        }
        patchForces(sum, "");

        // Continue without damping
        const symmTensorField s0(solid.sigma().primitiveField());
        const label nMore = label(std::ceil(frac*step));
        solid.setDamping(0);
        scalar vmax = 0;
        for (label k = 0; k < nMore; ++k)
        {
            cycle();
            vmax = max(vmax, maxSpeed());
        }
        record(ke, se);
        const scalar ds =
            gMax(mag(solid.sigma().primitiveField() - s0)());
        if (Pstream::master())
        {
            sum() << "contSteps " << nMore << nl << "contVmax " << vmax << nl
                << "contMaxStressChange " << ds << nl
                << "contKEratio " << ke/se << nl;
        }
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "mode settle or fixed" << exit(FatalIOError);
    }

    Info<< "End" << nl << endl;
    return 0;
}
