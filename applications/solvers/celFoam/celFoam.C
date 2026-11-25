/*---------------------------------------------------------------------------*\
Application
    celFoam

Description
    Explicit Coupled Eulerian-Lagrangian solver for large-deformation solids.

    Thin driver (WP 2.1): the Eulerian state and the steps of the explicit
    cycle live in library celEulerian (class eulerianSolid). Each cycle:
    stable time step, Lagrange step (WP 2.3), remap (WP 2.4), write. The
    two together are checked by tests/test_eulerian_cycle.sh (WP 2.4b);
    rigid bodies from STL move on prescribed paths (WP 4.1); fixed bodies
    are cut-cell walls with friction (WP 4.3a, decision 0047); moving
    bodies are coupled in WP 4.3b (see docs/work-packages.md). Version 2
    (decision 0125): with constant/gasProperties the pore gas (library
    celFluid, class poreGas) follows the solid; F.1 computes the gas
    fractions and permeability after each remap, F.2 solves the pore
    pressure with the solid frozen (decision 0127); F.3 couples the two
    both ways (effective stress, gas carried by the grains, decision 0128).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "rigidBodies.H"
#include "bodyCells.H"
#include "immersedWalls.H"
#include "poreGas.H"
#include <cmath>

// Print total mass and momentum with enough digits to compare runs
static void printDiagnostics(const Foam::eulerianSolid& solid)
{
    using namespace Foam;

    const label oldPrecision = Info().precision(16);

    Info<< "Total mass      = " << solid.totalMass() << nl
        << "Total momentum  = " << solid.totalMomentum() << endl;
    if (solid.freeSurface())
    {
        // WP 3.3 and 3.7a: mass out through outlets, in through inlets
        Info<< "Discharged mass = " << solid.discharged() << nl
            << "Fed mass        = " << solid.charged() << endl;
    }

    Info().precision(oldPrecision);
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Explicit Coupled Eulerian-Lagrangian solver (WP 2.1 driver)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Rigid bodies from constant/rigidBodies (WP 4.1); none if absent
    rigidBodies bodies(runTime);

    // Bodies act on the material as cut-cell walls (WP 4.3a, decision
    // 0047); moving bodies push it and the cut cells follow them each step
    // (WP 4.3b, decision 0053). Moving bodies need freeSurface yes; without
    // it they only move (kinematics, as in WP 4.1)
    autoPtr<bodyCells> cut;
    if (bodies.size() && anyMoving(bodies) && !solid.freeSurface())
    {
        Info<< "Moving rigid bodies without freeSurface yes: not coupled to "
            << "the material (decision 0053)" << nl << endl;
    }
    else if (bodies.size())
    {
        cut.reset(new bodyCells(mesh, bodies));
        setImmersedWalls(solid, cut.ref(), bodies);
    }

    // Pore gas (WP F.1, decision 0126): null without constant/gasProperties
    // or with model none, and then nothing below changes (v1 bit for bit)
    autoPtr<poreGas> gas(poreGas::New(solid));
    if (gas && gas->coupled())
    {
        // WP F.3 (decision 0128): sigma is the effective stress and the gas
        // pressure acts on the material
        // WP F.4 (decision 0129): with a free surface the solid sees the
        // void pressure on the faces between material and void
        if (gas->frictionOnly())
        {
            // WP F.8g (decision 0133): the pore pressure only lowers the
            // friction; sigma stays the total stress
            solid.setSupportPressure(gas->p());
        }
        else if (solid.freeSurface())
        {
            solid.setPorePressure(gas->p(), gas->pFace());
        }
        else
        {
            solid.setPorePressure(gas->p());
        }
        // WP F.7 (decision 0132): the reaction of the pore gas inertia
        if (gas->gasForcePtr())
        {
            solid.setGasForce(*gas->gasForcePtr());
        }
    }

    // Reaction output: force of the material on each body, every step
    autoPtr<OFstream> forceFile;
    if (cut && Pstream::master())
    {
        mkDir(runTime.globalPath()/"postProcessing");
        // A restarted run (WP 6.2) appends to the file of the first run
        const bool restart = runTime.value() > 0;
        forceFile.reset
        (
            new OFstream
            (
                runTime.globalPath()/"postProcessing"/"bodyForces.dat",
                IOstreamOption(),
                restart ? IOstreamOption::APPEND_APP : IOstreamOption::NO_APPEND
            )
        );
        forceFile().precision(10);
        if (!restart)
        {
            forceFile() << "# t";
            forAll(bodies, b)
            {
                forceFile() << ' ' << bodies[b].name() << "_Fx "
                    << bodies[b].name() << "_Fy " << bodies[b].name() << "_Fz";
            }
            forceFile() << "  (force of the material on the body, N)" << nl;
        }
    }

    // Stable time step at the start, printed once. Recomputed every cycle
    // below because density and stress will change once the physics is in.
    const scalar deltaT0 = solid.stableDeltaT();
    const label oldPrecision = Info().precision(16);

    Info<< "Smallest cell size  minDelta = "
        << gMin(cbrt(mesh.V().field())) << nl
        << "P-wave speed        c = "
        << Foam::sqrt(solid.material().waveModulus()/solid.rho0())
        << nl
        << "Stable time step    deltaT = " << deltaT0 << nl << endl;

    Info().precision(oldPrecision);

    printDiagnostics(solid);

    Info<< nl << "Starting time loop" << nl << endl;

    // Time step (decision 0081): the stable step; with adjustableRunTime
    // shortened so that steps land on the write times and on endTime. The write
    // times are counted from time zero, not from the start time, so a run
    // restarted from a write time takes the same steps as the run that wrote
    // it (OpenFOAM's own adjustment counts from the start time, which
    // changes the step in the last digits and lets restarts part). Never
    // longer than the stable step (OpenFOAM rounds the number of steps to
    // the nearest, which can lengthen it by up to half).
    const bool adjustWrite =
        runTime.controlDict().get<word>("writeControl") == "adjustableRunTime";
    const scalar writeInterval =
        runTime.controlDict().get<scalar>("writeInterval");
    const scalar endTime = runTime.endTime().value();

    while (runTime.run())
    {
        const scalar t = runTime.value();
        scalar dtStep = solid.stableDeltaT();
        scalar span = -1;
        if (adjustWrite && writeInterval > 0)
        {
            const scalar k =
                std::floor(t/writeInterval*(1 + 1e-12) + 1e-9) + 1;
            span = min(endTime, k*writeInterval) - t;
        }
        if (span > 0)
        {
            const scalar n = std::ceil(span/dtStep*(1 - 1e-12));
            dtStep = span/max(n, scalar(1));
        }
        runTime.setDeltaT(dtStep, false);
        ++runTime;

        Info<< "Time = " << runTime.timeName()
            << "  deltaT = " << runTime.deltaTValue() << nl;

        const scalar dt = runTime.deltaTValue();

        solid.lagrangeStep(dt);     // WP 2.3

        moveBodies(solid, bodies);  // WP 4.1, 4.3d: bodies to t

        solid.remap(dt);            // WP 2.4

        if (cut)
        {
            // WP 4.3b: cut cells at the new pose (moving bodies only)
            moveImmersedWalls(solid, cut.ref(), bodies);

            if (forceFile)
            {
                forceFile() << runTime.value();
            }
            forAll(bodies, b)
            {
                const vector f(-solid.bodyForce(b));
                if (forceFile)
                {
                    forceFile() << ' ' << f.x() << ' ' << f.y() << ' ' << f.z();
                }
            }
            if (forceFile)
            {
                forceFile() << nl;
            }
        }

        if (gas)
        {
            gas->correctFractions();    // WP F.1: phiS, phiG, k
            gas->solve(dt);             // WP F.2, F.3: pore pressure
        }

        if (runTime.write())
        {
            printDiagnostics(solid);
            solid.writeState();     // WP 6.2: restart state (decision 0081)
            bodies.write();
            forAll(solid.bodyNames(), b)
            {
                Info<< "Force on body " << bodies[b].name() << " = "
                    << -solid.bodyForce(b) << " N" << nl;
            }
            if (solid.movingWalls())
            {
                Info<< "Moving bodies: " << solid.nHandOver()
                    << " hand-overs from covered cells, mass lost "
                    << solid.lostMass() << " kg" << nl;
            }
            if (cut)
            {
                cut->write();       // alphaBody, bodyFaceOpen, ... (WP 4.2)
            }
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
