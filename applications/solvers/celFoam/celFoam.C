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
    bodies are coupled in WP 4.3b (see docs/work-packages.md).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "rigidBodies.H"
#include "bodyCells.H"
#include "immersedWalls.H"

// Print total mass and momentum with enough digits to compare runs
static void printDiagnostics(const Foam::eulerianSolid& solid)
{
    using namespace Foam;

    const label oldPrecision = Info().precision(16);

    Info<< "Total mass      = " << solid.totalMass() << nl
        << "Total momentum  = " << solid.totalMomentum() << endl;

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

    // Fixed bodies act on the material as cut-cell walls (WP 4.3a,
    // decision 0047). Moving bodies are not coupled yet (WP 4.3b).
    autoPtr<bodyCells> cut;
    bool anyMoving = false;
    forAll(bodies, b)
    {
        anyMoving = anyMoving || bodies[b].moving();
    }
    if (bodies.size() && !anyMoving)
    {
        cut.reset(new bodyCells(mesh, bodies));
        setImmersedWalls(solid, cut.ref(), bodies);
    }
    else if (anyMoving)
    {
        Info<< "Moving rigid bodies: not coupled to the material yet "
            << "(WP 4.3b)" << nl << endl;
    }

    // Stable time step at the start, printed once. Recomputed every cycle
    // below because density and stress will change once the physics is in.
    const scalar deltaT0 = solid.stableDeltaT();
    const label oldPrecision = Info().precision(16);

    Info<< "Smallest cell size  minDelta = "
        << gMin(cbrt(mesh.V().field())) << nl
        << "P-wave speed        c = "
        << Foam::sqrt(solid.material().waveModulus()/gMin(solid.rho().primitiveField()))
        << nl
        << "Stable time step    deltaT = " << deltaT0 << nl << endl;

    Info().precision(oldPrecision);

    printDiagnostics(solid);

    Info<< nl << "Starting time loop" << nl << endl;

    while (runTime.run())
    {
        runTime.setDeltaT(solid.stableDeltaT());
        ++runTime;

        Info<< "Time = " << runTime.timeName()
            << "  deltaT = " << runTime.deltaTValue() << nl;

        const scalar dt = runTime.deltaTValue();

        solid.lagrangeStep(dt);     // WP 2.3

        bodies.move();              // WP 4.1: prescribed motion to t


        solid.remap(dt);            // WP 2.4

        if (runTime.write())
        {
            printDiagnostics(solid);
            bodies.write();
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
