/*---------------------------------------------------------------------------*\
Application
    celFoam

Description
    Explicit Coupled Eulerian-Lagrangian solver for large-deformation solids.

    Thin driver (WP 2.1): the Eulerian state and the steps of the explicit
    cycle live in library celEulerian (class eulerianSolid). Each cycle:
    stable time step, Lagrange step, remap, write. The Lagrange step
    (WP 2.3) and remap (WP 2.4) are still empty; contact comes in WP 4.x
    (see docs/work-packages.md).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"

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

        // WP 4.x: Eulerian-Lagrangian contact forces

        solid.remap(dt);            // WP 2.4

        if (runTime.write())
        {
            printDiagnostics(solid);
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
