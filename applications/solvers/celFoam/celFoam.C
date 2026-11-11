/*---------------------------------------------------------------------------*\
Application
    celFoam

Description
    Explicit Coupled Eulerian-Lagrangian solver for large-deformation solids.

    Work package 0.1 skeleton: reads one material, sets the stable explicit
    time step from the elastic P-wave speed and the smallest cell, and steps
    time. The Lagrange step, contact and remap are added in later packages
    (see docs/work-packages.md).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Explicit Coupled Eulerian-Lagrangian solver (skeleton, WP 0.1)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // Stable time step: dt = maxCo * minDelta / c, with c the P-wave speed
    const scalarField delta(cbrt(mesh.V().field()));
    const scalar minDelta = gMin(delta);
    const scalar deltaT = maxCo*minDelta/waveSpeed.value();

    runTime.setDeltaT(deltaT);

    Info<< "Smallest cell size  minDelta = " << minDelta << nl
        << "P-wave speed        c = " << waveSpeed.value() << nl
        << "Stable time step    deltaT = " << deltaT << nl << endl;

    Info<< "Starting time loop" << nl << endl;

    while (runTime.run())
    {
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl;

        // Placeholders, filled by later work packages:
        //   2.x  stress update and Lagrange step
        //   4.x  Eulerian-Lagrangian contact
        //   2.4  remap

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
