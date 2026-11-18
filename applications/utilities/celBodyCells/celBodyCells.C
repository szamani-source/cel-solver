/*---------------------------------------------------------------------------*\
Application
    celBodyCells

Description
    Write which cells and faces the rigid bodies of constant/rigidBodies
    cover, at each selected time (bodies moved to that time), without running
    any physics (WP 4.2): alphaBody, bodyFaceOpen, wallArea, wallVelocity,
    mixNeighbour. Prints the covered volume per time and the search cost.

Usage
    celBodyCells [-time ...] [-parallel]
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "rigidBodies.H"
#include "bodyCells.H"
#include "clockTime.H"

int main(int argc, char *argv[])
{
    argList::addNote("Cells and faces covered by rigid bodies (WP 4.2)");
    timeSelector::addOptions(true, false);

    #include "setRootCase.H"
    #include "createTime.H"
    instantList times = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    rigidBodies bodies(runTime);
    bodyCells bc(mesh, bodies);

    const label oldPrecision = Info().precision(16);
    for (const instant& t : times)
    {
        runTime.setTime(t, 0);
        bodies.move();

        clockTime clock;
        bc.update();
        const scalar cost = clock.timeIncrement();

        const scalar covered =
            gSum(bc.alphaBody().primitiveField()*mesh.V().field());
        Info<< "Time = " << runTime.timeName()
            << "  covered volume = " << covered
            << "  cells searched = " << returnReduce(bc.nSearched(), sumOp<label>())
            << "  cut cells = " << returnReduce(bc.nCut(), sumOp<label>())
            << "  search time = " << returnReduce(cost, maxOp<scalar>()) << " s"
            << endl;
        bc.write();
    }
    Info().precision(oldPrecision);

    Info<< "End" << nl << endl;
    return 0;
}
