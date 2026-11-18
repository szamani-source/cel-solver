/*---------------------------------------------------------------------------*\
    setImmersedWalls (WP 4.3a): see immersedWalls.H
\*---------------------------------------------------------------------------*/

#include "immersedWalls.H"

void Foam::setImmersedWalls
(
    eulerianSolid& solid,
    bodyCells& cells,
    const rigidBodies& bodies
)
{
    cells.update();
    wordList names(bodies.size());
    scalarList friction(bodies.size());
    forAll(bodies, b)
    {
        names[b] = bodies[b].name();
        friction[b] = bodies[b].friction();
    }
    solid.setImmersedWalls
    (
        cells.alphaBody(),
        cells.faceOpen(),
        cells.wallVelocity(),
        cells.cellBody(),
        names,
        friction
    );
    solid.setMovingWalls(anyMoving(bodies));
}


bool Foam::anyMoving(const rigidBodies& bodies)
{
    forAll(bodies, b)
    {
        if (bodies[b].moving())
        {
            return true;
        }
    }
    return false;
}


void Foam::moveImmersedWalls
(
    eulerianSolid& solid,
    bodyCells& cells,
    const rigidBodies& bodies
)
{
    if (!solid.movingWalls())
    {
        return;
    }
    cells.update();
    solid.moveImmersedWalls
    (
        cells.alphaBody(),
        cells.faceOpen(),
        cells.wallVelocity(),
        cells.cellBody()
    );
}
