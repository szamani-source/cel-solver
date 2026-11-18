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
    scalarList adhesion(bodies.size());
    forAll(bodies, b)
    {
        adhesion[b] = bodies[b].adhesion();
    }
    solid.setBodyAdhesion(adhesion);
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


void Foam::moveBodies(const eulerianSolid& solid, rigidBodies& bodies)
{
    const scalar t = solid.mesh().time().value();
    const label nCoupled = solid.bodyNames().size();
    List<vector> F(bodies.size(), Zero);
    List<vector> M(bodies.size(), Zero);
    scalarList dt(bodies.size(), 0);
    forAll(bodies, b)
    {
        const rigidBody& body = bodies[b];
        if (!body.free())
        {
            continue;
        }
        // Force and moment of the material on the body: minus those of the
        // body on the material, at the pose of the Lagrange step
        if (b < nCoupled)
        {
            F[b] = -solid.bodyForce(b);
            M[b] = -solid.bodyMoment(b, body.centre());
        }
        F[b] += body.mass()*solid.gravityAt(t);
        dt[b] = t - body.time();
    }
    bodies.move();
    forAll(bodies, b)
    {
        bodies[b].accelerate(F[b], M[b], dt[b]);
    }
}
