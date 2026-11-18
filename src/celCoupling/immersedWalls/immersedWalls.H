/*---------------------------------------------------------------------------*\
Function
    Foam::setImmersedWalls

Description
    Hand the rigid bodies' cut-cell geometry (bodyCells, WP 4.2) and their
    wall friction (constant/rigidBodies, friction <deg>) to the Eulerian
    solid, which then treats them as cut-cell walls (WP 4.3a, decision 0047).
    Updates the geometry at the bodies' current pose first. If any body
    moves, the solid is told so (setMovingWalls; needs freeSurface yes).

    moveImmersedWalls (WP 4.3b, decision 0053): after the remap of each
    step, with the bodies moved to the new time, recompute the cut cells and
    carry the material onto them. Does nothing if no body moves.

    moveBodies (WP 4.3d, decision 0058): move the bodies to the current
    time, then change each free body's velocity by the force and moment of
    the material in the Lagrange step just taken (at the old pose) plus
    gravity. Prescribed bodies just move. Deformable bodies (WP 5.3,
    decision 0072) get the force of each wall cell on the nodes of the
    nearest surface triangle and step their elements to the current time.

    The cycle with moving bodies:
        solid.lagrangeStep(dt); moveBodies(solid, bodies); solid.remap(dt);
        moveImmersedWalls(solid, cells, bodies);

SourceFiles
    immersedWalls.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_immersedWalls_H
#define Foam_immersedWalls_H

#include "eulerianSolid.H"
#include "bodyCells.H"

namespace Foam
{

void setImmersedWalls
(
    eulerianSolid& solid,
    bodyCells& cells,
    const rigidBodies& bodies
);

//- True if any body has a prescribed motion
bool anyMoving(const rigidBodies& bodies);

void moveImmersedWalls
(
    eulerianSolid& solid,
    bodyCells& cells,
    const rigidBodies& bodies
);

void moveBodies(const eulerianSolid& solid, rigidBodies& bodies);

} // End namespace Foam

#endif
