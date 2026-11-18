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

    The cycle with moving bodies:
        solid.lagrangeStep(dt); bodies.move(); solid.remap(dt);
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

} // End namespace Foam

#endif
