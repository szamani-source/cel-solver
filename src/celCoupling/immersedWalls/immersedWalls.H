/*---------------------------------------------------------------------------*\
Function
    Foam::setImmersedWalls

Description
    Hand the rigid bodies' cut-cell geometry (bodyCells, WP 4.2) and their
    wall friction (constant/rigidBodies, friction <deg>) to the Eulerian
    solid, which then treats them as cut-cell walls (WP 4.3a, decision 0047).
    Updates the geometry at the bodies' current pose first.

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

} // End namespace Foam

#endif
