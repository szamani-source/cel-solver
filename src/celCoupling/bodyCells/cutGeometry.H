/*---------------------------------------------------------------------------*\
Namespace
    Foam::cutGeometry

Description
    Exact geometry of a closed triangulated surface cut by a convex cell
    (WP 4.2, decision 0042), built bottom up so that each level is exact:

      edge   covered length: crossings of the segment with the triangles,
             entering where the triangle's outward normal faces against the
             segment, leaving where it faces along it
      face   covered area:   Green's theorem over the boundary of the covered
             part, which is the triangles' cut lines clipped to the face
             (oriented with the body on their left) plus the covered parts
             of the face's edges
      cell   covered volume: the divergence theorem over the boundary of the
             covered part, which is the triangles clipped to the cell plus the
             covered part of each face (x - c).n is constant on a planar face)

    The surface must be closed with outward normals (STL convention). An edge
    or face without crossings is wholly inside or outside; the caller decides
    which with an inside test. Degenerate contact (a surface exactly through
    a mesh edge or vertex) is treated by merging crossings closer than a
    relative tolerance; results there are exact only to that tolerance.

SourceFiles
    cutGeometry.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_cutGeometry_H
#define Foam_cutGeometry_H

#include "pointField.H"
#include "labelList.H"
#include "triSurface.H"
#include "DynamicList.H"

namespace Foam
{
namespace cutGeometry
{

//- Covered fraction of segment a-b. Returns -1 when the segment crosses no
//  triangle (then wholly inside or outside: the caller tests).
scalar edgeCovered
(
    const point& a,
    const point& b,
    const triSurface& s,
    const labelUList& tris
);

//- Covered area of a planar convex polygon (points in right-hand order
//  about its unit normal e). edgeFrac[i]: covered fraction of edge i
//  (pts[i] to pts[i+1]), each 0 to 1. Returns -1 when no triangle cuts the
//  polygon and no edge is partly covered.
scalar faceCovered
(
    const UList<point>& pts,
    const vector& e,
    const UList<scalar>& edgeFrac,
    const triSurface& s,
    const labelUList& tris
);

//- Covered volume of a convex cell given by its faces: planes (unit outward
//  normal n[f] through point p[f]), covered area of each face cA[f], cell
//  centre c. Triangles are clipped to the cell.
scalar cellCovered
(
    const UList<vector>& n,
    const UList<point>& p,
    const UList<scalar>& cA,
    const point& c,
    const triSurface& s,
    const labelUList& tris
);

} // End namespace cutGeometry
} // End namespace Foam

#endif
