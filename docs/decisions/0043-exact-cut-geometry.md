# 0043: Exact cut-cell geometry, built from edges up, searched in the body frame

Date: 2026-11-17 (WP 4.2)

## Context

Decision 0042 needs, per cell and step, the body fraction, the open fraction of each face,
and the wall face (area times normal) of every cut cell. The wall face is taken as minus
the sum of open face areas times their outward normals, so it is only as good as the face
fractions, and a wrong cell fraction moves mass that is not there.

## Decision

1. Exact geometry from the bottom up (`src/celCoupling/bodyCells/cutGeometry`):
   - an edge's covered length from its crossings with the STL triangles (entering where
     the outward normal faces against the edge, leaving where it faces along it; equal
     crossings of the same kind merged, for edges through a shared triangle edge);
   - a face's covered area by Green's theorem over its covered part, whose boundary is the
     triangles' cut lines clipped to the face plus the covered parts of its edges (only the
     edge fractions are needed, since all of an edge lies on one line);
   - a cell's covered volume by the divergence theorem: the triangles clipped to the cell
     plus the covered part of each face, where (x - c).n is constant.
   An edge, face or cell that no triangle reaches is decided by an inside test
   (`triSurfaceSearch::calcInside`) or by its faces. No sampling, no tolerances beyond
   round-off and the merging of equal crossings.
2. Each body's STL is searched in its own frame with one octree for the whole run; the mesh
   cells near the body are mapped into that frame, so moving bodies need no rebuild.
3. Only cells whose bounding box meets the body's are looked at. Bodies must be closed
   STLs with outward normals and must not overlap each other.
4. The mixing neighbour of a cut cell (decision 0042 item 2) is across the face with the
   largest open area facing away from the body; -1 when that face is a boundary or a
   processor face (mixing across processors is left to WP 4.3a).

## Consequences

- `tests/test_body_cells.sh`: a faceted sphere's summed covered volume equals the STL's own
  volume to 4e-16 at 10 and 20 cells per diameter; a cube moved by one cell gives the same
  fractions shifted by one cell to 3e-16; a turning cube keeps its volume to 4e-16; 2
  processors match serial to 9e-16. A 30 deg plane matches the exact fractions to 2e-7,
  which is the STL's single precision, not the method.
- Cost: about 6e-5 s per searched cell (sphere at 20 cells per diameter: 9261 cells, 0.55 s).
  Fixed walls need it once; a moving body every step would need a cheaper update (only
  the cells near the surface, or reuse between steps) before WP 4.3b.
- A surface passing exactly through a mesh edge or vertex relies on the merging of equal
  crossings (relative 1e-10); not tested.
- Not done (moved to WP 4.3a): warnings for STL features thinner than about two cells and
  for facets much larger than a cell on curved surfaces (decision 0042 addendum).
