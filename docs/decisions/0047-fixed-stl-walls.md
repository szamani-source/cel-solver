# 0047 Fixed STL walls as cut-cell walls in the C++ solver (WP 4.3a)

Date: 2026-11-18. Status: accepted.

## Context

Decision 0042 chose cut-cell walls for STL bodies; decision 0043 (WP 4.2) gives the exact
geometry per cell and face. WP 4.3a puts the walls into the C++ Lagrange step and remap,
for fixed bodies with Coulomb friction, and checks them on the WP 3.1 Janssen tank.

## Decision

1. **Open volumes and areas.** Material fills the open volume `V (1 - alphaBody)` of each
   cell and crosses only the open part of each face: face forces, `phiL`, the velocity
   gradient and the Hancock predictor use the open area. Rho stays the material density
   (mass over open volume). A cell with an open fraction below 1e-6 is covered: it
   holds no material, and its faces count as closed for both neighbours.
2. **Wall face.** A wall cell has one more face, area vector minus the sum of the open
   face areas times their outward normals (recomputed in `eulerianSolid` from the
   effective open fractions, so a uniform stress gives exactly no force, and including
   cells whose face lies on a body but which are not cut themselves). It is solved with
   the mesh-wall Riemann solution of decision 0030 (`wallFace`: no tension, Coulomb cap
   with the body's `friction`), in the body's frame, from the cell state plus its
   Hancock increment. The body's force is minus the sum of these face forces.
3. **Ghosts.** Covered cells next to material hold the mean of their open neighbours'
   state, stress with the hydrostatic jump (`equilibriumJump`, now in its own header):
   the zero-gradient mirror ghost of a mesh wall. Reconstructions and remap limiters need
   no special case near a wall.
4. **Small cells** (open fraction below 0.5) are mixed with a neighbour (the largest open
   face to a cell that is not small; through up to 3 small cells in corners), on this
   processor. Mixing after the remap alone (decision 0042 item 2) was not enough: where
   the wall passes close to a mesh face, a sliver cell with a large open face took the
   neighbour's load on a tiny mass and volume, and its strain (an O(1) error even after
   mass-weighted mixing) set whole mesh columns sliding (the 10 deg tank never settled,
   walls carried half their load, 28 % error). Now the group (target plus its small
   cells) is one cell in the Lagrange step: summed force over summed mass, summed
   velocity gradient over summed volume, in the predictor and the update. Each keeps its
   own volume change for its density; after the remap the group's mass, momentum,
   stress (less the hydrostatic part) and history are averaged.
5. **Leak report.** Mass that enters through each body's wall faces (only where the
   material pulls away from it: the face is then a free surface, like an open top) is
   summed per body (`bodyMassIn`). No flux crosses a closed face, so covered cells hold
   no material by construction.
6. **Scope.** Without `freeSurface` only (STL walls with void and outlets: WP 4.3c).
   `celFoam` couples fixed bodies (no `velocity` or `omega`); moving bodies are not
   coupled yet (WP 4.3b). Mixing across processors is not done (small cells without a
   partner on their processor are counted and reported). Cut cells use the cell centre,
   not the open part's centroid, for reconstructions.
7. **Warnings** (decision 0042 addendum, in `bodyCells`): facets longer than 4 cells
   meeting at 10 to 45 deg (a curved surface faceted at the mesh scale); cut cells with no
   face neighbour more covered than themselves (a body thinner than about two cells).

## Consequences

- `tests/test_stl_walls.sh` (about 3 min): the WP 3.1 tank with STL side walls cutting
  cells at 1/4 and 3/4 matches the mesh-wall tank within 0.24 %; the whole tank (walls,
  bottom, smooth lid as free top) turned 10 deg against the mesh, gravity with it, within
  2.1 %; material through the walls 1e-11 of the total; weight balanced by the wall forces
  to 0.06 %.
- The tilted tank needs a lid: without void, the mesh outside the tank must be covered.
  With the free surface (WP 4.3c) the lid goes.
- Mesh-boundary wall code is unchanged; without bodies every path is the old one bit for
  bit (all other tests pass unchanged).
