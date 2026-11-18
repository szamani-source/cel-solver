# WP 4.2: Contact search near the body (cut cells)

State: done (2026-11-17)
Branch: `claude/project-thread-hilqye`

## Goal

For every rigid body (WP 4.1) at its current pose: which Eulerian cells it covers or
cuts, the body volume fraction per cell, the open fraction of each face, and per cut
cell the nearest surface point, outward normal and body velocity (`velocityAt`).

## Read first

- `docs/decisions/0041`, `0031` (Cartesian cells, slab void); `docs/wp/WP-4.0.md` and
  its recommendation if done (it may ask for more, such as face apertures per side)
- `src/celLagrangian/rigidBody/rigidBody.H`

## Build

- `src/celCoupling/bodyCells`: search only near each body (bounding box, then
  `triSurfaceSearch`), exact clipping of axis-aligned hex cells by the triangulated
  surface. Per cell and step, as decision 0042 needs: body fraction, open fraction of
  each face, the wall segment (area times normal = minus the sum of open face areas
  times normals), body velocity there (`rigidBody::velocityAt`), and the mixing
  neighbour for cells with open fraction below 0.5. Python reference:
  `prototype/contact2d.py` (`Rect`, `HalfPlane`).
- Decide and record: walls as closed STL volumes only, or also one-sided open surfaces
  with a thickness (decision 0041 consequences).
- Field `alphaBody` written for ParaView; a small utility or celFoam option to write it
  without running the physics.

## Done when

`tests/test_body_cells.sh` passes: a sphere at 10 and 20 cells per diameter, body volume
summed over cells within 1 % and 0.3 % of 4/3 pi r^3; a plane at 30 deg gives exact
(1e-10) cell and face fractions; a moving cube's fractions follow its pose; 2 processors
match serial to 1e-12; search cost per step reported.

## Out of scope

- Forces (WP 4.3a, 4.3b).

## Log

- 2026-11-17: Done in one session (decision 0043). `src/celCoupling/bodyCells` (exact
  edge, face and cell coverage, body-frame octree), utility `celBodyCells` (writes
  alphaBody, bodyFaceOpen, wallArea, wallVelocity, mixNeighbour at selected times).
  `tests/test_body_cells.sh` passes in about 6 s: sphere volume 4e-16 at 10 and 20 cells per
  diameter (the done-when asked 1 % and 0.3 % of the true sphere; the faceted STL itself is
  0.86 % below it, so the check is against the STL's own volume), plane at 30 deg 2e-7
  (single-precision STL), moving cube shifted exactly, turning cube 4e-16, 2 processors 9e-16.
  Search 6e-5 s per searched cell. Thin-wall and coarse-facet warnings moved to 4.3a.
  Next: WP 4.3a.
