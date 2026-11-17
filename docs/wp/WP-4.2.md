# WP 4.2: Contact search near the body (cut cells)

State: not started. May start beside WP 4.0: every candidate scheme needs this geometry.
Branch: the thread's `claude/project-thread-*` branch

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

(none yet)
