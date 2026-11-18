# WP 4.3a: Fixed STL walls with friction

State: done (2026-11-18). Needs WP 4.0 (scheme, shown to Sia) and WP 4.2.
Branch: the thread's `claude/project-thread-*` branch

## Goal

Material in the Eulerian mesh rests on and slides along fixed STL walls with Coulomb
friction and no leakage, by the scheme WP 4.0 recommends. This is the wall interaction
Sia asked for (2026-11-17): silo walls and hopper inserts that do not follow the mesh.

## Read first

- `docs/decisions/0030` (mesh-boundary walls), `0041`, `0042`, `0043`; `prototype/contact2d.py`;
  `tests/test_silo_settle.sh` (the Janssen tank)

## Build

- From WP 4.2 (decision 0043): `bodyCells` gives the geometry; mixing across processor
  faces is not done there; add the thin-wall and coarse-facet warnings (decision 0042
  addendum).

- Wall friction per body in `constant/rigidBodies` (`friction <deg>`, as `walls` in
  materialProperties), the coupling at cut cells, wall force output per body.

## Done when

`tests/test_stl_walls.sh` passes (10 mm, coarse): the WP 3.1 frictional tank with its
side walls given as STL boxes aligned with the mesh gives Janssen within 2 % of the
mesh-boundary result; the same tank tilted 10 deg against the mesh (gravity tilted with
it) within 5 %; material inside walls below 1e-6 of the total; wall force balances the
weight to 1 %.

## Out of scope

- Moving bodies and forces on them (WP 4.3b); discharge (WP 4.5).

## Log

- 2026-11-18: Done in one session (decision 0047). Cut-cell walls in `eulerianSolid`
  (`setImmersedWalls`, `eulerianSolidImmersed.C`; open areas and volumes in the Lagrange
  step and remap), `friction <deg>` per body in `constant/rigidBodies`, glue
  `src/celCoupling/immersedWalls`, celFoam and celSiloCheck couple fixed bodies, per-body
  wall force and leak report, thin-wall and coarse-facet warnings in `bodyCells`.
  `tests/test_stl_walls.sh` passes (about 3 min): aligned STL walls (cells cut at 1/4 and
  3/4) 0.24 % from the mesh-wall tank; tank turned 10 deg with gravity 2.1 %; leak 1e-11;
  force balance 0.06 %. Learned: mixing small cells only after the remap left sliver
  cells (wall close to a mesh face) straining wildly and the tilted tank never settled
  (28 %); small cells now act as one cell with their partner in the Lagrange step.
  Not done: free surface with STL walls (split to WP 4.3c, needed by 4.5), mixing
  across processors, parallel run of the test. Tutorial `tutorials/stlTank`. Next:
  WP 4.3c.
