# WP 4.3b: Moving bodies

State: done (2026-11-18). Needs WP 4.3a, 4.3c.
Branch: the thread's `claude/project-thread-*` branch

## Goal

A rigid STL body with prescribed motion (`velocity`, `omega` in `constant/rigidBodies`,
WP 4.1) pushes the material through its cut cells (decision 0042) and reports the force
the material puts on it (the reaction), with no material lost or passing through it.
Needed for WP 4.4 (footing and pile) and for inserts that move.

## Read first

- `docs/decisions/0042`, `0047`, `0051`; WP 4.3a and 4.3c logs;
  `src/celEulerian/eulerianSolid/eulerianSolidImmersed.C`, `remapFreeSurface`

## Build

- Each step: Lagrange step at the old pose (the wall face already moves with the body
  velocity), bodies moved, remap onto the open volume the body leaves (old open volume
  plus the wall's swept volume), then the cut cells recomputed at the new pose and the
  material carried over: mass, material volume and momentum kept per cell, a cell the
  body now covers hands its content to an open neighbour, a cell the body uncovers
  starts empty.
- Moving bodies need `freeSurface yes` (an uncovered cell must be able to be empty).
- Reaction output: per-body force in celFoam's log and `postProcessing/bodyForces.dat`.
- Drivers (`celFoam`, `celVoidCheck`) couple moving bodies through one celCoupling call.

## Done when

`tests/test_stl_push.sh` passes: a footing pushed into a T0 layer and pulled back out.
Mass kept to 1e-10 with nothing lost in hand-overs (no leakage); no void under the
footing while it pushes; material density within 2 % of the solid density; momentum
balance of the material (body, wall and gravity forces) within 1 %; after pulling out,
the force on the footing falls below 1 % of its peak (no tension).

## Out of scope

- Bodies whose motion follows from the forces (free bodies), deformable bodies (Phase 5).
- Moving bodies without the free surface; mixing across processors.
- Bearing capacity against theory (WP 4.4).

## Log

- 2026-11-18: Done in one session (decision 0053). `moveImmersedWalls` (eulerianSolid and
  celCoupling), remap target `Vt`, reaction file in celFoam, celVoidCheck with moving
  bodies and wall-force, lost-mass and density columns. `tests/test_stl_push.sh` passes
  (about 1 min) with the T0 rule set. Learned: dropping the small negative amounts that
  overdrawn group members hold when the body covers them created mass (3e-5); clipping
  group members before mixing squeezed the cells under the footing. Plain Mohr-Coulomb
  squeezes a nearly full corner cell (1.39 rho0); recorded, not fixed. celFoam checked by
  hand with a moving tank wall in `stlSilo` (2.3 s, 1462 hand-overs, nothing lost).
  Next: WP 4.4.
