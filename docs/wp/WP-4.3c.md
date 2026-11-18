# WP 4.3c: STL walls with the free surface and outlets

State: not started. Needs WP 4.3a. Split from 4.3a (2026-11-18).
Branch: the thread's `claude/project-thread-*` branch

## Goal

The cut-cell walls of decision 0047 work with `freeSurface yes`: void, slab rules,
shear cap and rule 9 at cut cells, so the silo (WP 4.5) can have STL walls and a free top.

## Read first

- `docs/decisions/0031`, `0047`; WP 4.3a log; `src/celEulerian/eulerianSolid/eulerianSolidImmersed.C`

## Build

- Remove the FatalError in `setImmersedWalls` for `freeSurface`; alpha per open volume,
  slab sides (`lean_`) and `touches()` at cut faces, `remapFreeSurface` with open areas,
  small-cell mixing of alpha.

## Done when

`tests/test_stl_void.sh` passes: a block falling onto a 10 deg STL floor stops with no
material through it (1e-6); the WP 4.3a tank with a void top and no lid within 2 % of the
lid version; lake at rest in an STL tank keeps speeds below 1e-10 m/s.

## Out of scope

- Moving bodies (WP 4.3b); discharge (WP 4.5).

## Log

(none yet)
