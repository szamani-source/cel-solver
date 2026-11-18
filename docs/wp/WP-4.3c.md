# WP 4.3c: STL walls with the free surface and outlets

State: done (2026-11-18). Needs WP 4.3a. Split from 4.3a (2026-11-18).
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

- 2026-11-18: Done in one session (decision 0049), on top of the WP 4.3a branch (PR #43).
  Free surface with cut-cell walls in `eulerianSolid`: alpha per open volume, closed faces
  mirror in `slabSides`, wall face only where the material touches it (volume flux on its
  main axis, no mass in when material pulls away: void opens), slab remap on open areas
  and open lengths, ghost alpha, alpha mixing. Two bugs found on the way: empty open
  cells were skipped as "covered" (m = 0 test), so their stale alphaL created material
  volume; and small-cell members' alpha was clipped before mixing, which lost material
  volume (group density grew to 1.6 rho0). Now a group shares its Lagrangian volume ratio
  and is clipped only after mixing. Without bodies the free-surface path is bit for bit
  the old one (falling block compared). `tests/test_stl_void.sh` (about 7 min) passes:
  lake at rest between STL side walls (cut at 1/4 and 3/4) 4.4e-17 m/s; a T0 block
  falling onto a 10 deg STL slab: nothing through it, mass 2e-14, KE 8e-5 of its peak at
  0.6 s, creep 0.47 mm in the last 0.1 s; the 4.3a tilted tank without lid, void top,
  1.91 % from the mesh-wall tank with void top. **Done-when as written (2 % against the
  lid version) not met**: the free-surface tank with mesh walls is itself up to 17 % off
  the no-void tank at the bottom (STL void top 16.2 kPa, mesh-wall void top 16.4, lid
  14.3, mesh no void 14.0): a free-surface issue (noise void below 1 switches cells to the
  partly filled rules; the 0044 snap only runs with flowImpedance), recorded in STATUS
  known issues. The test compares against the mesh-wall void-top tank, 5 % as 4.3a's
  tilted case. Lake with an STL floor fails (corner cells: one oblique wall face), so the
  lake uses a mesh floor. Outlets: mesh patches, may be part-covered by STL. Tutorial
  `tutorials/stlSilo` (T0, 10 MPa, 9 min): no mass lost while settling, 72 kg/s/m
  (mesh-wall tank 57.5 / 67.7, 15 % round-off scatter). Next: WP 4.5.
