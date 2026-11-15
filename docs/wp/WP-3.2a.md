# WP 3.2a: C++ void and free surface (split from WP 3.2)

State: done (2026-11-15)
Branch: `claude/project-thread-85hl11`

## Goal

One material plus void in celFoam's cycle (eulerianSolid), checked on the prototype's
free-surface cases (WP 3.2p, decision 0021). The collapse with the shear cap is WP 3.2b.

## Read first

- `docs/decisions/0021`, `0029`, `0030`; `prototype/void2d.py` (`void_sweep`, `verify()`)
- `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C`, `eulerianSolidRemap.C`

## Build

- `freeSurface` switch and `alpha` field in `eulerianSolid`; slab sides, homogenised
  cells, free faces and void closure in the Lagrange step; slab remap in
  `remapFreeSurface` (decision 0031).
- Test driver `tests/void/voidCheck` (`celVoidCheck`), checks `tests/void/void_checks.py`,
  `tests/test_void.sh`.

## Done when

`tests/test_void.sh` passes: lake at rest with a half-full top cell (speeds <= 1e-8,
no stress in void), block translating 50 cells (mass 1e-10, stress 1e-6 Pa, edges within
3 cells, travel 0.1 %), free fall (g within 0.5 %) and landing, cohesionless block
pulled apart (void within 5 % of 2 v t, no tension), and 2 processors match serial.

## Out of scope

- Face shear cap (rule 8), top layers spreading sideways (rule 9), column collapse at
  E = 0.1, 1, 10 MPa: WP 3.2b.
- Outlet (rule 10): WP 3.3.

## Log

- 2026-11-15: Done in one session; every check passed at the first run apart from two
  test-script slips (floor force sign, parallel stress scale). Numbers in decision 0031.
