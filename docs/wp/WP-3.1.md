# WP 3.1: C++ gravity and wall boundaries (smooth, frictional)

State: done (2026-11-15)
Branch: `claude/project-thread-ihzn3r`

## Goal

celFoam's cycle (eulerianSolid) settles a full tank under gravity with smooth and
Coulomb-frictional walls and matches the prototype of WP 3.1p. First step of the C++
silo track that settles rules 12 and 13 (decision 0029) with one code.

## Read first

- `docs/decisions/0020`, `0029`, `0018`; `prototype/silo2d.py` (checks in `verify()`)
- `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C`

## Build

- Gravity (`constant/g`, `gravityRamp`), damping (`damping`) and walls
  (`materialProperties` `walls`) in `eulerianSolid`; hydrostatic reconstruction and wall
  Riemann solution in the Lagrange step (decision 0030).
- Test driver `tests/silo/siloCheck` (`celSiloCheck`), checks `tests/silo/silo_checks.py`,
  `tests/test_silo_settle.sh`.

## Done when

`tests/test_silo_settle.sh` passes: the WP 3.1p checks on the C++ cycle (lake at rest
below 1e-10 m/s; smooth column bottom within 2 percent, K0 within 5 percent; frictional
tank Janssen within 10 percent, K in 0.27 to 0.45, wall ratio within 5 percent where it
slips, force balance 1 percent; settled and steady), and 2 processors match serial.

## Out of scope

- Void and free surface (WP 3.2), the outlet (3.3), the discharge (3.4).
- Face shear cap (decision 0023): needed before the collapse runs in C++ (WP 3.2).

## Log

- 2026-11-15: Done. All checks pass in 57 s (numbers in decision 0030). Lake at rest was
  exact first time. Found a one-cell shear layer at frictional walls: the Hancock
  predictor ignored the wall reaction; fixed by using the wall Riemann solution in the
  predictor. Two test-script slips cost time: sorting cells by rounded centres (use
  floor of x / cell) and a master-only reduction in the driver (hung in parallel).
  Next: WP 3.2 (C++ void and free surface, with the shear cap) and the collapse with
  neither rule 12 nor 13; then 3.3, 3.4 and the E sweep of the discharge (decision 0029).
