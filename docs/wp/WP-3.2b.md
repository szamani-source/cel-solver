# WP 3.2b: C++ column collapse: shear cap and rule 9 (split from WP 3.2)

State: done (2026-11-15)
Branch: `claude/project-thread-z7j1lb`

## Goal

Add the face shear cap (rule 8, decision 0023) and rule 9 (top layers spread sideways,
decision 0024) to the C++ cycle with void (WP 3.2a), then run the granular column
collapse at E = 0.1, 1 and 10 MPa with neither rule 12 nor 13 (decision 0029) and record
run-out and whether it stops, against the prototype and decision 0019.

## Read first

- `docs/decisions/0023`, `0024`, `0025`, `0029`, `0031`; `prototype/void2d.py`
  (`shear_cap`, `layer_across`, `column_collapse`)
- `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C`

## Build

- `materialProperties` switches `shearCap` and `spreadLayers` (default off);
  `capShear`, rule 9 in `slabSides`, slip velocities through `syncCoupled`.
- `celVoidCheck` gained `maxCo` (variable time step), `writeTimes`, and KE and capped
  face counts in `history.dat`.
- `tests/test_collapse.sh`, `tests/collapse/collapse_checks.py`.

## Done when

`tests/test_collapse.sh` passes: run-outs at 0.1, 1, 10 MPa recorded (a = 0.5) with mass
to 1e-10 and no stress in void; the 1 MPa runs (a = 0.5 and 2) within decision 0019's
bands, top kept, fitted slope at most 37 deg; 2 processors match serial at t = 0.1 s.

## Log

- 2026-11-15: Done in one session. Numbers and the parallel round-off finding in decision
  0032. Earlier C++ tests (void, settle, Lagrange step, smokeBox, skeleton) pass unchanged.
