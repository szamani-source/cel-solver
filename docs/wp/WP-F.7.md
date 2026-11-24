# WP F.7: Gas through and around grains; more drag laws

State: done (decision 0132)
Branch: assigned by the thread

## Goal

Gas momentum in pores and void with implicit drag on both phases, and the user drag laws beside the default `ergunWenYu`: `ergun`, `wenYu`, `diFelice`, `syamlalOBrien`, `beetstra` (Sia: several laws for users, `ergunWenYu` for testing).

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- New `dragModel` types, implicit in both momentum updates; reference values added to `tools/gas_reference.py`

## Done when

- `tests/test_drag_laws.sh`: each law matches `tools/gas_reference.py` to 1e-10 over gas fraction 0.35 to 1 and Re 1e-3 to 1e3. `tests/test_fluidisation.sh` (with `ergunWenYu`): U_mf against Wen-Yu within 10 %; dilute settling against Richardson-Zaki within 10 %; F.2 to F.5 results within 1 % of Stage A.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Validation pack (Sia, 2026-11-24): this package delivers pack test 2, onset of fluidisation: raise the airflow until the bed lifts, compare with the published minimum fluidisation velocity (Wen and Yu), as `tests/validation/gas/<name>/` with an `Allrun` that runs it on 1, 2 and 4 processors, checks that they match and writes the run times to `times.csv`.

## Out of scope

- Collapse validation (F.8).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: done (decision 0132). Drag laws ergun, wenYu, diFelice, syamlalOBrien,
  beetstra (`tests/test_drag_laws.sh`, celDragTable); pore gas momentum in gasFlow
  (poreInertia, inertialGasDrag in poreGasFlow.C); inflow as Function1;
  `tests/test_fluidisation.sh` (relaxation, onset -2.2 % of Wen-Yu, settling +7.6 % of
  Richardson-Zaki, Stage A 0.1 %, 2/4 processors 7e-11 up to the onset); pack test 2
  `tests/validation/gas/fluidOnset`. Learned: the generic slip bracket assumed K grows
  with u (Syamlal-O'Brien does not: 18 % error); boundary faces need the inertia too
  (5 % fast relaxation without); the explicit reaction of the pore gas inertia destabilises
  a lifting bed (left out, `inertiaForce`); a lifted bed parts serial from parallel
  (round-off at the Mohr-Coulomb apex). Next: F.8.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
