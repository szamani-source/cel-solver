# WP F.7: Gas through and around grains; more drag laws

State: not started
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

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
