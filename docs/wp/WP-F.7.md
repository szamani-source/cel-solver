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

## Out of scope

- Collapse validation (F.8).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
