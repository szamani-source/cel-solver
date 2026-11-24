# WP F.7: Gas through and around grains, Di Felice drag

State: not started
Branch: assigned by the thread

## Goal

Gas momentum in pores and void with implicit Di Felice drag on both phases.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- Drag model selectable (ergun, diFelice), implicit in both momentum updates

## Done when

- `tests/test_fluidisation.sh`: U_mf against Wen-Yu within 10 %; dilute settling against Richardson-Zaki within 10 %; F.2 to F.5 results within 1 % of Stage A.

## Out of scope

- Collapse validation (F.8).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
