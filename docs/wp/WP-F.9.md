# WP F.9: Restart, MPI, tutorial, manual, v2.0.0

State: not started
Branch: assigned by the thread

## Goal

Release version 2.0.0 with the gas model.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- Restart of p and U_g, MPI checks, `tutorials/gasHopper`, manual chapter, VERSION, release notes, CHANGELOG

## Done when

- `tests/test_gas_parallel.sh`: restart and 2/4 processors match serial to 1e-12 on the F.5 and F.8 cases; `tests/test_manual.sh` passes with the gas chapter.

## Out of scope

- Temperature, liquids (out of v2).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
