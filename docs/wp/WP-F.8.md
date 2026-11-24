# WP F.8: Collapse in gas (end of Stage B)

State: not started
Branch: assigned by the thread

## Goal

Initially fluidised dam break against Roche et al. (2010); read the paper first (the HAL page refused automated reading in F.0). Fall-back: Rondon et al. (2011).

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- Case, mesh shown to Sia, checkpointed run

## Done when

- `tests/manual/test_fluidised_collapse.sh`: run-out and front speed within 20 % of the paper; dry collapse with gas off unchanged.

## Out of scope

- Release (F.9).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
