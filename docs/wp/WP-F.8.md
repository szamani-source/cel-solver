# WP F.8: Collapse in gas (end of Stage B)

State: in progress, split (F.8b, decision 0133)
Branch: claude/project-thread-ndfemq (stacked on F.7)

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

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

## Out of scope

- Release (F.9).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: read the Roche (2008) abstract (full text refused by HAL, AGU, hal.science); case `tests/gas/collapse_cases.py`, mesh shown to Sia, `tests/manual/test_fluidised_collapse.sh` (checkpointed). Found and fixed (option `surfaceCarry`) the free surface sucking gas through the grains; expanded bed made soft with the cap. Result: run-out 1.95 h0 (-66 %), front -70 %: FAIL; dry identical to main, gas mass 1e-14; parallel parts within 0.002 s and blew up; cost 2.3x. Learned: the defluidised base's floor friction holds the flow (no floor friction: 7.0 h0); void flow near a thin front is not robust. Next (F.8b): full paper, then floor/gas supply, robust void flow and drag guard, parallel from an apex start, 5 mm mesh.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
