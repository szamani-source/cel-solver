# WP F.4: Free surface, void and outlet for the pore gas

State: not started
Branch: assigned by the thread

## Goal

Void cells at the boundary pressure, partly filled cells with k from the filled part, the free surface and the outlet as pressure faces.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- poreGas treatment by alpha; rules for alpha below the 1e-4 void threshold as in v1

## Done when

- `tests/test_gas_venting.sh`: over-pressured column with a free top vents within 3 % of the closed form; v1 tests test_void, test_outlet unchanged; outlet face holds p0.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

## Out of scope

- Gas flow inside the void (F.6).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
