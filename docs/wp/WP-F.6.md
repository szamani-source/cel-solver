# WP F.6: Gas momentum in the void

State: not started
Branch: assigned by the thread

## Goal

A gas velocity field solved in the void (no grains yet): low-Mach, implicit pressure, time step set by the solid.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- poreGas -> gasFlow model, U_g field, PISO-like loop on fvMatrix

## Done when

- `tests/test_gas_flow.sh`: plane channel against Poiseuille within 1 %; still air in a closed box stays below 1e-10 m/s; time step unchanged from the solid one.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

## Out of scope

- Drag on grains (F.7).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
