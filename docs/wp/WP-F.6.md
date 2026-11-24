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

## Out of scope

- Drag on grains (F.7).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
