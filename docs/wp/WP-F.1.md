# WP F.1: Gas fractions and the celFluid skeleton

State: not started
Branch: assigned by the thread

## Goal

A `src/celFluid` library (class `poreGas`) built by Allwmake and called by celFoam when `constant/gasProperties` exists. It computes phi_s = alpha rho / rho_p, phi_g and the Ergun permeability per cell and writes them.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- `src/celFluid/` (Make, poreGas.H/.C), `constant/gasProperties` (rhoP, d, T, p0, model none|poreGas)
- celFoam: construct and call `correctFractions()` after the remap

## Done when

- `tests/test_gas_fractions.sh`: settled column, phi_s equals the set packing within 0.5 %; k equals `tools/gas_reference.py` to 1e-12; every v1 test unchanged bit for bit with no gasProperties; 2 processors match serial.

## Out of scope

- Any pressure solve (F.2).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
