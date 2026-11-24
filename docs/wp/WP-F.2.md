# WP F.2: Pore pressure in a rigid bed

State: not started
Branch: assigned by the thread

## Goal

An implicit pore-pressure equation (decision 0125, item 4) with the solid frozen: gas mass balance with Darcy-Ergun flux, `fvMatrix` solve per step, pressure and flux boundaries.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- `poreGas::solve(dt)`, fields p, qGas; boundary types fixedPressure and sealed
- Ergun inertial term as a Picard correction of k

## Done when

- `tests/test_pore_pressure.sh`: 1D permeameter pressure drop against Ergun within 1 % at two flow rates; step at a face against erfc within 2 % at three times; gas mass in a sealed box to 1e-12.

## Out of scope

- Any force on the solid (F.3).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
