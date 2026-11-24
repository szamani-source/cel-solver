# WP F.1: Gas fractions and the celFluid skeleton

State: done (decision 0126)
Branch: claude/project-thread-0rpfg3

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

- 2026-11-24: done. `src/celFluid` (poreGas), celFoam calls it after the remap;
  `tests/test_gas_fractions.sh` passes (about 80 s): packing 0.048 %, k 7e-16 against the
  script, solid fields byte for byte with no file / `model none` / `model poreGas`,
  2 processors 1e-13. In code phi_s = rho/rhoP because rho is already alpha times the
  material density. The 2-processor check runs plain Mohr-Coulomb: with the T0 rule set
  the solid itself parts by 1e-5 (known issue, decision 0080). `drag` keyword read,
  only ergunWenYu accepted until F.2. Next: F.2 (pore pressure, dragModel selector).

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
