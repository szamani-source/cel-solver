# WP F.2: Pore pressure in a rigid bed

State: done (PR, decision 0127)
Branch: claude/project-thread-bc5lwa

## Goal

An implicit pore-pressure equation (decision 0125, item 4) with the solid frozen: gas mass balance with Darcy-Ergun flux, `fvMatrix` solve per step, pressure and flux boundaries.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- `poreGas::solve(dt)`, fields p, qGas; boundary types fixedPressure and sealed
- Run-time-selected `dragModel` (decision 0125, item 7) with `ergunWenYu` as the default; k = phi_g^2 / K, the inertial term as a Picard correction
- Tests run with `ergunWenYu`

## Done when

- `tests/test_pore_pressure.sh`: 1D permeameter pressure drop against Ergun within 1 % at two flow rates; step at a face against erfc within 2 % at three times; gas mass in a sealed box to 1e-12.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Validation pack (Sia, 2026-11-24): this package delivers pack test 1, packed-bed pressure drop against Ergun (air blown up through a fixed sand column), as `tests/validation/gas/<name>/` with an `Allrun` that runs it on 1, 2 and 4 processors, checks that they match and writes the run times to `times.csv`.

## Out of scope

- Any force on the solid (F.3).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: done (decision 0127). `poreGas::solve(dt)` (gauge p, Newton + Anderson,
  conservative finish), `gasDragModel` selector with `ergunWenYu`, inflow patches,
  `tests/test_pore_pressure.sh` passes (Ergun 1e-6, Wen-Yu branch, erfc 0.07 %, sealed
  non-orthogonal box mass to round-off, 2/4 processors 3e-11), validation pack test 1
  `tests/validation/gas/ergunBed` (0.008 %, 1/2/4 processors). Cost 1.09x (quiet) to 3.6x
  (sharp blob) of gas off. Learned: solve in gauge pressure; `fvMatrix::residual()` is
  wrong across processor faces for this use; Picard alone is far too slow with Ergun
  inertia. Next: F.3 (grain velocity, volume change, effective stress).

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
