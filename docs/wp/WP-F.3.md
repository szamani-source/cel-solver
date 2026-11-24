# WP F.3: Effective stress and volume-change source

State: done (decision 0128)
Branch: assigned by the thread

## Goal

Two-way coupling: constitutive models see sigma' = sigma + p I, momentum sees sigma' - p I; solid volume change and velocity enter the gas balance; implicit enough to keep the v1 step with fine grains.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- Hooks in eulerianSolid (stress update, momentum), poreGas source from u_s
- Possibly split into F.3a (one-way) and F.3b (two-way) if it grows

## Done when

- `tests/test_consolidation.sh`: gas-filled column loaded at the top, degree of consolidation against Terzaghi with c = k/(mu (e/p0 + 1/M)) within 3 %; column at rest stays at p0 to 1e-6 relative; 0.1 mm grains run at the v1 step without blow-up.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

Validation pack (Sia, 2026-11-24): this package delivers pack test 3, pore-pressure decay in a suddenly loaded fine-powder column against the exact diffusion (Terzaghi) solution, as `tests/validation/gas/<name>/` with an `Allrun` that runs it on 1, 2 and 4 processors, checks that they match and writes the run times to `times.csv`.

## Out of scope

- Void and free surface (F.4).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: done (decision 0128). Effective stress (sigma is sigma'; cell force
  -sum(p_f S_f)), gas carried with the grains matched to the remap's face masses,
  surface `loads`, `coupled yes|no`. `tests/test_consolidation.sh` passes (Terzaghi
  within 0.8 % from T = 0.05, at rest 2e-17, 2/4 processors 5e-13, mass 1e-15); pack
  test 3 `tests/validation/gas/terzaghiColumn`. Cost 1.4 to 1.5x on 2,000 cells, 2.5x on 50
  cells. Learned: the sudden load starts elastic waves the quasi-static theory leaves out
  (compare from T = 0.05 with damping, or make the column tall); F.2's inner relTol 0.1
  made 8 to 10 outer iterations per step, 1e-3 plus the error estimate makes 3. Next: F.4.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
