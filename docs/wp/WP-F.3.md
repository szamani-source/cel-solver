# WP F.3: Effective stress and volume-change source

State: not started
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

## Out of scope

- Void and free surface (F.4).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
