# WP 5.4: Projectile benchmark

State: done (2026-11-18). Needs WP 4.3d, 5.3b.
Branch: the thread's `claude/project-thread-*` branch

## Goal

Target problem T3: a deformable projectile into soil, free and coupled both ways.

## Read first

- `docs/decisions/0058`, `0072`, `0073`, `0074`; `tests/projectile/projectile_checks.py`

## Done when

`tests/test_projectile.sh`: a steel slug fired into Tresca clay stops within 5 % of the
same slug as a rigid free body (the reference run), both within 30 % of the Prandtl
energy estimate, at rest, mass kept, density within 3 %.

## Log

- 2026-11-18: Done (decision 0074). Steel 23.75 mm against rigid 23.84 mm (-0.3 %).
  Both 26 % short of the energy estimate (32.3 mm), worse at lower speed: known issue.
  A free deformable body now starts with its `velocity`.
