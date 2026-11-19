# WP 6.3b: Validation report

State: done (2026-11-19). Split from WP 6.3 (scaling is 6.3a).
Branch: the thread's `claude/project-thread-*` branch

## Goal

One page listing every benchmark: reference, result, tolerance, status, the test that
reproduces it, and the known issues.

## Read first

- `docs/results/validation-report.md`; the decision notes it cites

## Done when

`docs/results/validation-report.md` is in the repository and shared with Sia.

## Out of scope

- New runs: the numbers come from the packages' tests and decision notes.
- Accuracy fixes (after the solver is complete, Sia 2026-11-18).

## Log

- 2026-11-19: Written from decisions 0019 to 0081 and the merged and open PRs (#55 wedge
  and fed silo, #57 projectile). Keep it current: each later benchmark adds a row.
