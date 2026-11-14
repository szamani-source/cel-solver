# Status

Last updated: 2026-11-14

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 0, project setup.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Done | Repository szamani-source/cel-solver; CI builds celFoam in opencfd/openfoam-dev:2512 and test_smokeBox passes (2026-11-14) |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Draft | `docs/spec/target-problems.md` waits for Sia's answers |

## Waiting on Sia

- Answers to the open questions in `docs/spec/target-problems.md`.
- Confirm GPL-3.0 (decision 0003).

## Next

1. Sia answers the open questions in `docs/spec/target-problems.md` (closes 0.3).
2. Start WP 1.1 (brief: `docs/wp/WP-1.1.md`) and WP 1.4 (brief: `docs/wp/WP-1.4.md`). Both are Python and need neither GitHub nor OpenFOAM.
3. WP 0.4: fast OpenFOAM setup for fresh cloud sessions (a from-source build takes about an hour).

## Known issues

- None yet.
