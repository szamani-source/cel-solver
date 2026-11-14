# Status

Last updated: 2026-11-14

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 0, project setup.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | In progress | celFoam builds against v2512 and test_smokeBox passes locally (2026-11-11); CI has not run because no GitHub repository exists yet |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Draft | `docs/spec/target-problems.md` waits for Sia's answers |

## Waiting on Sia

- GitHub account and an empty repository (needed to push and run CI).
- Answers to the open questions in `docs/spec/target-problems.md`.
- Confirm GPL-3.0 (decision 0003).

## Next

1. Sia creates a GitHub repository; push this code and confirm CI passes (closes 0.1).
2. Sia answers the open questions in `docs/spec/target-problems.md` (closes 0.3).
3. Start WP 1.1 (brief: `docs/wp/WP-1.1.md`) and WP 1.4 (brief: `docs/wp/WP-1.4.md`). Both are Python and need neither GitHub nor OpenFOAM.
4. WP 0.4: fast OpenFOAM setup for fresh cloud sessions (a from-source build takes about an hour).

## Known issues

- None yet.
