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
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Done, PR open | `prototype/constitutive/`; `tests/test_proto_constitutive.sh` passes in CI on PR #2 (2026-11-14) |

## Waiting on Sia

- Answers to the open questions in `docs/spec/target-problems.md`.
- Confirm GPL-3.0 (decision 0003).

## Next

1. Sia answers the open questions in `docs/spec/target-problems.md` (closes 0.3).
2. WP 1.1 (brief: `docs/wp/WP-1.1.md`) is in progress in its own thread. WP 1.4 is done once its PR merges; the C++ ports (2.2, 2.5a, 2.5b) use its results as reference.
3. WP 0.4: fast OpenFOAM setup for fresh cloud sessions (a from-source build takes about an hour).

## Known issues

- Mohr-Coulomb prototype has no separate tension cut-off below the apex (decision 0004).
