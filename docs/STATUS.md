# Status

Last updated: 2026-11-14

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 0, project setup. Phase 1 prototypes started.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Done | Repository szamani-source/cel-solver; CI builds celFoam in opencfd/openfoam-dev:2512 and test_smokeBox passes (2026-11-14) |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Draft | `docs/spec/target-problems.md` waits for Sia's answers |
| 1.1 | 1D elastic Lagrange-plus-remap prototype | Done, in review (PR #3) | `prototype/lr1d.py`, `tests/test_proto_1d_elastic.sh` passes (order 1.88/1.83, monotone, conservative). Cell-centred instead of staggered: decision 0005 |
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Done | `prototype/constitutive/`; `tests/test_proto_constitutive.sh` passes in CI; PR #2 merged (2026-11-14) |

## Waiting on Sia

- Answers to the open questions in `docs/spec/target-problems.md`.
- Confirm GPL-3.0 (decision 0003).
- Review the WP 1.1 pull request, including decision 0005 (cell-centred Lagrange step instead of the staggered grid in the brief).

## Next

1. Sia answers the open questions in `docs/spec/target-problems.md` (closes 0.3).
2. After WP 1.1 merges, WP 1.2 (1D plasticity and history remap) builds on `prototype/lr1d.py` and the WP 1.4 models. The C++ ports (2.2, 2.5a, 2.5b) use WP 1.4 results as reference.
3. WP 0.4: fast OpenFOAM setup for fresh cloud sessions (a from-source build takes about an hour).

## Known issues

- WP 1.1 prototype: smooth-wave order is limiter-bound at about 1.8. It passes at the test
  setting (0.25 c streaming) but gives 1.78 at 0.5 c. Only matters if later packages need
  a sharper order; an extremum-preserving limiter would fix it.
- Mohr-Coulomb prototype has no separate tension cut-off below the apex (decision 0004).
