# Status

Last updated: 2026-11-14

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 1 prototypes, with Phase 0 tooling (0.4) and the first Phase 2 library (2.2) running alongside.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Done | Repository szamani-source/cel-solver; CI builds celFoam in opencfd/openfoam-dev:2512 and test_smokeBox passes (2026-11-14) |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Draft | `docs/spec/target-problems.md` waits for Sia's answers |
| 0.4 | Fast OpenFOAM setup for fresh cloud sessions | Next | Brief: `docs/wp/WP-0.4.md` |
| 1.1 | 1D elastic Lagrange-plus-remap prototype | Done | `prototype/lr1d.py`, `tests/test_proto_1d_elastic.sh` passes (order 1.88/1.83, monotone, conservative). Cell-centred instead of staggered: decision 0005. PR #3 merged (2026-11-14) |
| 1.2 | 1D plasticity and history remap | Next | Brief: `docs/wp/WP-1.2.md` |
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Done | `prototype/constitutive/`; `tests/test_proto_constitutive.sh` passes in CI; PR #2 merged (2026-11-14) |
| 2.2 | Constitutive library interface in C++, elastic and J2 | In progress | `src/celConstitutive/` (plain C++ kernel + OpenFOAM wrapper, decision 0006); `tests/test_constitutive_cpp.sh` kernel part passes locally (error 7e-16 vs Python); wrapper part waits on CI. PR open |

## Waiting on Sia

- Answers to the open questions in `docs/spec/target-problems.md` (closes 0.3). Not
  blocking 1.2, 0.4 or 2.2; needed before the Phase 2 benchmarks (2.6).
- Confirm GPL-3.0 (decision 0003). Not blocking development; needed before the code is
  shared with anyone.

## Next

Started 2026-11-14, one thread each, all independent of each other:

1. WP 1.2: 1D plasticity and history remap (`docs/wp/WP-1.2.md`).
2. WP 0.4: fast OpenFOAM setup (`docs/wp/WP-0.4.md`).
3. WP 2.2: constitutive library in C++, elastic and J2 (`docs/wp/WP-2.2.md`): PR open.

After those: 1.3 (2D prototype and formulation decision) follows 1.2; 2.5a and 2.5b
(Drucker-Prager and Mohr-Coulomb in C++) follow 2.2; 2.1 and 2.3 follow 0.4 and 2.2.

## Known issues

- WP 1.1 prototype: smooth-wave order is limiter-bound at about 1.8. It passes at the test
  setting (0.25 c streaming) but gives 1.78 at 0.5 c. Only matters if later packages need
  a sharper order; an extremum-preserving limiter would fix it.
- Mohr-Coulomb prototype has no separate tension cut-off below the apex (decision 0004).
