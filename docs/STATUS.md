# Status

Last updated: 2026-11-14

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 1 (2D prototype) and Phase 2 (C++ material models, Eulerian library skeleton) in parallel.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Done | Repository szamani-source/cel-solver; CI builds celFoam in opencfd/openfoam-dev:2512 and test_smokeBox passes (2026-11-14) |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Draft | `docs/spec/target-problems.md` waits for Sia's answers |
| 0.4 | Fast OpenFOAM setup for fresh cloud sessions | Done | `tools/setup_openfoam.sh` installs the CI image's OpenFOAM from registry layers; setup + build + tests in 66 s (decision 0007). PR #5 merged (2026-11-14) |
| 1.1 | 1D elastic Lagrange-plus-remap prototype | Done | `prototype/lr1d.py`, `tests/test_proto_1d_elastic.sh` passes (order 1.88/1.83, monotone, conservative). Cell-centred instead of staggered: decision 0005. PR #3 merged (2026-11-14) |
| 1.2 | 1D plasticity and history remap | Done | `prototype/lr1d_plastic.py`, `tests/test_proto_1d_plastic.sh` passes (wave speeds and HEL within 0.2 %, history edges 2 cells). Decision 0008 accepted. PR #6 merged (2026-11-14) |
| 1.3a | 2D plane-strain prototype, elastic waves | Next | Brief: `docs/wp/WP-1.3a.md` |
| 1.3b | 2D rotation, plasticity and formulation decision | Blocked on 1.3a | Brief: `docs/wp/WP-1.3b.md` |
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Done | `prototype/constitutive/`; `tests/test_proto_constitutive.sh` passes in CI; PR #2 merged (2026-11-14) |
| 2.1 | Eulerian library skeleton and solver fields | Done | `src/celEulerian` (class `eulerianSolid`), celFoam is a thin driver; `tests/test_eulerian_skeleton.sh` passes serial and on 2 processors (decision 0009). PR #8 |
| 2.2 | Constitutive library interface in C++, elastic and J2 | Done | `src/celConstitutive/` (plain C++ kernel + OpenFOAM wrapper, decision 0006); `tests/test_constitutive_cpp.sh` passes in CI (kernel matches Python to 3e-15; wrapper and celFoam checks pass). PR #4 merged (2026-11-14) |
| 2.5a | Drucker-Prager in C++ | Next | Brief: `docs/wp/WP-2.5a.md` |
| 2.5b | Mohr-Coulomb in C++ | Next | Brief: `docs/wp/WP-2.5b.md` |

## Standing instruction from Sia (2026-11-14)

Proceed without asking for permission, and merge each pull request automatically once
CI is green. Sia is away; report results, do not wait for review. Still ask before
anything that cannot be undone outside this repository.

## Waiting on Sia

- Answers to the open questions in `docs/spec/target-problems.md` (closes 0.3). Not
  blocking current work; needed before the Phase 2 benchmarks (2.6).
- Optional: paste the setup line from `docs/wp/WP-0.4.md` into the cloud environment's
  setup script so every session starts with OpenFOAM installed.
- Confirm GPL-3.0 (decision 0003). Not blocking development; needed before the code is
  shared with anyone.

## Next

Started 2026-11-14, one thread each, independent of each other:

1. WP 1.3a: 2D plane-strain prototype, elastic waves (`docs/wp/WP-1.3a.md`).
2. WP 2.5a: Drucker-Prager in C++ (`docs/wp/WP-2.5a.md`).
3. WP 2.5b: Mohr-Coulomb in C++ (`docs/wp/WP-2.5b.md`).
4. WP 2.1: Eulerian library skeleton and solver fields (`docs/wp/WP-2.1.md`): done.

After those: 1.3b follows 1.3a. 2.3 (stress update with objective rate) and 2.4 (remap)
follow the 1.3b formulation decision (2.1 is done; each fills one `eulerianSolid` method); briefs not written yet.

## Known issues

- WP 1.1 prototype: smooth-wave order is limiter-bound at about 1.8. It passes at the test
  setting (0.25 c streaming) but gives 1.78 at 0.5 c. Only matters if later packages need
  a sharper order; an extremum-preserving limiter would fix it.
- WP 1.2 prototype: an impulsive start (velocity jump at a face) overshoots by 2.7
  percent of the jump at the impact plane for the first 5 to 10 cells of wave travel
  (1.7 percent elastic). Matters for contact (WP 4.x); revisit in WP 1.3 (decision 0008).
- WP 1.2 prototype: history is remapped with superbee, stress with MC (decision 0008).
- Mohr-Coulomb prototype has no separate tension cut-off below the apex (decision 0004).
