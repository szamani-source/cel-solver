# Status

Last updated: 2026-11-14

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 2 (the C++ Eulerian solver: 2.3 done, 2.4 running) and Phase 3 (silo discharge, Sia's primary target, decision 0016) in parallel.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Done | Repository szamani-source/cel-solver; CI builds celFoam in opencfd/openfoam-dev:2512 and test_smokeBox passes (2026-11-14) |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Done | Sia named the primary target (2026-11-14): 2D silo discharge, cohesionless, phi 35 deg. `docs/spec/target-problems.md` (T0). Tank dimensions are placeholders |
| 0.4 | Fast OpenFOAM setup for fresh cloud sessions | Done | `tools/setup_openfoam.sh` installs the CI image's OpenFOAM from registry layers; setup + build + tests in 66 s (decision 0007). PR #5 merged (2026-11-14) |
| 1.1 | 1D elastic Lagrange-plus-remap prototype | Done | `prototype/lr1d.py`, `tests/test_proto_1d_elastic.sh` passes (order 1.88/1.83, monotone, conservative). Cell-centred instead of staggered: decision 0005. PR #3 merged (2026-11-14) |
| 1.2 | 1D plasticity and history remap | Done | `prototype/lr1d_plastic.py`, `tests/test_proto_1d_plastic.sh` passes (wave speeds and HEL within 0.2 %, history edges 2 cells). Decision 0008 accepted. PR #6 merged (2026-11-14) |
| 1.3a | 2D plane-strain prototype, elastic waves | Done | `prototype/lr2d.py` (symmetric dimensional splitting), `tests/test_proto_2d_elastic.sh` passes: P and S waves at 0 and 45 degrees, orders 1.77 to 1.87, speeds within 0.002 %; matches `lr1d.py` bit-for-bit; mirror-exact. PR #10 merged (2026-11-14). The 45 degree overshoot was a jagged test step (decision 0013) |
| 1.3b | 2D rotation, plasticity and formulation decision | Done | `prototype/lr2d_plastic.py` (Jaumann via Hughes-Winget, WP 1.4 models vectorised in `constitutive/batch.py`); `tests/test_proto_2d_rotation.sh` passes (disc 0.03 %, simple shear 3e-7, 1D match 1e-11). **Formulation: `docs/decisions/0013-eulerian-formulation.md`**: cell-centred, unsplit in C++, Jaumann/Hughes-Winget, accept start-up overshoot, MC + superbee for history, no nodal solver |
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Done | `prototype/constitutive/`; `tests/test_proto_constitutive.sh` passes in CI; PR #2 merged (2026-11-14) |
| 2.1 | Eulerian library skeleton and solver fields | Done | `src/celEulerian` (class `eulerianSolid`), celFoam is a thin driver; `tests/test_eulerian_skeleton.sh` passes serial and on 2 processors (decision 0010). PR #8 |
| 2.2 | Constitutive library interface in C++, elastic and J2 | Done | `src/celConstitutive/` (plain C++ kernel + OpenFOAM wrapper, decision 0006); `tests/test_constitutive_cpp.sh` passes in CI (kernel matches Python to 3e-15; wrapper and celFoam checks pass). PR #4 merged (2026-11-14) |
| 2.3 | Lagrange step in C++ | Done | `eulerianSolidLagrange.C`: MUSCL-Hancock face states, acoustic face Riemann solver, Hughes-Winget, boundaries from the U patch type, default maxCo 0.4 (decision 0018). `tests/test_lagrange_step.sh` passes: wave orders 1.9 to 2.3, rotation 3e-13, shear 3e-7, reflections 1e-4, 2/4 processors 3e-15 (2026-11-14) |
| 2.4 | Remap in C++, single material | Next | Brief: `docs/wp/WP-2.4.md` |
| 2.4b | Lagrange step plus remap, integration checks | Blocked on 2.3 and 2.4 | Brief: `docs/wp/WP-2.4b.md` |
| 2.5a | Drucker-Prager in C++ | Done | `kernel/celDruckerPrager`, wrapper `druckerPrager` (`fit` word option, decision 0009); 14 new reference paths match Python to 2e-15 with apex returns; celFoam selects it. Merged on green CI (2026-11-14) |
| 2.5b | Mohr-Coulomb in C++ | Done | `kernel/celMohrCoulomb`, Jacobi eigen-solver `kernel/celEigen` (decision 0011); `tests/test_constitutive_cpp.sh` matches Python on 10 MC paths to 1e-14, edge and apex returns required. Tension cut-off split to 2.5c |
| 2.5c | Mohr-Coulomb tension cut-off (Python, then C++) | Done | Optional `sigma_t` / `sigmaT` (Rankine planes, active-set return, decision 0015). Python checks pass (strengths, caps, unloading, 4 x 10,000 tension-biased increments, associated return = independent closest point to 1e-15); C++ matches Python on 10 new paths to 2e-14 with every tension return exercised; celFoam reads `sigmaT`. Merged on green CI (2026-11-14) |
| 3.0 | Silo benchmark definition and reference values | Next | Brief: `docs/wp/WP-3.0.md` |
| 3.1p | Prototype gravity, walls and settling | Next | Brief: `docs/wp/WP-3.1p.md` |
| 3.2p | Prototype free surface (material plus void) | Blocked on 3.1p | Brief: `docs/wp/WP-3.2p.md` |
| 3.3p | Prototype silo discharge (target T0) | Blocked on 3.2p | Brief: `docs/wp/WP-3.3p.md` |

## Standing instruction from Sia (2026-11-14)

Proceed without asking for permission, and merge each pull request automatically once
CI is green. Sia is away; report results, do not wait for review. Still ask before
anything that cannot be undone outside this repository.

## Waiting on Sia

- Optional: real tank dimensions and material data for the silo problem (placeholders
  are in `docs/spec/target-problems.md`).
- Optional: paste the setup line from `docs/wp/WP-0.4.md` into the cloud environment's
  setup script so every session starts with OpenFOAM installed.
- Confirm GPL-3.0 (decision 0003). Needed before the code is shared with anyone.
- Abaqus or LS-DYNA access for reference runs? Not needed for the silo problem.

## Next

WP 2.3 (Lagrange step in C++) is done (decision 0018). Running: WP 2.4 (remap in C++).
After it merges: WP 2.4b (`docs/wp/WP-2.4b.md`).

Silo track (decision 0016), one thread each:

1. WP 3.0: reference values (`docs/wp/WP-3.0.md`). Start now.
2. WP 3.1p: prototype gravity, walls, settling (`docs/wp/WP-3.1p.md`). Start now.
3. WP 3.2p: prototype free surface (`docs/wp/WP-3.2p.md`). After 3.1p.
4. WP 3.3p: prototype silo discharge, the first end-to-end result for Sia
   (`docs/wp/WP-3.3p.md`). After 3.2p.

Then the C++ counterparts 3.1 to 3.4 (after 2.4b; briefs not written yet) and 3.5
(cohesion).

## Known issues

- WP 1.1 prototype: smooth-wave order is limiter-bound at about 1.8. It passes at the test
  setting (0.25 c streaming) but gives 1.78 at 0.5 c. Only matters if later packages need
  a sharper order; an extremum-preserving limiter would fix it.
- WP 1.2 prototype: an impulsive start (velocity jump at a face) overshoots by 2.7
  percent of the jump at the impact plane for the first 5 to 10 cells of wave travel
  (1.7 percent elastic). Matters for contact (WP 4.x); revisit in WP 1.3 (decision 0008).
- WP 1.2 prototype: history is remapped with superbee, stress with MC (decision 0008).
- Mohr-Coulomb tension cut-off (WP 2.5c) is optional: without `sigmaT` cohesive
  material still carries tension up to c cot(phi). The vectorised 2D prototype
  (`prototype/constitutive/batch.py`) and Drucker-Prager have no cut-off.
- 2D prototype: a step at 45 degrees overshoots 2.7 percent during its first cell of
  travel, then below 0.3 percent (decision 0013; the 13 percent in decision 0012 was a
  jagged test step). Like the 1D impact start-up, accepted; contact (WP 4.3) re-measures.
- C++ Lagrange step (WP 2.3): stable Courant number 0.65 in 2D, 0.45 to 0.50 in 3D
  on uniform meshes (default maxCo 0.4); not measured on graded or distorted meshes. The
  rotating-disc check needs the remap and moves to WP 2.4b (decision 0018).
- `stableDeltaT` uses cbrt(V): 2D one-cell-thick meshes must be as thick as their
  in-plane cells (decision 0018, consequences).
- Until the remap (WP 2.4) is merged, celFoam runs the Lagrange step with no remap:
  the printed mass drifts with the Lagrangian density.
- WP 1.3a test takes about 5 minutes (400 x 400 meshes for the 45 degree convergence runs);
  WP 1.3b test about 1.5 minutes.
