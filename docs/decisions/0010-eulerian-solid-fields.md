# 0010 Eulerian solid: one class owns the state, fields default when absent (WP 2.1)

Date: 2026-11-14. Status: proposed (WP 2.1), confirm by merging the PR.

## Context

WP 2.1 moves the Eulerian state out of `celFoam` into `src/celEulerian` so that WP 2.3
(stress update) and WP 2.4 (remap) each fill one method. Several small choices were
needed that later packages might question.

## Decision

- One class, `Foam::eulerianSolid`, owns `materialProperties`, the `constitutiveModel`
  and the cell-centred fields `rho`, `U`, `sigma`, `epsP`, `epsPEq`. `celFoam` only calls
  `stableDeltaT()`, `lagrangeStep(dt)`, `remap(dt)` and the diagnostics.
- History is stored as two fields (`epsP` symmTensor with tensor shear components,
  `epsPEq` scalar), matching the wrapper's convention (decision 0006). The kernel's
  `apex` flag is not stored: it is an output of the last increment, not state that
  must survive a step. If WP 2.3 needs it per cell, add a field then.
- Every field is `READ_IF_PRESENT`; absent fields are created (rho from
  `materialProperties`, the rest zero) with `zeroGradient` patches. Boundary conditions
  that mean something physically (walls, symmetry) are left to WP 2.3.
- Stable time step: `maxCo * min over cells of cbrt(V)/sqrt(M/rho)` with
  `M = K + 4G/3` from the model. The brief wrote `min(cbrt(V)) / max c`; the per-cell
  form equals it on a uniform mesh and is less conservative on graded meshes.
  Recomputed every cycle.
- The tutorial keeps its fields in `0.orig/` (OpenFOAM convention) so `Allrun` can run
  `setFields` without changing committed files. Tests that copy the tutorial without
  restoring `0/` thereby check the all-defaults path.
- Mass, momentum and the start time step are printed with 16 digits so tests can
  compare runs to 1e-12.

## Consequences

- `test_smokeBox` and `test_constitutive_cpp` were not changed; they now run with every
  field defaulted (no `0/` directory), which is a valid case.
- Parallel decomposition needs no special code yet: diagnostics use `gSum`/`gMin`.
  WP 2.3/2.4 must keep processor patches in mind when they add face loops.
