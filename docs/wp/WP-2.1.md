# WP 2.1: Eulerian library skeleton and solver fields

State: done (pending CI and merge)
Branch: `wp-2.1-eulerian-skeleton`

## Goal

`src/celEulerian` exists as a library holding the Eulerian state, and `celFoam` runs the
explicit cycle through it with the physics steps still empty. Later packages (2.3 stress
update, 2.4 remap) then fill in one function each instead of reshaping the solver.

## Read first

- `CLAUDE.md` (architecture, "Setting up OpenFOAM in a new session": run
  `tools/setup_openfoam.sh` first)
- `applications/solvers/celFoam/`, `src/celConstitutive/constitutiveModel/constitutiveModel.H`
- `docs/decisions/0005`, `0006`, `0008`

## Build

- `src/celEulerian/eulerianSolid/eulerianSolid.{H,C}`: a class owning the cell-centred
  fields `rho` (volScalarField), `U` (volVectorField), `sigma` (volSymmTensorField),
  `epsP` (volSymmTensorField), `epsPEq` (volScalarField), all read from the time directory
  if present, else created with defaults (`rho` from `materialProperties`, the rest zero),
  and written at output times. It holds the `constitutiveModel`.
- Methods, each a clearly named stub that later packages fill:
  `stableDeltaT()` (implemented: `maxCo * min(cbrt(V)) / max wave speed`, using the
  model's wave modulus and the cell density), `lagrangeStep(dt)` (WP 2.3),
  `remap(dt)` (WP 2.4), plus `totalMass()`, `totalMomentum()` diagnostics (implemented,
  parallel-safe with `gSum`).
- `celFoam.C` becomes a thin driver: construct `eulerianSolid`, loop
  `setDeltaT(stableDeltaT())`, `lagrangeStep`, `remap`, write; print mass and momentum
  each write. The time step is recomputed every cycle (it will change once stress and
  density change).
- `tutorials/smokeBox` updated: fields written, a non-uniform initial `rho` or `U` allowed.
- `tests/test_eulerian_skeleton.sh`; update `Allwmake`, `Allwclean`, `src/README.md`.
- Do not choose the Lagrange-step discretisation here: WP 1.3b decides it.

## Done when

`tests/test_eulerian_skeleton.sh` passes in CI:

- `celFoam` runs `smokeBox` to the end; the written time directory contains `rho`, `U`,
  `sigma`, `epsP`, `epsPEq`.
- Printed time step equals the hand calculation (as `test_smokeBox`), and total mass
  equals `rho * volume` to 1e-12 relative.
- The same case decomposed on 2 processors (`decomposePar`, `mpirun -np 2 celFoam -parallel`;
  allow root in the container) prints the same time step, mass and momentum.
- Existing tests still pass unchanged (`test_smokeBox` may need its expected output
  adjusted only if the printed text moved; never loosen a numeric check).

## Out of scope

- Any physics in `lagrangeStep` or `remap`. More than one material (Phase 3).

## Log

2026-11-14 (session 1):

Done:
- `src/celEulerian/eulerianSolid/eulerianSolid.{H,C}`: class `eulerianSolid` owning
  `materialProperties`, the constitutive model and the fields `rho`, `U`, `sigma`,
  `epsP`, `epsPEq` (read if present, else defaults, written at output times).
  `stableDeltaT()`, `totalMass()`, `totalMomentum()` implemented (parallel-safe);
  `lagrangeStep(dt)` and `remap(dt)` are empty stubs for WP 2.3 and 2.4.
- `celFoam` is a thin driver: each cycle `setDeltaT(stableDeltaT())`, `lagrangeStep`,
  `remap`, write; prints total mass and momentum (16 digits) at start and each write.
- `tutorials/smokeBox`: fields moved to `0.orig/` (adds `rho`), `setFieldsDict`
  (half the box rho 2000, U (1 0 0)), `decomposeParDict` (2 subdomains), `Allrun`
  restores `0/` and runs `setFields`.
- `tests/test_eulerian_skeleton.sh`: (A) no `0/`, all defaults; (B) setFields serial;
  (C) same case on 2 processors. Checks End, five fields written, start deltaT equals
  the hand value to 1e-12, mass 1800 / 1900 and momentum (1000 0 0) to 1e-12, and
  serial and parallel print identical time steps.
- Decision 0009 records the choices (per-cell time step formula, `0.orig`, defaults
  with `zeroGradient` patches, `apex` not stored).

Local result: all six tests pass (OpenFOAM from `tools/setup_openfoam.sh`). Mass errors
1.7e-14 (serial A) and 4e-15 to 8e-15 relative (B, C).

Next: WP 2.3 fills `lagrangeStep` once WP 1.3b has chosen the discretisation; WP 2.4
fills `remap`.
