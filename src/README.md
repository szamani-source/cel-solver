# src

Libraries are added here by later work packages, one folder per library:

- `celConstitutive/` (WP 2.2): material models. `kernel/` is plain C++ (no OpenFOAM, tests with
  g++); `constitutiveModel/` and `constitutiveModels/` are the OpenFOAM wrapper with run-time
  selection. Layout and conventions: decision 0006.
- `celEulerian/` (WP 2.x, 3.x): Eulerian fields, Lagrange step, remap, volume fractions.
- `celLagrangian/` (WP 4.x, 5.x): rigid bodies and explicit hex8 FE.
- `celCoupling/` (WP 4.x): Eulerian-Lagrangian contact search and penalty forces.
