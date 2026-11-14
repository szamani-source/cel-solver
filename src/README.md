# src

Libraries are added here by later work packages, one folder per library:

- `celConstitutive/` (WP 2.2): material models. `kernel/` is plain C++ (no OpenFOAM, tests with
  g++); `constitutiveModel/` and `constitutiveModels/` are the OpenFOAM wrapper with run-time
  selection. Layout and conventions: decision 0006.
- `celEulerian/` (WP 2.1, then 2.x, 3.x): Eulerian fields, Lagrange step, remap, volume fractions.
  `eulerianSolid/` (WP 2.1) owns the fields `rho`, `U`, `sigma`, `epsP`, `epsPEq` and the
  constitutive model, with one method per step of the explicit cycle (decision 0010).
  The Lagrange step (WP 2.3) is in `eulerianSolidLagrange.C` (decision 0016).
- `celLagrangian/` (WP 4.x, 5.x): rigid bodies and explicit hex8 FE.
- `celCoupling/` (WP 4.x): Eulerian-Lagrangian contact search and penalty forces.
