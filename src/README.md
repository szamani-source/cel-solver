# src

Libraries are added here by later work packages, one folder per library:

- `celConstitutive/` (WP 2.2): material models. `kernel/` is plain C++ (no OpenFOAM, tests with
  g++); `constitutiveModel/` and `constitutiveModels/` are the OpenFOAM wrapper with run-time
  selection. Layout and conventions: decision 0006.
- `celEulerian/` (WP 2.1, then 2.x, 3.x): Eulerian fields, Lagrange step, remap, volume fractions.
  `eulerianSolid/` (WP 2.1) owns the fields `rho`, `U`, `sigma`, `epsP`, `epsPEq` and the
  constitutive model, with one method per step of the explicit cycle (decision 0010).
  The Lagrange step (WP 2.3) is in `eulerianSolidLagrange.C` (decision 0018).
- `celLagrangian/` (WP 4.1, then 5.x): rigid bodies from STL with prescribed motion
  (`rigidBody`, `rigidBodies`, decision 0041), later explicit hex8 FE.
- `celCoupling/` (WP 4.2, then 4.3a, 4.3b): `bodyCells`, exact cells and faces covered by
  each rigid body (decision 0043); the cut-cell wall coupling comes in 4.3a (decision 0042).
