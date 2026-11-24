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
- `celFluid/` (WP F.1, then F.2 to F.9, version 2): gas in the pores and the void
  (decisions 0125, 0126, 0127). `poreGas/` reads `constant/gasProperties`, writes the gas
  fractions `phiS`, `phiG`, the permeability `k`, and solves the pore pressure `p` (gauge)
  and flux `qGas` each step (F.2); `gasDragModels/` holds the run-time selected drag laws
  (`ergunWenYu`, the default; more in F.7).
