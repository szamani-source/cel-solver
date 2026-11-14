# WP 2.5b: Mohr-Coulomb in C++

State: not started
Branch: `wp-2.5b-mohr-coulomb`

## Goal

Mohr-Coulomb in `src/celConstitutive`, matching the Python reference to round-off and
selectable at run time (`type mohrCoulomb;`). This is the main model for sand and other
bulk solids, so its corner and apex handling must be exact.

## Read first

- `docs/decisions/0004-constitutive-prototype-conventions.md`, `0006-constitutive-library-layout.md`
- `prototype/constitutive/models.py`: `_mc_principal_planes`, `mohr_coulomb`, `mc_yield`
- `src/celConstitutive/kernel/celVonMises.{H,C}` and
  `constitutiveModels/vonMises/` as the pattern to copy
- `tests/gen_constitutive_reference.py`, `tests/constitutive/compareReference.C`,
  `tests/test_constitutive_cpp.sh`, `docs/wp/WP-2.2.md` log

## Build

- `kernel/celMohrCoulomb.{H,C}`: return mapping in sorted principal stress space as in
  the Python model: main plane, then the compression or extension edge, then the apex;
  non-associated flow. Needs a symmetric 3x3 eigen-decomposition in plain C++ (no
  OpenFOAM, no external library): write a small Jacobi or analytic solver in the kernel
  and unit-test it, including repeated eigenvalues (hydrostatic and triaxial states).
- `constitutiveModels/mohrCoulomb/mohrCoulombModel.{H,C}` wrapper; add to `Make/files`
  and to `makeModel`.
- Reference paths added to the generator for `mohr_coulomb`: triaxial compression and
  extension (edge returns), simple shear, uniaxial strain load-unload, 1000 random
  increments (fixed seed) that include edge and apex returns, for at least two parameter
  sets: frictional-cohesive (phi 30, psi 5, c 10 kPa) and cohesionless sand (phi 35,
  psi 0, c 0).
- Tension cut-off: the Python reference has none (decision 0004). Do not add one here
  unless the reference gets it first; note in the log whether it is still wanted.

## Done when

- `tests/test_constitutive_cpp.sh` passes with the new paths: stress and history match
  the Python reference to a relative 1e-10 on every path (if repeated eigenvalues make
  the principal directions non-unique, compare stresses, which are unique; record any
  path that needs a looser tolerance and why). Existing paths unchanged.
- Library and wrapper check build and pass in CI.

## Out of scope

- Drucker-Prager (2.5a, running in parallel: expect small merge conflicts in `Make/files`,
  `makeModel`, the generator and `docs/STATUS.md`; whoever merges second resolves them).
- Hardening, softening. Use in the solver loop.

## Log

(none yet)
