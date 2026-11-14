# WP 2.5a: Drucker-Prager in C++

State: not started
Branch: `wp-2.5a-drucker-prager`

## Goal

Drucker-Prager in `src/celConstitutive`, matching the Python reference to round-off and
selectable at run time (`type druckerPrager;`).

## Read first

- `docs/decisions/0004-constitutive-prototype-conventions.md`, `0006-constitutive-library-layout.md`
- `prototype/constitutive/models.py`: `dp_coefficients`, `drucker_prager`, `dp_yield`
- `src/celConstitutive/kernel/celVonMises.{H,C}` and
  `constitutiveModels/vonMises/` as the pattern to copy
- `tests/gen_constitutive_reference.py`, `tests/constitutive/compareReference.C`,
  `tests/test_constitutive_cpp.sh`, `docs/wp/WP-2.2.md` log

## Build

- `kernel/celDruckerPrager.{H,C}`: cone return and apex return as in the Python model;
  non-associated flow (psi separate from phi); the three fits to Mohr-Coulomb
  (`compression`, `extension`, `plane_strain`). The kernel's `Parameters` map holds
  numbers only: choose how to pass `fit` (for example a separate constructor argument or
  string option map) and record it in the brief log; keep it usable from the dictionary
  as a word (`fit compression;`).
- `constitutiveModels/druckerPrager/druckerPragerModel.{H,C}` wrapper; add to `Make/files`
  and to `makeModel`.
- Reference paths added to the generator for `drucker_prager`: triaxial compression and
  extension, simple shear at constant normal stress, uniaxial strain load-unload, 1000
  random increments (fixed seed) that include apex returns, for at least two parameter
  sets: frictional-cohesive (phi 30, psi 5, c 10 kPa) and cohesionless (phi 35, psi 0,
  c 0), and for each of the three fits on one path.
- `tutorials/smokeBox` variant or extra check showing `celFoam` selects `druckerPrager`.

## Done when

- `tests/test_constitutive_cpp.sh` passes with the new paths: stress and history match
  the Python reference to a relative 1e-10 on every path; existing paths unchanged.
- Library and wrapper check build and pass in CI.

## Out of scope

- Mohr-Coulomb (2.5b, running in parallel: expect small merge conflicts in `Make/files`,
  `makeModel`, the generator and `docs/STATUS.md`; whoever merges second resolves them).
- Hardening, softening, caps. Use in the solver loop.

## Log

(none yet)
