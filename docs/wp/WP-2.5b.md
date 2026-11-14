# WP 2.5b: Mohr-Coulomb in C++

State: done (PR open, merges when CI is green)
Branch: `claude/project-thread-ie0bvq` (thread branch; the brief planned `wp-2.5b-mohr-coulomb`)

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

Newest first.

### 2026-11-14 (session 1)

Done:
- `kernel/celEigen.{H,C}`: symmetric 3x3 eigen-solver by cyclic Jacobi, eigenvalues
  sorted descending (decision 0011). `tests/constitutive/eigenCheck.C`: distinct,
  hydrostatic, both triaxial repeats, gap 1e-9, zero, pure shear, wide range, tiny values
  and 10000 random spectra (a quarter with two equal, a quarter with three equal); worst
  error 5.9e-15 relative, tolerance 1e-13.
- `kernel/celMohrCoulomb.{H,C}`: return mapping copied step for step from `mohr_coulomb`
  in `models.py`. `updateReturn` also reports the return type (for the test).
  `plasticCorrector` added to the kernel with the same signature and body as WP 2.5a's.
- `constitutiveModels/mohrCoulomb/` wrapper (`type mohrCoulomb; E nu phi c psi`),
  `Make/files`, `makeModel`.
- Generator: model names in the reference file are now the C++ names (`PY_MODELS` maps
  them to the Python functions, as in WP 2.5a). Ten new paths: uniaxial strain (load,
  unload, pull into tension), triaxial compression, triaxial extension, simple shear and
  the 1000 random increments, for phi 30 / psi 5 / c 10 kPa and phi 35 / psi 0 / c 0.
  Existing paths unchanged (same numbers).
- `compareReference.C` counts returns on Mohr-Coulomb paths and requires: compression
  edge on triaxial compression, extension edge on triaxial extension, apex on uniaxial
  strain, and an edge and an apex return on the random paths.
- `wrapperCheck.C`: selects `mohrCoulomb` and checks an isochoric axial compression lands
  on the compression edge (sxx = syy) with f = 0.

Result: `tests/test_constitutive_cpp.sh` passes locally with OpenFOAM loaded. All 25
paths match; Mohr-Coulomb worst errors 5.4e-15 (stress) and 1.1e-14 (history), tolerance
1e-10; no path needed a looser tolerance. Returns seen: random cohesive 98 plane, 6
compression edge, 2 extension edge, 7 apex; random sand 75 / 4 / 6 / 14.

Learned:
- Uniaxial strain from 100 kPa never yields in loading with nu = 0.3 (elastic K0 0.43 is
  above the active ratio), so the uniaxial path pulls into tension to reach the apex.
- Tension cut-off: still wanted for cohesive soils (cohesive MC carries tension up to
  c cot(phi)); split out as WP 2.5c, which needs a Python reference first.

Merged main after WP 2.5a (PR #7) landed: kept both models in `Make/files`, `makeModel`,
the generator, `compareReference.C` (DP apex check and MC return check) and
`wrapperCheck.C`. `plasticCorrector` was already identical. After the merge all 39
reference paths, the eigen check and the wrapper checks pass locally.

Next: CI green, then merge.
