# WP 2.2: Constitutive library interface in C++, elastic and J2

State: not started
Branch: `wp-2.2-constitutive-library`

## Goal

`src/celConstitutive`: the C++ material-model library both the Eulerian and Lagrangian
sides will call, with linear elastic and von Mises (J2) ported from the Python reference
and matching it to round-off. Drucker-Prager and Mohr-Coulomb follow in 2.5a and 2.5b
using the interface built here.

## Read first

- `docs/decisions/0004-constitutive-prototype-conventions.md` (conventions to keep)
- `docs/constitutive-models.md`
- `prototype/constitutive/models.py` (`linear_elastic`, `von_mises`), `driver.py`,
  `check_models.py`
- `CLAUDE.md` architecture section

## Build

Proposed design, to be confirmed in a decision note (change it if you find a good reason):

- **Kernel in plain C++ (no OpenFOAM headers):** `src/celConstitutive/kernel/`. Each model
  is a function or small class with the same contract as the Python reference: strain
  increment, old stress, history, parameters in; new stress and history out. Stress and
  strain as 6 doubles in the Voigt order of decision 0004. History as a fixed-size struct
  (`epsP[6]`, `epsPEq`, `apex`), because it will be stored per cell and remapped.
  Reason: it unit-tests in seconds with plain `g++`, with no OpenFOAM installed, and the
  FE library can call it directly.
- **Thin OpenFOAM wrapper:** `src/celConstitutive/` with a `constitutiveModel` base class
  using OpenFOAM run-time selection (`type vonMises;` in `constant/materialProperties`),
  converting `symmTensor` to and from the kernel's arrays (mind the component order:
  OpenFOAM `symmTensor` is xx, xy, xz, yy, yz, zz). Built by `Allwmake` into
  `libcelConstitutive.so`.
- **Reference data:** a script `tests/gen_constitutive_reference.py` that runs the Python
  models on fixed strain paths (uniaxial strain load-unload, triaxial, simple shear, 1000
  random increments with a fixed seed) and writes them as text; a small C++ test program
  reads the same paths and compares.
- `tests/test_constitutive_cpp.sh`: builds the kernel test with `g++`, runs the comparison.
  If OpenFOAM is available it also builds the library; CI covers that part.

## Done when

- `tests/test_constitutive_cpp.sh` passes: C++ elastic and von Mises stresses and history
  match the Python reference on all paths to a relative 1e-10.
- `libcelConstitutive.so` builds in CI and `celFoam` links it and prints the selected model
  for the `smokeBox` tutorial. Existing tests still pass unchanged.

## Out of scope

- Drucker-Prager (2.5a), Mohr-Coulomb (2.5b), calling the model in the solver's time loop
  (2.3), objective stress rates (2.3).

## Note

OpenFOAM is not installed in a fresh cloud session (WP 0.4 fixes that). Do the kernel and
its tests locally with `g++`; rely on CI for the OpenFOAM wrapper build rather than
spending an hour compiling OpenFOAM.

## Log

(none yet)
