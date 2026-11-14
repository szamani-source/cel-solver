# 0006 Constitutive library layout: plain C++ kernel plus thin OpenFOAM wrapper (WP 2.2)

Date: 2026-11-14. Status: proposed (WP 2.2), confirm by merging the PR.

## Context

The material models are shared by the Eulerian FV solid and the Lagrangian FE library,
must match the Python reference of WP 1.4 to round-off, and must be testable in a fresh
cloud session where OpenFOAM is not installed.

## Decision

- `src/celConstitutive/kernel/`: plain C++17, no OpenFOAM headers, namespace
  `cel::constitutive`. Abstract `Model` with one call
  `update(deps, sigOld, hOld, sigNew, hNew)`; parameters fixed at construction from a
  name-to-value map; `makeModel(name, params)` factory. Stress and strain are `double[6]`
  in the Voigt order of decision 0004 (engineering shear strain). History is a fixed-size
  struct `History { epsP[6], epsPEq, apex }`, so it can be stored per cell and remapped.
- `src/celConstitutive/constitutiveModel/`: OpenFOAM class `constitutiveModel` with
  run-time selection (`type` keyword), one small derived class per model that only lists
  its parameter names. It takes `symmTensor` and converts the order
  (OpenFOAM xx, xy, xz, yy, yz, zz to Voigt xx, yy, zz, xy, yz, zx). Its strain input is
  the tensor strain increment (e.g. `symm(gradU)*dt`); the wrapper doubles the shear
  terms for the kernel. Both parts build into `libcelConstitutive.so`.
- Kernel tests replay strain increments written by the Python models
  (`tests/gen_constitutive_reference.py`), so the Python code stays the single reference
  and the reference data is regenerated each run instead of being committed.
- One virtual call per point per step. Its cost is small next to the rest of the update;
  revisit only if profiling says so.

## Consequences

- Kernel changes are checked in seconds with `g++` alone; only the wrapper needs
  OpenFOAM (CI, or a session with WP 0.4's setup).
- The error measure is relative to the largest value on each path (stress or plastic
  strain), not per component, because components cross zero.
- `apex` is a `bool`. When history is remapped (WP 2.4) it will need a rule (for example
  majority or "any donor at the apex"); decide then.
- Drucker-Prager (2.5a) and Mohr-Coulomb (2.5b) add one kernel class, one wrapper class
  and their reference paths to the generator; nothing else changes.
