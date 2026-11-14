# 0011 Mohr-Coulomb in C++: Jacobi eigen-solver, return labels, no tension cut-off (WP 2.5b)

Date: 2026-11-14. Status: proposed (WP 2.5b), confirm by merging the PR.

## Context

The Mohr-Coulomb return mapping works in sorted principal stresses, so the kernel needs a
symmetric 3x3 eigen-decomposition without OpenFOAM or an external library (decision
0006). Repeated eigenvalues are the normal case at the edges (triaxial states) and the
apex (hydrostatic), exactly where the model must be exact.

## Decision

- Eigen-solver: cyclic Jacobi rotations (`kernel/celEigen.{H,C}`), not the closed-form
  trigonometric solution. Jacobi is accurate to round-off relative to the matrix norm
  for repeated and nearly repeated eigenvalues and always returns orthonormal vectors;
  the closed form loses digits there. Cost (a few sweeps of 3 rotations) is small next
  to the rest of a cell update. Unit-tested by `tests/constitutive/eigenCheck.C`.
- The algorithm copies `mohr_coulomb` in `models.py` step for step (main plane, then the
  compression or extension edge, then the apex; the same tolerances), with 1x1 and 2x2
  solves written out. With repeated eigenvalues the principal directions differ from
  numpy's, but the returned stress is unique, and the reference test compares stress.
- `MohrCoulomb::updateReturn` also reports which return was taken (elastic, plane,
  compression edge, extension edge, apex). The reference test uses it to require that
  each path exercises the return it is there for, so a path that silently stops reaching
  an edge or the apex fails.
- Kernel name `mohrCoulomb`; dictionary parameters `E nu phi c psi` (degrees, Pa), all
  required. Constructor rejects c < 0, phi outside [0, 90), psi outside [0, phi], and
  phi = 0 with c = 0.
- `plasticCorrector` (history from the corrector, `_finish` in `models.py`) is a shared
  kernel function, with the same signature and body as WP 2.5a's, so the two packages
  merge without a semantic conflict.
- No tension cut-off. The Python reference has none (decision 0004); the brief says not to
  add one before the reference does.

## Consequences

- Cohesive material can carry tension up to the apex c cot(phi) (17.3 kPa for
  c = 10 kPa, phi = 30). For sand (c = 0) the apex is at zero stress, so no cut-off is
  needed. A Rankine cut-off is still wanted for cohesive soils; it needs a Python
  reference first (new small package, see `docs/work-packages.md`).
- Tresca (phi = 0) has no apex; a failed corner return throws, as in Python.
