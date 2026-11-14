# WP 1.3b: 2D rotation, plasticity and the formulation decision

State: done (2026-11-14)
Branch: `claude/project-thread-f3dbw1` (thread branches replace the `wp-` names)

## Goal

The 2D prototype handles rotating, yielding material, and a decision note fixes the
formulation the C++ solver will implement (WP 2.3, 2.4).

## Read first

- `docs/wp/WP-1.3a.md` log, `prototype/lr2d.py`
- `docs/decisions/0005`, `0008`; known issues in `docs/STATUS.md`

## Build

- Objective stress update (Jaumann rate, or an incrementally objective rotation such as
  Hughes-Winget) using the spin from the face velocities.
- Stress update through the WP 1.4 models (von Mises, Drucker-Prager), with history
  remapped as in WP 1.2 (superbee for history, MC for the rest).
- `tests/test_proto_2d_rotation.sh`
- `docs/decisions/NNNN-eulerian-formulation.md` answering, with numbers from the tests:
  1. Cell-centred Godunov (as prototyped) or node-based staggered (as Abaqus) for 3D?
  2. Direction-split or unsplit?
  3. Which objective rate?
  4. Impulsive start-up overshoot (decision 0008): accept, elastic-plastic Riemann
     solver, or artificial viscosity at contact?
  5. Per-field limiter choice for the C++ remap.
  6. Oblique fronts: the WP 1.3a scheme overshoots 3 to 14 percent on a step P-wave at
     45 degrees (decision 0012). Accept, add a multidimensional (nodal) Riemann solver,
     or go staggered?

## Done when

`tests/test_proto_2d_rotation.sh` passes:

- A pre-stressed elastic disc in rigid rotation (velocity field prescribed, one full
  turn): stress tensor rotates with the material; error in the principal stresses below
  2 percent at 100 cells across the disc; no stress appears outside the disc beyond the
  smeared edge.
- Simple shear of a von Mises block to shear strain 1: shear stress stays at the yield
  value within 1 percent, no oscillation (the classic objective-rate check).
- 2D plastic step wave along x reproduces the WP 1.2 result on the same resolution.
- The decision note exists and STATUS.md points to it.

## Out of scope

- More than one material, void, contact, 3D.

## Log

### 2026-11-14

Done:
- Wrote the checks first (`lr2d_plastic.py --verify`, run by
  `tests/test_proto_2d_rotation.sh`), then `prototype/lr2d_plastic.py`: the WP 1.3a split
  scheme with the full 6-component stress and WP 1.4 history per cell, Hughes-Winget
  rotation (Jaumann) from each sweep's half of the spin, mass-weighted remap of velocity,
  stress and history (MC, history superbee), and a prescribed-velocity mode for the
  rotation and shear checks. `prototype/constitutive/batch.py` vectorises the
  linear elastic, von Mises and Drucker-Prager models (match `models.py` to 4e-15).
  All checks pass; the test takes about 1.5 minutes.
- Disc in rigid rotation, 100 cells across, one turn: principal stresses within
  0.027 %, direction within 0.0006 deg, stress outside below 0.015 % (6-cell band each
  side of the edge; the edge is 10 cells wide, 10-90 %). Interior is exact to 1e-12.
- Simple shear to gamma 1: von Mises shear stress within 3e-7 of Y/sqrt(3), no drop;
  elastic matches the Jaumann closed form G sin(gamma) within 2.1e-4 G; Drucker-Prager
  at p = -50 kPa stays on the cone.
- Plastic step along x on a 1000 x 2 mesh matches `lr1d_plastic.py` to 1e-11.
- Impact at 45 degrees (von Mises, Drucker-Prager, elastic) adds at most 0.08 % over
  the same impact along x after start-up.
- Decision note `docs/decisions/0013-eulerian-formulation.md` answers the six questions.

Learned:
- The 13 % overshoot at 45 degrees in WP 1.3a (decision 0012) was mostly a jagged initial
  step: `x + y` from cell centres in floating point puts cells of one diagonal on both
  sides of the threshold. With integer indices: 2.7 % in the first cell of travel, 0.22 %
  after 10. `lr2d.check_step` fixed and the 45 degree step is now checked (after a
  10-cell start-up) in `test_proto_2d_elastic.sh`.
- In a split scheme the constitutive model must be skipped where a sweep's strain
  increment is exactly zero; otherwise a remapped stress slightly outside the yield
  surface is returned early and a 1D problem on the 2D mesh differs from 1D by 5e-4.
- An unsplit proxy (sweep increments added) is clean at 45 degrees but needs CFL about
  0.3 in 2D; it is not usable with plasticity (returns per sweep, then adds).

Next: WP 2.3 (C++ stress update: unsplit face Riemann solver, Hughes-Winget, one model
call per cell per step) and WP 2.4 (remap with per-field limiters), per decision 0013.
Both briefs still to write.
