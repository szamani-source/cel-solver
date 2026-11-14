# WP 1.3b: 2D rotation, plasticity and the formulation decision

State: blocked on WP 1.3a
Branch: `wp-1.3b-2d-rotation-decision`

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

(none yet)
