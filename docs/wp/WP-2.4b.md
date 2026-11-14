# WP 2.4b: Lagrange step plus remap, integration checks

State: blocked on WP 2.3 and WP 2.4
Branch: `wp-2.4b-integration-checks`

## Goal

The full C++ cycle (Lagrange step, remap) reproduces the 2D prototype checks, closing
the risk recorded in decision 0013 ("the unsplit C++ scheme is not the scheme
prototyped").

## Read first

- `docs/decisions/0013`, `0014`; logs of `docs/wp/WP-2.3.md` and `docs/wp/WP-2.4.md`
- `tests/test_proto_2d_rotation.sh`, `prototype/lr2d_plastic.py`

## Build

- Tutorials and `tests/test_eulerian_cycle.sh`. No new solver features; fix what the
  checks expose (each fix with a test), or record it as a known issue with numbers.

## Done when

`tests/test_eulerian_cycle.sh` passes in CI:

- Elastic P and S waves in material streaming at 0.25 c, along x and diagonally:
  speeds within 1 percent, order 1.7 or better, no new extrema above 1 percent after
  10 cells of travel, conservation to 1e-12.
- Pre-stressed disc in rigid rotation, one turn, 100 cells across: principal stresses
  within 2 percent (prototype: 0.03 percent), nothing beyond the smeared edge.
- 1D plastic step wave (von Mises) on the 3D mesh: precursor and plastic wave speeds and
  the HEL within 2 percent of the analytic values (WP 1.2 criteria).
- 45 degree elastic-plastic impact: overshoot after start-up within 1 percent.
- A decision note updating decision 0013 item 2 with the measured stable Courant
  numbers and the cost per cell per step.

## Out of scope

- Void, more than one material, contact, benchmarks against experiments (2.6).

## Log

(none yet)
