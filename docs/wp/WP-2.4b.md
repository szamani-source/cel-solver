# WP 2.4b: Lagrange step plus remap, integration checks

State: done (2026-11-14)
Branch: `claude/project-thread-504odx`

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

### 2026-11-14: done in one session

Built `tests/eulerian/cycleCheck` (`celCycleCheck`: initial U and sigma from `0/`,
fixed dt, Lagrange step then remap per cycle, records totals and wave-frame ranges
every step, optional prescribed rotation with U held), `tests/eulerian/cycle_checks.py`
and `tests/test_eulerian_cycle.sh` (about 3 minutes). Every "done when" item passes;
numbers and the stability study are in decision 0021.

The checks exposed two bugs, both fixed:

- Prescribed motion (`solveMomentum no`) used the acoustic face velocity, whose
  stress-jump term blew up the rotating disc. Now the mean of the face states.
- Remap limiter divided by -0 for a subnormal gradient (1D plastic wave on a 3D mesh,
  floating-point exception). Now skips extrapolations below 1e-150.

Learned: the WP 2.3 "second crossing not above the first" stability rule is noise
(about 0.5 percent from crossing to crossing); the full cycle is stable to maxCo 0.8 in
2D and at least 0.7 in 3D. Known issues recorded (S step at 45 degrees overshoots
0.9 percent with streaming, 1.4 without; slow bounded creep of a 45 degree step's
overshoot over many crossings).

Next: WP 3.1 (C++ gravity and walls), after WP 3.1p.
