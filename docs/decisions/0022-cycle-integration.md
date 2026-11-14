# 0022 Lagrange step plus remap: integration results and two fixes (WP 2.4b)

Date: 2026-11-14. Status: accepted (Sia's standing rule: proceed and merge on green CI).

Numbers 0019 to 0021 are taken by WP 3.0, WP 3.1p and WP 3.2p (void). This note updates decision 0013
item 2 ("unsplit in C++") with what the C++ cycle does, and closes the risk recorded
there ("the unsplit C++ scheme is not the scheme prototyped").

## Context

WP 2.3 (Lagrange step, decision 0018) and WP 2.4 (remap, decision 0017) were built and
tested separately. WP 2.4b runs them together, as `celFoam` does, on the 2D prototype
checks of WP 1.2 and 1.3b: `tests/test_eulerian_cycle.sh` with the driver
`tests/eulerian/cycleCheck` (`celCycleCheck`) and `tests/eulerian/cycle_checks.py`.
Material: rho 2000 kg/m^3, E 20 MPa, nu 0.3 (von Mises Y 10 kPa; Drucker-Prager
phi 30, c 5 kPa), default maxCo 0.4, material streaming at 0.25 c_p.

## Decision

### 1. Prescribed motion: face velocity without the stress-jump term

With `solveMomentum no` the Lagrange step used the acoustic Riemann face velocity,
which adds (stress jump)/(sum of impedances) to the mean of the two sides. With U held
fixed nothing balances that term. For the pre-stressed disc in rigid rotation (stress
jump 100 kPa at the edge, face velocity correction about 0.2 m/s against a rim speed of
0.5 m/s, at an acoustic Courant number near 36 because the time step follows the
material speed) the run blew up within the first turn (remap: negative density).
Prescribed motion now takes the mean of the two reconstructed face states (internal and
coupled faces), as the prototype's kinematic mode takes the prescribed face velocity.
For a linear velocity the reconstruction is exact, so the face velocity is exact. With
`solveMomentum yes` nothing changes. The WP 2.3 rotation and shear checks (uniform
stress) still pass unchanged.

### 2. Remap limiter: skip vanishing gradients

The 1D plastic step wave on a 1000 x 2 x 2 mesh stopped with a floating-point
exception in the remap after about 35 steps: round-off in the transverse velocity gave
a subnormal gradient, the relative test `|g.d| <= 1e-10 |g||d|` underflowed to `0 <= 0`
and let it through, and half of the smallest subnormal rounded to -0 before dividing.
A face is now skipped also when the extrapolated difference is below `ROOTVSMALL`
(1e-150), which cannot change any result. Lesson for later limiters: every ratio needs
an absolute floor as well as a relative one.

### 3. Results (all pass; the test takes about 3 minutes)

| Check | C++ cycle | Prototype | Criterion |
| --- | --- | --- | --- |
| P / S along x, streaming: speed (25, 50, 100 cells) | +0.075 % / +0.18 % worst | | 1 % |
| P / S along x: L1 order | 2.00, 2.08 / 1.80, 1.94 | 1.77 to 1.87 | 1.7 |
| P / S at 45 deg: speed | +0.026 % / +0.036 % worst | | 1 % |
| P / S at 45 deg: L1 order | 2.11, 2.11 / 2.04, 2.05 | | 1.7 |
| Step waves after 10 cells: P x, S x, P 45, S 45 | 0.000, 0.003, 0.24, 0.92 % | | 1 % |
| Mass, momentum (periodic, streaming) | 1.2e-13 | | 1e-12 |
| Disc, one turn, principal stresses inside | 0.11 % | 0.03 % | 2 % |
| Disc, beyond the 6-cell edge | 0.06 % | | 1 % |
| Plastic step, 3D mesh: precursor / plastic speed / HEL | +0.04 / +0.09 / -0.11 % | WP 1.2 | 2 % |
| 45 deg impact, extra over along x after start-up: VM / DP / elastic | 0.000 / 0.022 / 0.000 % | 0.000 / 0.082 / 0.000 % | 1 % |
| 45 deg impact start-up (reported) | 0.13 / 0.36 / 0.21 % | 2.7 % (elastic step) | |
| 2 processors vs serial (45 deg VM impact) | about 1e-13 | | 1e-12 |

The scheme prototyped (split) and the scheme built (unsplit) agree on every check; the
C++ cycle is at least as good except for the disc (0.11 against 0.03 percent, still far
inside 2 percent). Decision 0013 item 2 stands; nothing in items 1, 2 or 6 needs revisiting.

### 4. Stable Courant number of the full cycle (`cycle_checks.py WORKDIR --study`)

Step P-wave at 45 degrees (2D, 50 x 50) or along the body diagonal (3D, 20^3), streaming
at 0.25 c_p, 16 crossings; unstable means the velocity passes 1.5 times the jump.

| maxCo | 2D: start-up overshoot, verdict | 3D: start-up overshoot, verdict |
| --- | --- | --- |
| 0.4 (default) | 0.5 %, stable | 2.2 %, stable |
| 0.5 | 0.9 %, stable | 2.8 %, stable |
| 0.6 | 1.3 %, stable | 7.7 %, stable |
| 0.7 | 1.4 %, stable | 19 %, stable |
| 0.8 | 11 %, stable | not run |
| 0.9 | blows up | not run |

So the remap does not lower the Lagrange step's limit. WP 2.3's numbers (2D 0.65,
3D 0.50) used a strict "second crossing not above the first" rule, which this package
found to be noise: the largest value moves by about 0.5 percent from crossing to
crossing as the front shifts against the mesh, at any Courant number. The default
stays 0.4: above about 0.6 in 3D the start-up overshoot, not stability, is the limit.
The test checks 4 crossings at the default (largest below 1.05).

### 5. Cost

2D 100 x 100, linear elastic, serial, second order, cloud container: Lagrange step
8e-7 s, remap 1.4e-6 s, whole cycle 2.3e-6 s per cell per step. The remap is the
larger part (three least-squares gradients and limiters per field component); worth
profiling before the silo runs (WP 3.4) if run time matters.

## Consequences

- Known issue, not fixed here: a 45 degree S step overshoots by about 0.9 percent after
  start-up in streaming material and 1.3 to 1.4 percent without streaming (so it is the
  Lagrange step's limiting of primitive components, decision 0018 item 1, not the remap).
  The P step at 45 degrees is 0.24 percent. Characteristic limiting, as in the 1D
  prototype, would be the fix if contact (WP 4.3) or the silo needs it.
- Known issue: over many crossings the overshoot of a 45 degree P step creeps up and
  levels off (2D, maxCo 0.4: 0.5 percent at first, 2.7 percent after 16 crossings,
  5.7 percent after 35 to 40), with or without streaming. Bounded, so not an
  instability; it matters only for long runs of sharp elastic fronts.
- The C++ phase 2 risk list is closed; WP 3.1 to 3.4 (silo in C++) can build on the cycle.
