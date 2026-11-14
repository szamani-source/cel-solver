# WP 3.1p: Prototype gravity, walls and settling (Janssen)

State: done (2026-11-14)
Branch: `wp-3.1p-gravity-walls`

## Goal

The 2D Python prototype fills a closed tank with cohesionless Mohr-Coulomb material,
applies gravity, and settles to rest with the wall stresses Janssen predicts. This is
stage 1 of target T0 (`docs/spec/target-problems.md`).

## Read first

- `docs/spec/target-problems.md` (T0), `docs/decisions/0013`, `0016`
- `prototype/lr2d_plastic.py`, `prototype/lr2d.py`, `prototype/constitutive/batch.py`
- `docs/wp/WP-1.3b.md` log. If `docs/spec/silo-reference.md` exists (WP 3.0), use its
  Janssen values and tolerances; otherwise derive Janssen in the test and say so.

## Build

In a new `prototype/silo2d.py` that reuses `lr2d_plastic.py`:

- Gravity as a body force in the Lagrange step, consistent with the Godunov scheme
  (well-balanced: a column at rest in equilibrium must stay at rest to round-off, not
  drift; say how this is achieved, for example by reconstructing the stress minus the
  equilibrium gradient, or a source term in the Riemann states).
- Wall boundaries: rigid smooth wall (zero normal velocity, zero shear traction) and
  rigid frictional wall (zero normal velocity, shear traction limited to
  `tan(delta_w)` times the normal pressure, opposing sliding; sticking below the limit).
  No tension at a wall: material may separate (normal traction not positive).
- Mohr-Coulomb in the vectorised batch models if it is not there yet (c = 0, psi = 0,
  tension cut-off 0), checked against `models.py` on random increments.
- Settling: gravity ramped up over a set time, with a damping option (for example
  velocity-proportional damping or dynamic relaxation) that is switched off for the
  discharge; a "settled" measure (kinetic energy below a threshold).
- The tank is completely full in this package (no free surface: the top is a free
  boundary with zero traction, material height = domain height). Void is WP 3.2p.
- `tests/test_proto_silo_settle.sh`

## Done when

`tests/test_proto_silo_settle.sh` passes:

- Lake at rest: a column started in exact K0 equilibrium under gravity between smooth
  walls keeps velocities below 1e-10 m/s for 1000 steps.
- Smooth walls, settled from zero stress: vertical stress at the bottom within 2 percent
  of `rho g H`; horizontal stress consistent with the material: ratio nu/(1 - nu) =
  0.4286 within 5 percent (elastic, no yield; `docs/spec/silo-reference.md` section 2).
- Frictional walls (20 degrees), T0 geometry with a tall fill (H = 3 W) so saturation
  shows: mean vertical stress against depth follows Janssen within 10 percent below one
  width from the top, and the wall shear to normal stress ratio is `tan(delta_w)` within
  5 percent where the wall slips. Record the K that fits; it must lie in 0.27 to 0.45
  (between the wall-yield value 0.29 and the elastic value 0.43, decision 0019).
- Total vertical force balance: weight = bottom force + wall friction, within 1 percent.
- Settled state reached (kinetic energy criterion) and unchanged over a further 10 percent
  of the run with damping off.
- Existing tests still pass unchanged. The new test runs in under 5 minutes.

Record in the log: steps and wall-clock time to settle, the damping used, and the time
step (cost note in the spec).

## Out of scope

- Free surface and void (3.2p), the outlet (3.3p), cohesion, C++.

## Log

### 2026-11-14 (thread "Prototype gravity, walls and settling")

Done. `prototype/silo2d.py` (class `SiloBox`, a `PlasticBox2D` with its own sweep
`silo_sweep`), `tests/test_proto_silo_settle.sh` passes in about 170 s here. Decision 0020.
Branch was `claude/project-thread-ql1dhs` (thread branch), not `wp-3.1p-gravity-walls`.

- Mohr-Coulomb vectorised in `constitutive/batch.py`, `check_batch_mc()`: three cases
  (T0 c = 0 with sigma_t = 0; c = 5 kPa psi = 10; c = 5 kPa with sigma_t = 2 kPa) match
  `models.py` to 3e-14 with edge and apex returns exercised. For c = 0 the tension
  cut-off at 0 is redundant (the cone already implies s1 <= 0). `check_batch()` is
  unchanged, so the WP 1.3b test is unchanged.
- Well-balanced gravity: Lagrange step on s_yy - G (discrete weight profile per column).
  Lake at rest (0.1 x 0.6 m, 10 mm, K0 state, 1000 steps): 2.5e-18 m/s.
- Smooth walls, 0.1 x 1.2 m, 10 x 120 cells, from zero stress: settled in 2780 steps
  (dt 57.8 us, 0.161 s, 35 s wall clock); bottom s_yy 0.126 % from rho g H; s_xx/s_yy =
  0.4286 = nu/(1 - nu) everywhere below the top 10 percent (elastic, as expected: the
  active Rankine limit is 0.271).
- Frictional walls and bottom (20 deg), 0.4 x 1.2 m, 40 x 120 cells: settled in 2520
  steps (0.146 s, 109 s wall clock). K measured at the walls 0.425, best fit 0.411; mean
  s_v within 1.67 % of Janssen (measured K) below z = W (bottom 14.0 kPa against 21.1 kPa
  hydrostatic). Wall shear/normal (cells extrapolated to the wall) 0.985 to 1.018 of
  tan 20 deg on the slipping faces below z = W (the wall slips on 110 of 120 faces; the
  bottom corner sticks). Weight = bottom + wall friction to 0.08 %. Mass exact.
- Damping used: gravity ramp 0.1 s (1 - cos), damping alpha 250 1/s, settled at
  KE <= 1e-6 strain energy. Without damping for a further 10 % of the steps: stress
  change 1.5 % (frictional) and 0.24 % (smooth) of rho g H, KE stays about 1e-6.
- Tried: ramp 0.05 s with alpha 200 settles faster (0.08 s) but the column rebounds and
  the upper wall unloads to 0.76 to 0.98 of the friction limit, failing the 5 % check;
  alpha 500 with ramp 0.08 s passes but takes 3660 steps. 20 x 60 cells: Janssen 7 %
  with the measured K, wall ratio down to 0.53 near the bottom: too coarse.
- Cost note: at 5 mm cells the frictional tank would be about 8 times the work
  (4 times the cells, twice the steps), about 15 minutes in this Python.

Next: WP 3.2p (free surface) starts from `SiloBox`; replace the no-flux free top.

