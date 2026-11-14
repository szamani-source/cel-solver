# WP 3.1p: Prototype gravity, walls and settling (Janssen)

State: not started
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

(none yet)
