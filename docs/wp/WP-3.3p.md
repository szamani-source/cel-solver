# WP 3.3p: Prototype silo discharge (target T0)

State: blocked on WP 3.2p
Branch: `wp-3.3p-silo-discharge`

## Goal

The complete target problem T0 in the 2D prototype: fill, settle with the outlet closed,
open the outlet, discharge. This is the first end-to-end result to show Sia.

## Read first

- `docs/spec/target-problems.md` (T0), `docs/spec/silo-reference.md`
- `prototype/silo2d.py`; logs of `docs/wp/WP-3.1p.md` and `docs/wp/WP-3.2p.md`;
  the void decision note

## Build

- Outlet boundary: a segment of the bottom wall that is a rigid wall until `t_open`,
  then an outflow boundary (zero traction; material crossing it leaves the domain and is
  counted). Opening instantly, and over a short ramp as an option if the instant
  opening gives a start-up spike.
- Run script for T0 with the spec's geometry and material; outputs: mass in the tank and
  discharged mass against time, discharge rate, free-surface profile, velocity and
  plastic-strain fields at set times as PNG images, and one animation (GIF or MP4) of
  the discharge for Sia.
- Cost: measure the run time; if a 3 s discharge takes more than about 30 minutes, test
  a reduced stiffness E and show the discharge rate is unchanged within 2 percent, then
  use it (record in the spec's cost note).
- `tests/test_proto_silo_discharge.sh` (a coarse, short version for CI, under 10 minutes),
  and `docs/results/silo-prototype.md` with the figures and the comparison table.

## Done when

`tests/test_proto_silo_discharge.sh` passes (coarse mesh) and the full run is reported:

- Nothing leaves before `t_open`; mass in tank plus discharged mass constant to 1e-8.
- After a start-up transient the discharge rate is steady: variation below 10 percent
  while the bed is deeper than one tank width.
- Rate independent of head: runs with H = 0.6 m and H = 0.9 m give the same steady rate
  within 10 percent.
- Rate against outlet width D = 0.04, 0.06, 0.08 m: fitted exponent within 1.5 +/- 0.2,
  and the constant within the range in `docs/spec/silo-reference.md`.
- Flow pattern: stagnant zones in the bottom corners (velocity below 1 percent of the
  outlet velocity), a central channel, a crater in the free surface; figure in the
  results page.
- Mesh check: 5 mm and 2.5 mm cells give steady rates within 10 percent.
- Existing tests still pass unchanged.

Where a criterion fails for a physical reason that the literature confirms (for example
perfect plasticity over-predicting the rate), do not tune it away: record the numbers,
the source, and what constitutive change would address it.

## Out of scope

- Cohesion (next package after this), hopper bottoms, 3D, C++.

## Log

(none yet)
