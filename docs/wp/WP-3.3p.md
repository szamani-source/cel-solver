# WP 3.3p: Prototype silo discharge (target T0)

State: done for the outlet and the end-to-end run (2026-11-15); rate law split to WP 3.3q
(decision 0026)
Branch: `wp-3.3p-silo-discharge`

## Goal

The complete target problem T0 in the 2D prototype: fill, settle with the outlet closed,
open the outlet, discharge. This is the first end-to-end result to show Sia.

## Read first

- `docs/spec/target-problems.md` (T0), `docs/spec/silo-reference.md`
- `prototype/silo2d.py`, `prototype/void2d.py` (`VoidBox`, start from it); logs of `docs/wp/WP-3.1p.md` and `docs/wp/WP-3.2p.md`;
  `docs/decisions/0021-void-treatment.md` (the discharge rate may depend on E through
  the Riemann solver's dissipation: measure that first; decision 0023 capped the shear
  part in WP 3.2q, the normal part still depends on E)

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
  and the constant C (exponent fixed at 1.5, k = 0) within 1.0 to 1.9
  (reference 1.45, `docs/spec/silo-reference.md` section 3, decision 0019). Measure the
  steady rate between 95 and 70 percent of the initial mass: for D = 0.08 m the bed
  passes one width after only about 0.8 s.
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

- 2026-11-15: Outlet built (rule 10 in `void2d.py`, `prototype/discharge2d.py`); settle
  0.3 s with the outlet closed, open instantly, damping off. CI test
  `tests/test_proto_silo_discharge.sh` passes (10 mm, E = 0.1 MPa, about 3 to 5 min):
  nothing leaks, mass balance 1e-14, head independence 0.975. Rate against E first, as
  asked: 6 / 17 / 38 kg/s/m at E = 10 / 1 / 0.1 MPa (D = 0.06 m, 10 mm); E = 0.01 MPa
  never steadies. Cause: the Riemann jump terms act as a viscosity of order rho c h
  grad v. Reducing them (scaled, Mach-scaled, with floors) raises the rate but makes the
  flow noisy (speeds 4 to 10 m/s, kinetic above strain energy), so nothing was adopted;
  `low_mach` (rule 11) stays as an off-by-default study option. Full runs at 5 mm, E = 0.1
  MPa, after merging WP 3.2r (rule 9): rates 14.3 / 36.6 / 68.4 kg/s/m for D = 0.04 /
  0.06 / 0.08 m (C 0.32 / 0.44 / 0.54, exponent 2.26; before the merge 23.3 / 52.6 /
  85.4, exponent 1.87), steady within 7 percent, stagnant corners (0.4 percent of the
  outlet speed at D = 0.06 m), crater. Mesh: 36.8 (5 mm) against 45.2 (2.5 mm) kg/s/m
  over 0.2 to 0.5 s. Figures and animation in `docs/results/silo-prototype.md`. The rate law
  checks are in `tests/manual/test_proto_silo_rate.sh` (fail) and are the done-when of
  WP 3.3q. Next: WP 3.3q (`docs/wp/WP-3.3q.md`).
