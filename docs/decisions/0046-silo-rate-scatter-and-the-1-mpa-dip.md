# 0046 Silo rate scatter and the 1 MPa dip (WP 3.6f)

Date: 2026-11-18. Status: accepted (finding recorded; no solver change).

## Context

WP 3.6d (decision 0045) left the T0 discharge with mu(I) and `flowImpedance yes` about
20 % slower at 1 MPa than at 0.1 and 10 MPa, from single runs, and showed that a single
10 mm rate can move about 15 % with round-off. WP 3.6f measures each point as a band and
looks for the cause.

## Method

`tests/manual/test_discharge_scatter.sh` (`tests/discharge/scatter.py`): each E three
times with maxCo 0.40, 0.38, 0.36 (same physics, a different round-off history), 10 mm.
Suspects switched with temporary environment switches (not committed).

## Results (kg/s per m depth; three runs, mean and range)

| Variant | 0.1 MPa | 1 MPa | 10 MPa |
| --- | --- | --- | --- |
| As merged (decision 0045 snap) | 65.8 (65.5 to 66.3) | 52.0 (50.6 to 53.8) | 60.8 (56.5 to 67.7) |
| No nearly-full snap at all | 65.5 (64.6 to 66.3) | 60.5 (58.5 to 61.6) | 63.1 (62.4 to 63.8) |
| Snap tolerance 1e-2 instead of 1e-4 | 65.4 (63.9 to 66.2) | 49.8 (48.9 to 50.9) | |
| Snap tolerance 1e-7 | | 58.2 (53.9 to 63.3) | |
| No snap; void < 1e-4 counted as full in the Lagrange rules only (no volume change) | 63.2 (62.7 to 63.8) | 43.7 (42.7 to 44.6) | |
| Same with void < 1e-2 counted as full | 83.9 (83.6 to 84.3) | 67.0 (64.8 to 68.6) | |

Inclined layer (`tests/manual/test_inclined_layer.sh`): passes as merged; fails with
tolerance 1e-7 (1 MPa free surface +20 %, 0.1 MPa lid +78 %); without the snap it fails
as in decision 0044; passes with the rules-only variant (errors 3.6 / 0.9 / 0.9 %).

## Findings

1. The 1 MPa dip is real: the bands of 0.1 and 1 MPa are far apart (65 to 66 against 51
   to 54). Scatter is 1 to 6 % at 0.1 and 1 MPa and 18 % at 10 MPa with the snap.
2. The nearly-full rule (decision 0044, bounded in 0045) causes both the dip and the
   10 MPa scatter: without it the three E agree within 8 % (65.5 / 60.5 / 63.1) and
   every band is within 5 %, which would pass the 10 % check at 10 mm. But the inclined
   layer needs some such rule, and no tolerance tried serves both.
3. The rate is set by whether the slightly loosened cells above the outlet count as
   full (flow viscosity, second order, stress) or partly filled (void-first, no stress
   in the opening direction): moving that line moves the rate by up to 30 % at fixed E
   (43.7 to 67.0 at 1 MPa). The friction law is not the lever; the void closure is. At
   0.1 MPa the elastic strains are large, the void per step exceeds every tolerance, and
   the cells go partly filled at once; at 1 MPa they stay full; at 10 MPa they hover.
4. So the E independence of the no-snap runs rests on the same softening of loosened
   cells that decision 0029 found in the prototype, not on the physics. A discharge
   rate that does not depend on E needs the loosening above the outlet to be a
   material property: a dilatant, compressible mu(I)/phi(I) law (solid fraction as a
   function of I), already the preferred cure for the shear bands (decision 0038).

## Decision

- No solver change in this package: the decision 0045 rule stays (it keeps the layer
  right and fixed the 10 MPa collapse of the rate). The 1 MPa dip is a known issue.
- Silo rates are compared as bands of three runs from now on
  (`tests/manual/test_discharge_scatter.sh`).
- Proposed next package (for Sia to approve, since it changes the model):
  compressible mu(I)/phi(I), so loosened material is a state of the material, with
  the silo E sweep and the collapse shear-band test as its checks.

## Consequences

- The tank case (`tutorials/siloDischarge`) is usable at any E, with a rate uncertainty
  of about 20 % from the 1 MPa dip and the scatter.
- `outlet_checks.MAXCO` sets the step fraction in the outlet and discharge drivers.
