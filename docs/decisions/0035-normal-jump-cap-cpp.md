# 0035 Rule 12 (normal jump cap) in the C++ Lagrange step (WP 3.4b)

Date: 2026-11-16. Status: accepted as an option (`normalCap yes`, default off). WP 3.4b
closes with its done-when not met; the choice between a soft numerical E and mu(I) is
put to Sia.

## Context

Decision 0034: with the collapse code (shear cap, rule 9) the C++ discharge rate falls
with E (10 mm: 74.1 / 65.7 / 70.9 kg/s/m at 0.1 / 1 / 10 MPa; 5 mm: 78.0 / 55.0 at
0.1 / 1 MPa, -29 %). The cause named in decision 0028 is the normal jump term of the
acoustic Riemann solver, a numerical viscosity of order rho c h on plastic flow. In the
split prototype capping it at yield (rule 12) broke the column collapse (decision 0029).

## Decision

Rule 12 is added to `eulerianSolidLagrange.C` (`capNormal`) behind the materialProperties
switch `normalCap yes` (uses `phi` and `c`; default off, so every existing result is
unchanged bit for bit). At each contact face, after the shear cap, the normal jump term
zpL zpR (vnR - vnL)/(zpL + zpR) is scaled back by bisection (20 steps) until the face
stress (tangential block of the mean of the two cells' stresses, n row from the
traction) is no further outside Mohr-Coulomb than without the jump. One traction per
face (momentum exact) and one normal velocity (geometry and remap, mass exact). A full
cell strains with its own characteristic velocity at capped faces, deviatoric part
only: the extra strain symm(dt sum Sf dn / V) has its trace over the geometric
directions removed, so the volume stays tied to the density (decision 0028). Rule 13 is
not added: the unsplit step has no one-sweep void.

The test drivers take extra materialProperties lines from the environment variable
`CEL_EXTRA_MATERIAL` (`tests/void/void_checks.py`), so the collapse and discharge checks
run with the switch on unchanged:
`CEL_EXTRA_MATERIAL='normalCap yes;' tests/manual/test_discharge_stiffness.sh`.

## Result (rule 12 on; shear cap and rule 9 on as in decision 0034)

Discharge, T0 tank, D = 0.06 m, rate over 97 to 88 % of the mass at opening:

| Run | 0.1 MPa | 1 MPa | 10 MPa | Without rule 12 (0034) |
| --- | --- | --- | --- | --- |
| 10 mm, serial | 78.6 (C 0.95) | 73.8 (-6.1 %) | 84.8 (+7.9 %) | 74.1 / 65.7 / 70.9 |
| 5 mm, 3 processors | 82.6 (C 1.00) | 69.7 (-15.7 %) | not run | 78.0 / 55.0 |

Column collapse (`tests/collapse/collapse_checks.py`, 10 mm):

| Case | Run-out with rule 12 | Without (0032) | Prototype |
| --- | --- | --- | --- |
| a = 0.5, 0.1 MPa | 0.50, slope 19.3 deg | 0.55 | 0.55 |
| a = 0.5, 1 MPa (checked) | 0.35, slope 22.8, top kept: PASS | 0.35 | 0.40 |
| a = 0.5, 10 MPa | 0.25, slope 25.5 | 0.20 | 0.25 |
| a = 2, 1 MPa (checked) | 1.50, stops at 0.58 s: PASS | 1.40 | 1.8 |

- Unlike the split prototype (decision 0029), rule 12 does not break the C++ collapse:
  no run-away at any E, run-outs within 0.05 of the old ones, the top kept. The unsplit
  step does not need the normal jump for stability.
- At 10 mm the done-when is met (spread -6.1 / +7.9 % of the 0.1 MPa rate). At 5 mm it
  is not: the fall from 0.1 to 1 MPa halves (-29 % to -16 %) but stays above 10 %. So
  the normal jump was a large part of the E dependence, not all of it.
- Mass balance 2e-14 to 5e-14 and zero stress in empty cells in every run.
- Parallel: the discharge on two processors gives 79.11 against 78.60 kg/s/m serial
  (+0.64 %, within the 2 % of `tests/test_discharge.sh`). The collapse 2-processor check
  (1e-10 at 0.1 s) fails with the switch on: 1e-11 at 0.002 s, 4e-8 at 0.01 s
  (without the switch 3e-16 at 0.01 s). The bisection returns the cap factor on a 2^-20
  grid, so round-off differences between serial and parallel inputs move it by a whole
  grid step at faces near the threshold; the plastic flow then amplifies it (as for the
  shear cap later in the run, decision 0032). Not a sync error: the sync of the normal
  slips reuses the shear-slip path, and the discharge rates agree.

## Consequences

- `normalCap yes` is usable for T0 at 10 mm, but the 5 mm result says the rate would
  still drift with E on finer meshes. The switch stays off by default and in CI.
- The remaining E dependence at 5 mm is not the normal jump. Candidates not tried here
  (to keep the package small): the shear cap leaves the shear jump of still-elastic
  faces, and the outlet zone may still hold void. mu(I) adds rate dependence that is
  physical rather than numerical.
- Choice for Sia (decision 0034 item 4): run T0 at a soft numerical E = 0.1 MPa (rate
  C 0.95 to 1.0, collapse as the prototype; explicit codes use mass scaling the same
  way) or build mu(I) rate-dependent friction.
- If the switch is used in parallel, compare runs by rate, not cell by cell.
