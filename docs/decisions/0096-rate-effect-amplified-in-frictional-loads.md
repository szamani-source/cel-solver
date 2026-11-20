# 0096: The sand rate effect is amplified friction, not one scheme term (A.5 findings)

Date: 2026-11-19 (accuracy package A.5, thread "Fix sand footing bearing capacity")
Status: finding recorded; scheme change parked (see Decision).

## Context

Decision 0095 found the sand footing reading scales with rho c v (push speed times
impedance): psi 30 N_gamma 20 at 0.01 m/s, 11.9 at 0.00125 m/s, exact 7.653. A.5 looked
for the term that carries it, to make slowly yielding Mohr-Coulomb rate-free.

## Findings (psi 30, 10 mm cells, 0.01 m/s, implied N_gamma at s 0.02 to 0.05 B; base 20)

| Change (experiment only, not kept) | N_gamma |
| --- | --- |
| Face traction jump terms (normal and shear) times 0.25 | 17.7 to 18.2 |
| Normal jump only times 0.25 / shear only times 0.25 | 18.6 to 19.2 / 19.0 to 19.3 |
| Hancock predictor from the material model (plastic) instead of elastic | 16.5 to 17.1 |
| Body wall face impedance times 0.25 | 19.5 to 19.9 |
| Face velocity stress-jump term times 0.25 | 20.6 to 20.9 (higher) |
| All face terms times 0.25 | unstable |
| Velocity damping 50 1/s (a drag of about rho 50 v, 0.9 kPa/m) | 22.9 to 23.7 |
| Every length times 10 (B 1 m, same cells; rho c v / gamma B / 10) | still mobilising (E / gamma B 10 times lower) |
| Every length times 4 | 13.4 to 14.3 (as a quarter speed at B 0.1 m: 13.9) |

No single term carries the excess. The damping run shows why: a drag of about 90 Pa
over the 0.1 m mechanism raises q by about 2.6 kPa, some 30 times. Any velocity-
proportional resistance raises the mean stress in the mechanism and with it the
frictional strength, so the reading depends on rho c v / (gamma B) and not on the step
count or time. The static value needs that ratio near 0.01; at B 0.1 m that is a push
near 1e-4 m/s (days per run). The scale runs confirm the ratio is what matters.

The same effect probably drives the discharge rate's dependence on E (decision 0091:
-12 to -21 % from 0.1 to 1 MPa): there rho c v is of the order of the silo stresses.

## Decision

1. No scheme change in A.5. A rate-free scheme for slow frictional flow (lower
   dissipation for low-Mach flow, or an implicit quasi-static step) is a multi-session
   study, parked until Sia asks for it.
2. `tests/footing/sand_checks.py` gains `--scale K` (every length times K, same cells).
3. Guidance: quasi-static loads on frictional material are upper estimates unless
   rho c v is at most about 1 % of the stress level of the mechanism (gamma B for a
   footing). Real-scale problems (B of metres, slow loading) get much closer than lab
   scale.

## Consequences

- The sand footing done-when (decision 0057) stays failing; A.5 closes as a finding.
- The experiment patch is not in the repository; the table above records the effect of
  each term.
