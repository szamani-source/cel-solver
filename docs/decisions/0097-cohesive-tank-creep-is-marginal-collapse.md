# 0097: The cohesive tank's creep is a marginal plug collapse, slowed by the rate effect

Date: 2026-11-19 (accuracy package A.6, thread "Fix sand footing bearing capacity")
Status: finding recorded; no code change.

## Context

Decision 0054 (WP 3.5) found that a T0 tank stopped behind a cohesive arch keeps
creeping out at 0.5 to 2 kg/s/m. Sia asked to find where errors come from.

## Findings (tests/cohesion/cohesion_checks.py --case, H 0.3 m, 10 mm, T0 rules, 1 s open)

Rate over the last 0.15 s, kg/s/m:

| c (Pa) | D (mm) | E 1 MPa | Note |
| --- | --- | --- | --- |
| 500 | 40 / 60 / 80 / 100 | 0.07 / 0.19 / 0.85 / 0.89 | Jenike B_c 109 mm; 120 mm flows |
| 1000 / 2000 | 100 | 0.11 / 0.06 | B_c 218 / 435 mm |
| 500 | 100, E 0.25 / 10 MPa | 14.1 / 1.8 | |

At c 500 Pa, D 100 mm, a plug 0.13 m wide, wider than the outlet, slides down the full
height of the fill at about 7 mm/s. Two vertical bands near x = 0.125 and 0.275 are at
yield (mobilised Mohr-Coulomb ratio 1.00) from the floor to the top, and its underside
tears off at the outlet. So the creep is a collapse mechanism that exists at these
widths, not a leak. Its speed is set by the velocity-proportional resistance of the
scheme (decisions 0095, 0096). It dies away with reserve: the creep falls tenfold
at c 1000 Pa and is a trickle at 40 mm. Jenike's B_c is a criterion for the arch and
does not cover this plug, so a marginal tank can creep below it.

## Decision

1. No code change. The creep is the slow form of a collapse that the rate effect
   (decision 0096, parked) turns from an acceleration into a steady motion.
2. The cohesion test keeps its 5 % stop threshold; it is still a fair split between
   stop (under 1 % of Beverloo at 80 and 100 mm) and flow (above 25 %).

## Consequences

- A stopped cohesive tank within about 30 % of its critical width may creep; with a
  rate-free scheme it would likely either hold or flow. Treat stop/flow near the
  critical width as uncertain by that margin.
- Known issue "Cohesive tank creep" now points here.
