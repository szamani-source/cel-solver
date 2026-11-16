# 0034 C++ discharge rate across stiffness and in parallel (WP 3.4)

Date: 2026-11-16. Status: accepted (WP 3.4, done-when not met; follow-up WP 3.4b).
Measures the C++ tank of decision 0033 (outlet, rule 10b) with the collapse code of
decision 0032 (shear cap and rule 9 on, neither rule 12 nor 13, decision 0029) at
E = 0.1, 1 and 10 MPa, and on two processors.

## Context

Decision 0029 left one question for C++: with the code that runs the column collapse, is
the discharge rate independent of E? If not, a cap of the normal jump (rule 12) or void
spreading (rule 13) is needed, with the collapse as its check. The prototype without
those rules fell from 81 to 54 kg/s/m between 0.1 and 10 MPa (decision 0027).

## Result

Tank W 0.4 m, H 0.6 m, D 0.06 m, walls 20 deg, opens at 0.3 s; rate over 97 to 88 % of
the mass at opening (as decision 0033).

| Run | 0.1 MPa | 1 MPa | 10 MPa |
| --- | --- | --- | --- |
| C++, 10 mm, serial | 74.1 (C 0.89) | 65.7 (-11.3 %) | 70.9 (-4.4 %) |
| C++, 5 mm, 4 processors | 78.0 (C 0.94) | 55.0 (-29 %) | not run |
| Prototype, 10 mm, no rules 12/13 (0027) | 81 | 73 | 54 |

- The 10 mm numbers look nearly E-independent, but not monotone, and the 5 mm run shows
  a clear fall with E (-29 % at 1 MPa; steady over the window: 77.3 / 77.4 against
  55.2 / 54.5 on its two halves). The coarse mesh hides the dependence; the C++ code
  behaves as the prototype did without rules 12 and 13.
- The 10 mm rate does not depend on the window much (0.1 MPa 72 to 76, 1 MPa 61 to 68
  over six windows between 97 and 80 %), so the -11 % at 1 MPa is real, not noise.
- Two processors (simple split at x = 0.2 m, through the middle of the outlet):
  mass at opening 8e-16 from serial, rate 74.09 against 74.09 kg/s/m, mass balance
  7e-15. Unlike the collapse (decision 0032), the discharge rate is not sensitive to
  round-off switching.
- Mass balance 4e-14 and zero stress in empty cells in every run.

## Decision

1. WP 3.4 closes with its done-when not met: the rate depends on E with the collapse
   code. Not weakened: the E sweep stays a pass/fail check (within 10 %), run by hand as
   `tests/manual/test_discharge_stiffness.sh` (it fails now; about 8 minutes at 10 mm)
   and is the done-when of WP 3.4b.
2. `tests/test_discharge.sh` (in CI, about 50 s) checks what is met: mass balance, no
   stress in void, two processors against serial (rate within 2 %).
3. Any E check of the discharge is run at 5 mm as well as 10 mm (at least 0.1 and
   1 MPa, four processors, about 5 minutes): 10 mm alone hides the dependence.
4. WP 3.4b: try rule 12 (normal jump capped at yield) in the C++ Lagrange step behind a
   switch, with both the E sweep and `tests/test_collapse.sh` as checks. The prototype
   found no variant serving both (decision 0029); if C++ repeats that, the fallback to
   put to Sia is to run T0 at a soft numerical stiffness (E = 0.1 MPa, where the rate is
   C 0.9 to 0.94 and the collapse run-out is the prototype's), as explicit codes use
   mass scaling for quasi-static flow, or a rate-dependent friction (mu(I)).

## Consequences

- The C++ T0 tank runs serial and in parallel with the same rate; it is usable now at
  E = 0.1 MPa.
- Stiffer materials under-predict the discharge rate, as they under-predict the
  collapse run-out (decision 0032). One cause (the normal Riemann jump acting as a
  viscosity rho c h, decision 0028) is behind both.
