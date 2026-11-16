# WP 3.4: Silo discharge in celFoam across stiffness, serial and parallel

State: closed, done-when not met (2026-11-16)
Branch: `claude/project-thread-0v2c67`

## Goal

Measure the C++ tank discharge rate at E = 0.1, 1 and 10 MPa with the WP 3.3 outlet and
the collapse code (shear cap and rule 9 on, neither rule 12 nor 13), and check that a
run on two processors gives the serial rate. This settles whether rules 12 and 13 are
needed in C++ (decision 0029).

## Read first

- `docs/decisions/0029`, `0032`, `0033`; `tests/outlet/outlet_checks.py`

## Build

- `tests/test_discharge.sh`, `tests/discharge/discharge_checks.py` (reuses the WP 3.3
  set-up, only E changes). No solver change unless the rate depends on E.

## Done when

Rate at 1 and 10 MPa within 10 percent of the 0.1 MPa rate (now
`tests/manual/test_discharge_stiffness.sh`, NOT MET), and `tests/test_discharge.sh`
(met, in CI): mass balance 1e-10 and no stress in void in every run;
two processors (outlet split
across them) give the serial mass at opening to 1e-10 and the serial rate within 2 percent.

## Out of scope

- Mesh refinement, other outlet widths, rate-dependent friction (mu(I)).

## Log

- 2026-11-16: One session (decision 0034). 10 mm: 74.1 / 65.7 / 70.9 kg/s/m at
  0.1 / 1 / 10 MPa (1 MPa -11.3 %, fails the 10 %); 5 mm on four processors: 78.0 / 55.0
  at 0.1 / 1 MPa, so the rate does fall with E and 10 mm hides it. Two processors give
  the serial rate (74.09) and mass at opening to 8e-16. E sweep moved to
  `tests/manual/test_discharge_stiffness.sh` (fails); CI test keeps the met checks.
  Next: WP 3.4b (rule 12 in C++ with collapse and E sweep as checks).
