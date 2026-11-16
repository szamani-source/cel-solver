# WP 3.4b: C++ discharge rate independent of E (split from 3.4)

State: closed, done-when not met (2026-11-16)
Branch: `claude/project-thread-gcc4kh`

## Goal

Make the C++ T0 discharge rate independent of E (0.1, 1, 10 MPa within 10 %) with
rule 12 (normal Riemann jump capped at yield) while the column collapse still passes.

## Read first

- `docs/decisions/0028`, `0029`, `0032`, `0034`; `tests/manual/test_discharge_stiffness.sh`

## Done when

`tests/manual/test_discharge_stiffness.sh` passes at 10 mm and at 5 mm (0.1, 1 MPa) while
`tests/test_collapse.sh` still passes, both with the same switches.

## Out of scope

- Rule 13 (no one-sweep void in the unsplit step), mu(I).

## Log

- 2026-11-16: One session (decision 0035). Rule 12 in `eulerianSolidLagrange.C`
  (`capNormal`, switch `normalCap yes`, default off). Test drivers read extra
  materialProperties lines from `CEL_EXTRA_MATERIAL`. With it on: 10 mm 78.6 / 73.8 /
  84.8 kg/s/m at 0.1 / 1 / 10 MPa (met); 5 mm 82.6 / 69.7 at 0.1 / 1 MPa (-15.7 %, not
  met; was -29 %). Collapse passes its physics checks (0.50 / 0.35 / 0.25 at a = 0.5,
  1.50 at a = 2); its 2-processor bitwise check fails with the switch on (bisection
  grid, round-off). CI unchanged (switch off). Next: Sia chooses E = 0.1 MPa or mu(I).
