# 0027 Outlet faces pass void (rule 10b); rate-law checks at 5 mm with averaged corner speed

Date: 2026-11-15. Status: accepted for the prototype (WP 3.3q).

## Context

Decision 0026 left the T0 discharge rate at a third to a half of Beverloo (C 0.32 to
0.54 at 5 mm, exponent 2.26), dependent on E (38 / 17 / 6 kg/s/m at E = 0.1 / 1 / 10
MPa, 10 mm) and on the mesh (23 percent between 5 and 2.5 mm), and blamed the
dissipation of the acoustic Riemann solver. Every attempt to reduce the jump terms made
the flow noisy.

## Finding: the outlet was intermittently shut

Tracing one outlet cell step by step (E = 0.1 MPa, 10 mm, 0.25 s after opening) showed
its vertical stress jumping between 0.2 and 9.6 kPa from one step to the next while the
cells above held 0.1 to 1 kPa. In the y sweeps where it jumped, the outlet face had
velocity exactly 0. Cause, in `void2d.void_sweep`: when the outlet cell is partly filled
(al 0.97 to 0.9995, which happens every few steps as material stretches there), rule 2
leans its slab towards the fuller neighbour, which is the cell above, because the ghost
beyond an outlet is empty (rule 10). The slab then does not touch the outlet face, and
the boundary code gave a non-touching boundary face the closed-wall velocity 0. With the
face held still and material arriving from above, the cell's void closed and its
material was compressed by up to 3 percent in one sweep. The outlet row rang and held
the column up. The jump terms were not the main cause.

## Decision 1: rule 10b

An outlet face whose cell does not touch it takes the free one-sided velocity anyway (the
same relation as when it does), for the normal and the tangential system. The void at
the bottom of the cell passes out first; material leaves once the swept length exceeds
the void (`out_overlap`, unchanged). Walls keep velocity 0. Switch: `void2d.OUTLET_VOID`
(default on). Nothing else changed: closed outlets, walls and interior faces are as
before, so the void, collapse and settling tests are unaffected.

Result, E = 0.1 MPa, H = 0.6 m. 10 mm: rate from 97 to 88 percent of the mass (D = 0.10 m:
95 to 86); 5 mm: the rate-law check (`--rate-check`, 97 to 88 percent):

| D (m) | 10 mm: rate (kg/s/m), C | 5 mm: rate, C | Before (5 mm, 0026) |
| --- | --- | --- | --- |
| 0.04 | 33.9, 0.75 (4 cells across) | 44.7, 0.99 | 14.3, 0.32 |
| 0.06 | 81.4, 0.98 | 82.2, 0.99 | 36.6, 0.44 |
| 0.08 | 125.2, 0.98 | 135.1, 1.06 | 68.4, 0.54 |
| 0.10 | 174.0, 0.98 (95 to 86 percent) | | |

5 mm: exponent 1.59, fitted C 1.035, 0.1 s windows within 1.7 to 6.1 percent, head
independence 0.992 (H = 0.9 m: 81.5). 10 mm without the under-resolved 0.04 m outlet:
exponent 1.47 (0.06 to 0.10 m), C 0.98. D = 0.06 m at 10 against 5 mm: 81.4 against
82.2 (1 percent; 0026 had 23 percent between 5 and 2.5 mm).
The bottom outlet row now carries 0 to 0.2 kPa and the centre outlet speed is 1.2 m/s
(free fall from an arch of D/2 gives 0.77 m/s at the centre per S5's sqrt(2 g R); the
material is still accelerating as it leaves).

C of about 1.0 is below the mu(I) continuum value 1.45 and inside the 1.0 to 1.9 band at
5 mm. A lower C is expected: T0's friction (phi 35 deg, sin phi = 0.57) is higher than
the mu_s = 0.32 to 0.38 of the reference runs, and S2 reports slower discharge for higher
friction (`docs/spec/silo-reference.md` section 3). Not tuned.

## Decision 2: the rate-law checks run at 5 mm, in parallel

At 10 mm the narrowest T0 outlet is 4 cells wide and gives C 0.75 against 0.98 for the
others, so the 10 mm exponent (1.90) measures resolution, not the scheme. The checks of
the WP 3.3q done-when (bands unchanged) now run at 5 mm (8, 12 and 16 cells across), the
production mesh of the spec, with the four cases in parallel processes
(`discharge2d.py --rate-check`, `tests/test_proto_silo_rate.sh`). The outlet checks
(`--check`, CI) stay at 10 mm and gain a guard: C at D = 0.06 m at least 0.8 (0.98 with
rule 10b, about 0.3 without).

## Decision 3: stagnant corners are judged on the mean velocity

With damping off after opening, the soft bed (E = 0.1 MPa, c = 8.7 m/s) rings. In the
5 cm corner squares the instantaneous largest speed averaged 1.9 cm/s with peaks of 4 cm/s
(2 to 4 percent of the outlet speed, 10 mm), while the largest mean vertical speed over
0.4 s was 0.52 cm/s (0.6 percent). Making slow capped faces stick (`VSTICK` = 0.01 m/s)
changed nothing, so it is elastic vibration, not creep along faces. The corner check is
about flow, so it now uses velocities averaged over a further 0.3 s of discharge
(`averaged_pattern`), at least one period of the bed's first elastic mode (4 H / c_p =
0.28 s at H = 0.6 m, E = 0.1 MPa). A 0.1 s average still read 1.5 percent at 5 mm
(instantaneous 4.2). The instantaneous ratio is still printed.

## Still open (split to WP 3.3r)

- Stiffness: at 10 mm, D = 0.06 m, rates 81 / 73 / 54 kg/s/m at E = 0.1 / 1 / 10 MPa
  (at 10 MPa the 0.1 s windows vary by 11.5 percent). The done-when's 10 percent across
  E is not met. Re-measured with the outlet fixed, the dissipation knobs still do not
  help: rule 11 (low-Mach jump terms) gives 111 kg/s/m at 1 MPa but corner speeds 16
  percent of the outlet speed, and with a floor of 0.1 at 10 MPa 186 kg/s/m (C 2.25) with
  the corners moving faster than the outlet; halving only the normal or only the shear
  jump term (`DISS_N`, `DISS_T`, new study knobs) at 10 MPa raised the rate 5 to 10
  percent and left it as noisy. The E dependence is also not a coarse-mesh dissipation
  of order rho c h: at E = 1 MPa the rate falls with refinement (73 kg/s/m at 10 mm,
  61 at 5 mm) while at 0.1 MPa it does not move (81, 82). Cause unknown.
- Fine mesh: D = 0.06 m, rate over 0.2 to 0.5 s after opening: 81.5 kg/s/m at 5 mm,
  93.3 at 2.5 mm (14 percent, C 0.98 against 1.13; 10 against 5 mm agreed to 1
  percent). The done-when's 10 percent between 5 and 2.5 mm is not met; C stays in band.
- mu(I) was not needed to reach the band and was not built.

## Consequences

- The C++ outlet (WP 3.3) must give a boundary face next to a cell that does not touch
  it the free velocity at an outlet, not the wall velocity.
- Rule 11 (`low_mach`) stays off and in the code as a study knob.
- `tests/test_proto_silo_rate.sh` is in CI: 1407 s here on four cores shared with two other
  runs; expect 20 to 25 minutes, the longest test in CI.
