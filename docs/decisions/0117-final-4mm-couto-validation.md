# 0117: Final Couto validation on the 4 mm mesh, and checkpoints for the test driver

Date: 2026-11-23. Thread "Validate solver on 4 mm mesh".

## Context

Sia asked for a final verification and validation on the fine ("ideal") 4 mm
Couto mesh of decision 0115's addendum, against Couto et al. (2013) and Zheng
and Yu (2015). That mesh had 1.5 mm cells at the outlet edge, below Sia's
3-5 mm outlet rule. The runs take about 6 h on 4 cores, and cloud containers
were restarted three times during the first attempts (once 3 minutes after the
last tracked task ended, once with two tracked tasks running), which kills
the run. The test driver (tests/void/voidCheck) could not restart.

## Decisions

1. **Outlet columns.** The hopper block keeps the cylinder's 58 columns (the
   mesh stays conforming), but its outlet edge is graded evenly
   (`outlet_wall_ratio`, blockMesh `edgeGrading`); the top edge keeps the
   cylinder's wall grading and blockMesh blends the two down the cone.
   Outlet cells 3.1-3.2 mm, 4 mm at the walls and through the hopper,
   20,398 cells, checkMesh OK (non-orthogonality 34.9 deg, aspect 3.3).
   Knob `CEL_OUTLET_WALL_RATIO` (gap_case.py, hopper_pressure_checks.py).
2. **Checkpoints.** celVoidCheck writes the fields and the solid's restart
   state (decision 0081) every `checkpointEvery` seconds (voidCheckDict) and,
   started from a later time (`startFrom latestTime`), skips the fill, keeps
   history.dat up to that time and appends. The nonlocal fluidity `g` is now
   read when present. `tests/test_voidcheck_restart.sh`: the coarse Couto silo
   restarted at 4 ms matches the straight run bit for bit in every cell; the
   discharged-mass counter in history.dat differs in the 8th digit (round-off
   in the counter only, not used by the solver).
3. **Run script.** `tests/hopper/run_final_4mm.sh mc|dp`: checkpoints every
   0.05 s, resumes from the latest complete checkpoint when run again, keeps
   the latest two, and mirrors the case to /mnt/project-files every 5 min (a
   lost container is restored from the mirror). `writeFormat binary` so the
   restart is exact.
4. **What is compared** (both papers read, 2026-11-22): Couto's gauges measure
   the normal wall pressure (cylinder: horizontal; GA0 on the hopper wall:
   normal to the cone); Table 1's discharge rows are the **maximum** during
   discharge. Zheng and Yu Fig. 11 is one snapshot at 0.5 s after opening, the
   cylinder from sigma_xx and the hopper rotated to the wall normal. Ours: the
   same stresses (`hopper_case.wall_pressures`). Earlier decisions compared our
   discharge **mean** with Couto's maximum; the final comparison uses our
   maximum over the discharge writes against Couto, and t = 1.0 s (0.5 s after
   opening) against Zheng and Yu, with the mean reported too.

## Results (2026-11-23)

Both runs completed (28,200 s and 28,100 s on 4 cores, dt 3.50e-6 s fixed,
set by the axis cells at the outlet: in-plane 1.6 mm and hoop length 1.6 mm,
so 1.13 mm, over the P-wave speed 129.6 m/s at maxCo 0.4). No leak.
`tests/hopper/final_report.py OUT MC_CASE DP_CASE`:

| Quantity | Reference | MC | DP |
|---|---|---|---|
| Discharge rate, kg/s | Beverloo 110-131 | 118.1 (inside) | 118.2 (inside) |
| GA0 at rest, kPa | Couto 8.61-9.61 | 13.05 (+36 %) | 13.04 (+36 %) |
| GA1 at rest, kPa | Couto 4.00-4.49 | 4.10 (inside) | 4.11 (inside) |
| GA2 at rest, kPa | Couto 4.07-4.46 | 4.64 (+4 %) | 4.64 (+4 %) |
| GA3 at rest, kPa | Couto 3.67-4.16 | 4.04 (inside) | 4.04 (inside) |
| GA4 at rest, kPa | Couto 3.17-3.62 | 3.45 (inside) | 3.45 (inside) |
| GA5 at rest, kPa | Couto 2.78-3.24 | 2.79 (inside) | 2.79 (inside) |
| GA6 at rest, kPa | Couto 2.47-2.66 | 2.08 (-16 %) | 2.08 (-16 %) |
| GA7 at rest, kPa | Couto 1.39-1.69 | 1.30 (-6 %) | 1.30 (-6 %) |
| GA0 discharge max (mean), kPa | Couto max 18.65-24.77 | 29.04 (24.33); +17 % | 28.67 (23.89); +16 % |
| GA1 discharge max (mean), kPa | Couto max 7.14-8.23 | 5.92 (4.90); -17 % | 6.30 (5.16); -12 % |
| GA2 discharge max (mean), kPa | Couto max 5.81-6.45 | 5.20 (4.20); -11 % | 5.13 (4.28); -12 % |
| GA3 discharge max (mean), kPa | Couto max 4.42-5.59 | 5.01 (3.81); inside | 5.00 (3.82); inside |
| GA4 discharge max (mean), kPa | Couto max 3.21-4.01 | 4.56 (3.42); +14 % | 4.23 (3.37); +5 % |
| GA5 discharge max (mean), kPa | Couto max 2.88-3.72 | 3.51 (2.61); inside | 3.58 (2.64); inside |
| GA6 discharge max (mean), kPa | Couto max 3.35-3.63 | 2.88 (2.03); -14 % | 3.09 (2.00); -8 % |
| GA7 discharge max (mean), kPa | Couto max 1.43-1.74 | 1.95 (1.19); +12 % | 1.84 (1.35); +6 % |
| Peak at 0.5 s after opening, MC | Abaqus 27.5 kPa at z 0.46 m | 27.7 kPa at z 0.45 m | |
| Peak at 0.5 s after opening, DP | Abaqus 27.5 kPa at z 0.46 m | 24.5 kPa at z 0.44 m | |
| Largest speed after opening, m/s | | 9.90 | 9.85 |

Discharge writes used: MC 16; DP 16

`hopper_pressure_checks.py` (mean in discharge against Couto's maxima,
+-15 % bands of decision 0100) reports 4 failures per model: GA0 and GA6 at
rest, GA1 and GA2 in discharge. Compared as the papers report (maxima), GA2
moves inside 15 % and GA1 is 12 % (DP) and 17 % (MC) low.

## Decision on the results

Accepted as the solver's validation evidence and added to the manual's
Validation page, with the gaps stated: rate inside Beverloo for both models,
the Abaqus transition peak matched (MC 27.7 against 27.5 kPa), cylinder at
rest inside the measured range at GA1 to GA5, every discharge load cell within
17 % of Couto's maximum (Sia accepts about 20 % for now, 2026-11-18). Known
gap: GA0 at rest 36 % high with the gravity-on start (no filling history).
Compared with the 12 mm-row mesh (decision 0115, wall-gap start, means) the 4
mm cylinder does not raise GA1/GA2, so that deficit is not cylinder
resolution. The time step is set by the axis cells at the outlet; a wider
axis column there would lengthen it (next mesh).
