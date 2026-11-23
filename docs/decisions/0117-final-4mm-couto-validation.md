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

## Results

To be added when the runs end.
