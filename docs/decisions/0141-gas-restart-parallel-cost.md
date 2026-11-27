# 0141 Gas runs: restart state, parallel outcomes, cost, tutorials (WP F.9)

Date: 2026-11-27. Status: accepted (follows decisions 0080, 0125, 0130, 0133, 0134).

## Context

F.9's done-when (decision 0125 item 8): restart and 2/4 processors match serial on the
F.5 hopper and the F.8 dam break; gas on costs at most 1.5 times gas off per step and
4 processors run at least 2.5 times faster than serial; a tutorial and the manual.
F.9a (decision 0134) found no parallel error and the parting amplified round-off;
F.9b (decision 0135) set the range of validity. Sia (2026-11-27, through the
coordinator): wrap up without serial waiting; the gas tutorials go in
`tutorials/solidFluid/` (decision 0140).

## Decisions and results

1. **Restart state of the gas.** A restart rebuilt the gas from the written fields
   only. The F.5 hopper restarted at 0.2 s still matched (1.8e-13), but the F.8 dam
   break parted by order 1 (p 0.85, U 0.94 of their largest at 0.2 s after a restart
   at 0.1 s; a dry restart matched exactly). Cause: route B with `frictionOnly`
   carries the support ratio p/P with P the total mean stress of the last solve
   (`Pprev`), which a restart reset to the current P. celFoam now writes, with the
   fields at each write, `<time>/uniform/poreGas` (each processor): the gas totals
   (`netIn`, `voidIn`), the last step's length and pressure (the Newton first guess),
   `Pprev` and the lift check's clock (`poreGas::writeState`, read back in the
   constructor). Both cases now restart bit for bit (every written field 0 difference;
   `tests/test_gas_restart_mpi.sh`). Not restarted yet: the old face fluxes of
   `model gasFlow` (pore inertia, `qnOld`, and the void's `phiOld`): the first step
   after a restart of a gasFlow run differs at the level of one step's inertia.
   Gas mass totals no longer restart at 0 (the note of decision 0129).
2. **Parallel: compare cell by cell until a switch flips, outcomes after.** The hopper
   (plain Mohr-Coulomb) on 2 and 4 processors matches serial to 1e-10 (p, U) with gas
   mass to 1e-12. The dam break with the test's cuts (x = 0.1 m, and y = 0.1 m on 4)
   differs from serial by 2e-2 at the first write, 0.002 s, in the cells at the
   column's free face (not at the cut). Its start is degenerate: fully supported, the
   effective stress is zero, the apex of the yield cone, where round-off decides
   whether a cell yields. A serial run with E changed by 1e-13 departs as fast (U 0.25
   of its largest at 0.02 s, against 0.29 for 2 processors). With gas off the column
   parts too, later, and the parting starts in the version 1 free-surface rules. Cut across the column (1 x 2), each rule alone
   with gas off: plain Mohr-Coulomb 7e-13 at 0.04 s, plus mu(I) 2e-13; plus
   `flowImpedance` 1.5e-1, `fullBlend` 9e-6 (2e-4 at 0.02 s), `overfillStress` 5e-2.
   Step by step (`flowImpedance`): 1e-16 at the start, growing smoothly to 1e-9 by
   0.011 s, then 1.5e-3 in one step in surface cells at the top of the column (not at
   the cut): a switch (yield or impedance on/off) flipped by round-off. `parallelCheck`
   (decision 0134) on the dam break, 2 processors cut through the column, every fifth
   step to 0.4 s: largest processor-boundary mismatch 7.5e-16. So the parallel code is
   consistent and the column amplifies round-off, as decision 0080 found for the T0
   rules. Outcomes, 10 mm, to 9 sqrt(h0/g):

   | Run | Run-out from the gate | Front speed |
   | --- | --- | --- |
   | serial | 5.20 h0 | 2.140 m/s |
   | serial, E x (1 + 1e-13) | 5.30 | 2.138 |
   | serial, E x (1 - 1e-13) | 4.88 | 2.147 |
   | serial, E x (1 + 1e-12) | 5.73 | 2.155 |
   | 2 processors (cut at x = 0.1 m) | 5.50 | 2.062 |
   | 4 processors (also y = 0.1 m) | 5.10 | 2.114 |

   A change of 1e-13 in E moves the run-out by up to 8 %, so the run-out on 2 and 4
   processors is inside the round-off band of serial runs. The front speed on 2 and 4
   processors is 1 to 4 % below the serial runs, whose own spread is 0.8 %: not
   explained beyond round-off growth (the boundaries are consistent); four serial
   samples do not size that band well. So the done-when's 1e-10 holds for the hopper
   but not for the dam break: `tests/test_gas_restart_mpi.sh` checks the hopper to
   1e-10 and reports how long the dam break stays within 1e-10;
   `tests/manual/test_fluidised_collapse.sh` keeps its 1e-10 check over the whole run
   and fails it. Comparing the dam break's processor counts by outcome (run-out, front
   speed) instead would relax that test, so it is Sia's call (asked 2026-11-27).
3. **Cost and speed-up** (`tests/manual/test_gas_scaling.sh`: the F.5 hopper refined to
   64 columns, 8,704 cells, 0.2 mm grains, T0 rules, outlet open from the start, 200
   steps, idle 4-core machine): dry 45.5 ms per step, gas 88.2 ms (1.94 times: fails
   the 1.5 limit); gas on 2 processors 1.65 times serial, on 4 processors 2.96 times
   (passes 2.5, 2,176 cells per processor). The gas cost is the Newton iteration: the
   default `pTolerance` (1e-14 p0, 1e-9 Pa) takes 11.4 solves per step. With
   `pTolerance` 1e-6 / 1e-4 / 1e-2 Pa: 7.3 / 4.8 / 2.9 solves, 1.68 / 1.50 / 1.38 times
   dry; gas mass conserved to 5e-14 in every case (the conservative finish), pressures
   within 0.04 Pa of 49 Pa at 1e-4 (the T0 rules' round-off flips give the same order).
   The default is Sia's call (asked 2026-11-27); until then it stays 1e-9 Pa.
4. **Tutorials** in `tutorials/solidFluid/` (decision 0140 layout: Setup, Numerics, How
   to run, Expected results, figures):
   - `airImpededHopper`: the F.5 hopper with gas and dry. Default 0.2 mm grains: 54.6
     against 74.4 kg/(m s) dry (-27 %), about 2 minutes. At 0.1 mm (-64 %) the lift
     check warns from 0.23 s (the gas carries about 2 % of the grains above the outlet
     with the T0 rules; the F.9b validity test uses plain Mohr-Coulomb, 0.06 %), so
     0.2 mm is the default and the README says why.
   - `fluidisedDamBreak`: the F.8 case, fluidised and dry, 10 mm by default (5 mm, the
     validated mesh, by argument); `makeColumn.py` writes the start fields. Its 10 mm
     run is bit for bit the test case of `collapse_cases.py`.

## Consequences

- Every gas run restarts exactly; `uniform/poreGas` must be kept with the fields.
- The cell-by-cell parallel match holds for the hopper and for the first steps of the
  dam break; whether the dam break is compared by outcome after that is open.
- The 1.5 times cost target is met only with a looser `pTolerance`.
