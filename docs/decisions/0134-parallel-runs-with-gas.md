# 0134 Parallel runs with gas: why they part from serial (WP F.9a)

Date: 2026-11-26. Status: accepted (follows decisions 0080, 0130, 0132).

## Context

F.5 to F.8 found parallel runs with the pore gas parting from serial: the fluidisation
onset run by up to 1e-1 after the bed lifts (decision 0132), pack test 4 on 4
processors by 8e-10 in one surface cell (decision 0130), the fluidised dam break from
its first steps (decision 0133). Sia asked (2026-11-26) why parallelisation fails, and
whether the processor-boundary fields are treated correctly and every operation that
must be global is. F.9a answers that before the rest of F.9 (restart, speed, tutorial,
release).

## What was checked

1. **Global operations.** Every sum, maximum and minimum in `celFluid` and the gas
   coupling in `celEulerian` that feeds a decision or a total: the time step
   (`gMin`, identical step sequences on 1, 2 and 4 processors over 25,370 steps), the
   gas mass (`gSum`), the boundary inflow and outflow (`returnReduce`), the void
   exchange, the convergence test of the pressure iteration (`gMax` of the
   correction), the Anderson least-squares sums (`listCombineReduce`), the void
   velocity report. All are global. None was missing.
2. **Processor boundaries.** New debug option `parallelCheck yes` in
   `constant/gasProperties` (`parallelCheckEvery N`, default 1;
   `src/celFluid/poreGas/poreGasParallelCheck.C`). Before and after each gas solve it
   compares, on every processor cut:
   - 15 cell fields (gas: `p`, `Ug`, 1/A, H/A, e, k, phi_s, phi_g, void flag, q, the
     inertia reaction; solid: `alpha`, `rho`, `U`, `sigma`): the stored neighbour
     values against a fresh exchange. A field changed after its last
     `correctBoundaryConditions` (a stale halo) shows here;
   - 12 face fields (gas mass flux, face pressure, 1/A, H/A flux, old fluxes, inertia
     terms; solid: `phiL`, remap mass): this side's value against the other side's
     (fluxes with opposite sign). A face value worked out from different data on the
     two sides shows here.
   Results (each mismatch over the field's largest value): fluidised bed 1 x 4 and
   2 x 2 through the lift to 0.28 s, 1.2e-14 and 2.6e-14; hopper, `model poreGas`
   and `model gasFlow`, 2 x 2 with the outlet open to 0.5 s, 6.5e-16 and 2.5e-15.
   Round-off everywhere: no stale halo and no one-sided face.

## What the parting is

3. **Before lift: round-off, amplified by the gas.** The bed (8 x 40 cells, sand
   0.2 mm, `model gasFlow`) to 0.03 s, each difference over the field's largest value:

   | Pair | U | p |
   | --- | --- | --- |
   | dry bed in motion (gravity 11 m/s^2), serial vs 2 x 1 and 2 x 2 | 5e-15 | - |
   | gas, serial vs serial with p0 changed by 1e-15 | 2.2e-10 | 1.3e-12 |
   | gas, serial vs serial on a renumbered mesh (`renumberMesh`: only the order of the sums) | 2.5e-10 | 1.4e-12 |
   | gas, serial vs 2 x 2 | 1.8e-9 | 3.0e-11 |
   | gas, serial vs 1 x 4 | 1.7e-9 | 2.7e-11 |

   The solid alone keeps round-off (5e-15). The gas-grain coupling turns a change of
   1e-15 into 2e-10 of U within a few ms: the gas density takes the absolute pressure
   p0 + p, so each step rounds the gauge pressure by about eps p0 = 2e-11 Pa, and the
   bed is nearly still (U is 1e-4 m/s, the differences 1e-13 m/s). Parallel runs, whose
   sums and processor faces round differently every step, land about ten times above
   the serial perturbations, spread evenly over the bed (not at the cuts). Solving the
   pressure to 1e-14 Pa instead of 1e-12, its linear solves to 1e-12 instead of 1e-3,
   and the void gas velocity to 1e-14 changes none of these numbers, so it is not the
   solver tolerance either.
4. **After lift: the bed amplifies any difference.** A serial run with only
   `pTolerance` changed from 1e-12 to 1e-13 Pa (no parallel at all) parts from serial
   by 1e-2 at 0.21 s and O(1) at 0.22 s; 2 x 1 and 1 x 4 part at 0.24-0.26 s.
   Once the effective stress falls to zero the Mohr-Coulomb grains sit at the apex
   and a uniformly fluidised bed is unstable (it bubbles, Anderson and Jackson 1968):
   a difference of 1e-13 grows by about e^18 in 40 ms. Two serial runs cannot agree
   there either. `tests/manual/test_lift_sensitivity.sh` shows it.

## Decision

- No solver change: the parallel code is correct. The rule for parallel checks with
  gas is decision 0080's: cell by cell to round-off (1e-10) where the solution is
  stable (before lift, quiet beds, steady discharge); as a band of outcomes (U_mf,
  rates, run-out) after a switch or an instability.
- `parallelCheck` stays as a debug option for any new gas or coupling code: a stale
  halo or a one-sided face shows within one step. Default off (it costs about one
  extra exchange per field).
- `tests/test_gas_parallel.sh` (about 1 min): `parallelCheck` on the bed (2 x 2,
  1 x 4) and the gas-flow hopper (2 x 2, outlet open): every mismatch within 1e-12;
  serial and parallel bed to 0.03 s within 1e-10 (p) and 1e-8 (U) of their largest
  values (measured 3e-11, 1.8e-9), with the renumbered-serial band printed beside them;
  gas mass to 1e-12 (6e-16).

## Consequences

- The known issues of decisions 0130 (pack test 4 on 4 processors, 8e-10 in one
  nearly empty surface cell) and 0132 (parting after lift) are this round-off,
  amplified where a cell or the bed is at a switch: not parallel errors. Pack tests
  keep their 1e-10 rule only up to such a point and compare outcomes after it.
- The F.8 dam break (decision 0133, its own branch) starts at the Mohr-Coulomb apex
  and parts from its first steps; that looks like the same amplification but was not
  rechecked here. Running it with `parallelCheck yes` will tell.
- Still in F.9: restart of the gas state, cost per step (1.9-3.8x against the 1.5x
  target), tutorial, manual, release.
