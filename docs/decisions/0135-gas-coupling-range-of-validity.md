# 0135 Range of validity of the gas coupling; lift check (WP F.9b)

Date: 2026-11-26. Status: accepted (follows decisions 0125, 0127-0134).

## Context

After the bed in the fluidisation onset case lifts (inflow above U_mf), the state is
unphysical in every run, serial included (decision 0134 follow-up): pressure drop
swinging to 8-10 bed weights, gas at 60 m/s, one crash in the drag routine. Sia asked
(2026-11-26) how two-fluid (Euler-Euler) models handle it. Findings, from OpenFOAM's
twoPhaseEulerFoam source and the literature (report:
https://claude.ai/artifact/7AAVYPS4TLUHsL9Tgf7WKf):

- Two-fluid models carry a solids fraction phi_s free between packed and dilute. Drag,
  the gas equation (one in every cell) and the grain pressure (Johnson-Jackson
  frictional pressure diverging at the packing limit, kinetic theory in the dilute
  range) are all functions of it; drag is implicit in both phases (partial
  elimination); residual floors guard the dilute limit.
- Our grains stay at their packed density. A stretched bed opens void inside its cells
  (decision 0021 rule 4) instead of loosening, and the gas sees the permeability of the
  packed part (decision 0125 item 4) until a cell falls below alphaVoid and becomes open
  void. A diagnostic run of the onset bed (8 x 40, 1.4 U_mf, main at e74bf42) to 0.45 s:
  full cells at phi_s 0.620 throughout; the bed lifts as a plug on a gas gap; at 0.44 s
  the plug moves up at 2.6 m/s on a cushion of 34 bed weights; material phi_s reaches
  0.83 (no packing limit).

Sia's decision (2026-11-26): the scope is packed beds; a fluidised system is the job of
a two-fluid model. Accept that the gas coupling is valid only while the bed does not
expand, and state where it is valid. The two-fluid packages proposed in the report
(G.0-G.4) are not adopted.

## Decision

1. **Range of validity.** The gas coupling is valid while the grains' own stress
   carries their weight less the gas's share:

   | Situation | Valid | Evidence |
   | --- | --- | --- |
   | Pore pressure in a packed or consolidating bed | yes | pack tests 1, 3 |
   | Gas through a packed bed up to U_mf | yes | U_mf -3.05 % of Wen-Yu, pack test 2 |
   | Air-impeded discharge of fine powder | yes | pack test 4, lift check quiet (max 0.06 %) |
   | Gas in the void around the material | yes | decision 0131 |
   | Partial gas support, grains in contact | yes while the lift check is quiet | lift check |
   | Uniform settling of a dilute slab as a block | approximately | +7.6 % of Richardson-Zaki, decision 0132 |
   | Gas-supported (fluidised) dam break | not validated | run-out 60 % short, decision 0133 |
   | Bed above U_mf: lifted, expanded, bubbling, slugging; dilute suspensions; conveying | no | two-fluid model needed |

2. **Lift check** (`liftCheck yes` by default in `constant/gasProperties`,
   `src/celFluid/poreGas/poreGasLift.C`). Each step, each cell with grains gets the force
   the solid momentum uses, -alpha sum(p_f S_f) over its open faces, and is counted as
   carried when its upward part is at least `liftRatio` (1) times the cell's weight. The
   carried share of the grain mass is printed at write times ("Gas lift: ...") and a
   warning is given once it has stayed above `liftWarn` (1 %) for `liftTime` (0.01 s).
   The time filter keeps the opening of an outlet (a 3.4 % spike for a few steps in the
   fine hopper) from warning. Cost: one face loop per step. Immersed-wall faces are left
   out of the force (a wall cell may read slightly low).
3. **No drag guard now.** The crash in `inertialGasDrag::slip` happens only after lift,
   which is outside the range; the warning comes about 10 ms before the bed lifts.

## Results

`tests/test_gas_validity.sh` (about 1.5 min):

- Onset bed (8 x 40, ramp to 1.4 U_mf over 0.28 s): share 0 % at every write to 0.16 s,
  1.25 % at 0.18 s, 87 % at 0.20 s; the warning says carried since 0.189 s, between 0.9
  and 1.1 U_mf (0.180 to 0.220 s); 82 % at 0.26 s.
- Fine-powder hopper with gas (pack test 4 case, outlet opens at 0.15 s): share at most
  0.06 %, no warning.

## Consequences

- Users get a clear message when a run leaves the range; results after it are not to be
  trusted (decision 0134's "post-lift" numbers stand as examples).
- The F.8 dam break (fully gas-supported by design) will warn: it sits at the edge of
  the range and is reported as not validated.
- If fluidised beds become a goal, the report's G.0-G.4 plan (one gas equation on the
  cell-average phi_s, Johnson-Jackson pressure, partial elimination, residual floors,
  checked against twoPhaseEulerFoam) is the starting point.
