# WP F.6: Gas momentum in the void

State: done (decision 0131)
Branch: assigned by the thread

## Goal

A gas velocity field solved in the void (no grains yet): low-Mach, implicit pressure, time step set by the solid.

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- poreGas -> gasFlow model, U_g field, PISO-like loop on fvMatrix

## Done when

- `tests/test_gas_flow.sh`: plane channel against Poiseuille within 1 %; still air in a closed box stays below 1e-10 m/s; time step unchanged from the solid one.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

## Out of scope

- Drag on grains (F.7).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-24: done (decision 0131). `model gasFlow` in gasProperties: void cells solve
  Ug (momentum predictor, `nCorr` PISO corrections, default 1) and their own pressure in
  the poreGas Newton solve (rho_f rAU_f on void faces, H/A flux as a source); bed cells
  keep Darcy-Ergun and hold U_s + q/phi_g for the void. Code in
  `src/celFluid/poreGas/poreGasFlow.C`. `tests/test_gas_flow.sh` (about 1.5 min) passes:
  Poiseuille 0.25 % of u_max (mean +0.5 %), still box exactly 0, steps equal to the dry
  run, venting column 0.24 % of Terzaghi and 0.04 % of poreGas, outlet box discharge
  +0.06 % of poreGas; 2/4 processors 5e-11, gas mass 1e-15. Cost gas on / off 2.3x
  (outlet box), 4.0x (28-cell column). Learned: ddtCorr is needed (3 Pa checkerboard in
  the void without it); the void cell's pressure gradient must use the face pressures
  (pFace) or the cell above a bed accelerates without bound; low-Mach gauge pressure
  carries round-off of eps p0 Co, so parallel checks at tiny pressure drops need a small
  solid step. Next: F.7 (drag on grains, momentum in the pores, the other drag laws).

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
