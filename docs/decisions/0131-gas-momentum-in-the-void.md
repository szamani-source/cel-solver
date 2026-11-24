# 0131 Gas momentum in the void (WP F.6)

Date: 2026-11-24. Status: accepted (follows decisions 0125 items 6-9, 0129).

## Context

F.6 starts Stage B: a gas velocity in the void, low Mach, implicit pressure, time step
set by the solid (decision 0125 item 6). Until F.4 the void held the ambient pressure
p0. Sia's rules: ergunWenYu in every test; 2 and 4 processors match serial to 1e-10
with gas mass to 1e-12, cuts through the bed, the free surface and the outlet; report
the cost with gas on against off (F.9 target 1.5x).

## Decision

1. **One model switch, one pressure solve.** `model gasFlow;` in
   `constant/gasProperties` (needs a `freeSurface` solid). Void cells (alpha below
   `alphaVoid`) are no longer held at p0: they solve their own pressure in the same
   Newton solve as the pore gas, with face coefficient rho_f rAU_f on faces between
   void cells and the H/A mass flux as a source. Faces between bed and void keep the
   F.4 rule (the pore side's mobility over its half distance), now against the void
   cell's pressure. Open patches (fixedValue p) give the atmosphere. Gas mass stays
   conservative through the same conservative finish (no separate void bookkeeping).
2. **Momentum in void cells only.** Isothermal ideal gas,
   rho_o (Ug - Ug_o)/dt + div(m Ug) - Ug div(m) - div(mu grad Ug) = -grad p, predictor
   once per step, then `nCorr` pressure corrections (default 1; 2 gives the same
   results to 1e-11 in the tests and costs 25 % more). Cells with grains are not
   solved: they hold the grain velocity plus the pore flux over phi_g (moving walls for
   the void gas). Inside them the pore gas is Darcy-Ergun as in Stage A; F.7 replaces
   that with momentum plus implicit drag.
3. **Rhie-Chow with ddtCorr.** The face flux of H/A takes its time-derivative part
   from the old face flux, not the interpolated old cell velocity. Without it the
   venting column grew a 3 Pa odd-even pressure pattern in the void.
4. **Pressure gradient from the face pressures.** The cell pressure gradient of a void
   cell uses `pFace` (the void cell's own pressure on faces to the bed). With the
   linear face value the void cell above a bed felt half the pressure jump at the
   surface as a gradient that no face flux balances, and its velocity grew without
   bound (15 m/s after a few steps).
5. **Gas weight left out**, as in the pore gas (12 Pa/m against about 16 kPa/m for
   the grains). A balanced (hydrostatic) form can come with F.7 if needed.
6. **Defaults.** Ug: 0/Ug if present (a restart reads it), else no slip on walls,
   zeroGradient where p is fixedValue, constraint patches as they are. Convection
   `divUg` "Gauss upwind"; momentum solver `USolver` PBiCGStab/DILU to relTol 1e-10.

## Results

`tests/test_gas_flow.sh` (about 1.5 min on 4 cores):

- Plane channel, void only, from rest (20 cells across, Re 1.2): Ug within 0.25 % of
  Poiseuille's u_max, mean +0.50 % (the second-order wall error of the cell-centred
  scheme: 0.25 % / 0.5 % at 20 cells, checked by a 1D solve), cross flow 3e-12 m/s; the
  18,404 steps equal the run without gas step by step; gas mass 2e-16.
- Still air, closed box of non-orthogonal cells: |Ug| exactly 0.
- Venting column of F.4 with the void solved: within 0.24 % of Terzaghi and 0.04 % of
  the poreGas run (Stage A reproduced, decision 0125 item 6); void speed 4e-4 m/s.
  Mohr-Coulomb outlet box (material and gas leaving): discharged mass +0.06 % of
  poreGas, gas mass 2e-15, steps equal to the dry run.
- Parallel: channel 2 x 1 and 2 x 2, column 1 x 2 and 1 x 4 (cuts in the bed and at
  the surface), outlet box 2 x 1 and 2 x 2 (cuts through the outlet, bed and surface,
  material leaving): p within 5e-11 of the largest |p|, Ug and U within 7e-11 or
  better, gas mass 1e-15 on every count.
- Cost per step, gas on against off: outlet box (480 cells) gasFlow 2.3x, poreGas 2.0x;
  venting column (28 cells) 4.0x and 2.7x. Above F.9's 1.5x target.

Two things the tests had to allow for, both recorded in the check script:

- The gauge pressure of a low-Mach gas carries round-off of about eps p0 times the
  Courant number of the gas (the mass balance turns a 1e-16 flux error into
  p0 x 1e-16 per unit of flow through a cell). In the channel at 1e-3 Pa that is
  4e-10 of the pressure drop at the solid step of E = 100 Pa, so the channel's solid
  is E = 1e4 Pa (10x smaller step), which puts it at 5e-11.
- The nearly still bed of the column (largest U 4e-4 m/s) needs the pressure to 1e-14
  Pa (1e-16 of p_i) for its velocity to agree to 1e-10; at 1e-12 Pa it agreed to
  1.2e-10.
- The column's own time step differs from its dry run by about 1e-6 (median; the surface
  cell's density under the gas pressure; up to 5e-3 where a write time shortens a
  step): reported, not checked. The exact "same step" checks are the channel and the
  outlet box.

## Consequences

- gasFlow is a superset of poreGas: same fields and outputs plus `Ug`. `poreGas` cases
  are unchanged (F.1 to F.5 tests pass).
- Not yet: drag on grains and momentum in the pores (F.7), restart of the convection
  flux (the first step after a restart uses the flux of Ug, F.9), gas weight, partly
  open (cut) void cells use the full cell volume in the momentum.
- Cost: F.9 should profile gasFlow together with the F.4 items (one drag-law call per
  face, assembling once per iteration, skipping the non-orthogonal correction).
