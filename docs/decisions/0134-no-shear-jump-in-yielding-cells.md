# 0134 No elastic shear jump beside yielding cells (WP F.8j)

Date: 2026-11-26. Status: accepted as an opt-in (Sia approved the dry-case trial).

## Context

The dry dam break of WP F.8 ran 35-50 % further on 5 mm than on 10 mm cells, and its
run-out fell with E. The face Riemann solver damps shear with the elastic impedance
z_s = rho c_s on every face, also where the material flows plastically, where the
shear stress is fixed by the yield surface and no elastic shear wave carries it.

## Decision

`materialProperties yieldShearJump no` (Mohr-Coulomb or Drucker-Prager; default yes): a
cell is yielding when its stress is within 3 % in sin(phi) of the static Mohr-Coulomb
surface (partly filled cells too). At a face beside a yielding cell the yielding side
contributes only its mu(I) viscosity (min(z_s, 2 eta/delta), with flowImpedance) or
nothing; the sides combine in series; elastic cells keep z_s. The Hancock predictor
drops the elastic stress increment along the flow direction in yielding cells.

## Consequences

- Plain Mohr-Coulomb dry run-out 10 / 5 mm: 1.08 / 1.56 h0 before, 1.85 / 2.59 after.
  Run-out grows, but the mesh dependence (+40 %) stays: shear damping was not the
  whole cause. The rest is the normal jump in slow flow (decision 0135).
- Default unchanged, so every existing test is unchanged.
