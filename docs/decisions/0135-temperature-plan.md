# 0135 Temperature: an energy equation for grains and gas (plan, T series)

Date: 2026-11-26. Status: accepted (Sia chose "Go: T.0 to T.8 as planned",
2026-11-26 20:31 UTC). Plan doc: https://claude.ai/artifact/4BPFHQjq6ijFqRRpQ4sfj7.
References and values: `docs/spec/thermal-references.md`, `tools/thermal_reference.py`.

## Context

Sia (2026-11-25, 2026-11-26): once the gas phase is in, add temperature (the energy
equation), before the HTML manual is rebuilt and the version 2 release. The gas work
(decisions 0125 to 0134) gives gas fractions, a pore and void gas velocity and its mass
fluxes; the gas is isothermal at the case temperature.

## Decision

1. **Own library `src/celThermal`**, called from the celFoam driver after the gas
   solve; selected in `constant/thermalProperties` by
   `model none | equilibrium | twoTemperature` (`none` or no file: every v1 and v2
   result unchanged bit for bit).
2. **Equations**: the volume-averaged two-temperature model (Nield and Bejan; Kuipers et
   al. 1992), written in `docs/spec/thermal-references.md`. Closures: Gunn (1978) for
   grain-gas heat transfer with a_s = 6 phi_s/d; Zehner and Schluender (1970) stagnant
   bed conductivity split between the phases as in Kuipers et al.; optional plastic
   heating with a fraction beta (default off). `equilibrium` carries one mixture
   temperature (the local-equilibrium limit).
3. **Grain temperature moves with the material** through the remap, mass weighted like
   the cap's p_c (decision 0113), so energy is conserved exactly; void carries no
   grain energy; outlets take energy out with the mass; the solid inlet brings material
   in at a set temperature.
4. **Gas temperature is carried by the gas mass fluxes** already computed (Darcy flux
   in the bed, Ug in the void), so a uniform temperature stays uniform under any flow.
5. **Numerics**: conduction explicit (stable step about 20 s on 4 mm cells, checked
   every step, sub-cycled if ever needed); grain-gas exchange implicit and exact per
   cell (for 0.1 mm grains the gas relaxes in 1e-5 s, the solid step).
6. **Walls** in `0/T` and `0/Tg`: fixedValue, zeroGradient (adiabatic, the default),
   heat flux, convective (h, T_inf); an optional near-wall air gap (Sullivan and
   Sabersky's 0.085 d) for moving beds against walls.
7. **Gas density and viscosity with temperature** (rho_g = p/(R T_g), Sutherland mu and
   k) only in T.7: it changes the gas pressure solve, so it waits until the fluidised
   dam-break research thread (F.8) has finished. Before T.7 the gas properties stay at
   the case temperature.
8. **Tests**: one closed form or published case per package (table in
   `docs/work-packages.md`); every test also on 2 and 4 processors against serial
   (1e-10), energy conserved to 1e-12 where no heat crosses the boundary (decision 0125
   item 8 applied to temperature).
9. **T.6 validation case changed from the plan**: Sullivan and Sabersky (1975) could not
   be read (the thesis is a scan without text, the paper is paywalled). Adapa et al.
   (2024, arXiv 2403.19892) give measured bed conductivity and near-wall gap for
   flowing particles between heated plates, the slug-flow solution with the gap, and
   measured heat transfer coefficients; T.6 uses them.

## Consequences

- One field per cell (T) for `equilibrium`, two (T, Tg) for `twoTemperature`, plus
  their restart data. No new linear solve before T.7 (conduction explicit, exchange
  per cell).
- Out of scope: thermal expansion and stress, temperature-dependent friction or
  cohesion, radiation beyond the Zehner-Schluender form, moisture and drying, phase
  change, reactions.

## Implementation notes

- T.1 (2026-11-26): class `heatTransfer` in `src/celThermal`. Face flux
  k_f |S| [(T_N - T_P) delta + c . grad T] with OpenFOAM's over-relaxed
  `nonOrthDeltaCoeffs` and `nonOrthCorrectionVectors` and a Gauss gradient with
  linear face values, written out by hand (no `fvSchemes` entries needed). Two
  OpenFOAM details matter: a processor patch field holds the neighbour cell values,
  not face values, so the gradient interpolates them; and `fvPatch::delta()` is only
  the normal part on non-coupled patches, so Dirichlet faces use Cf - C for the
  non-orthogonal correction. With these a linear temperature on 27 deg cells is
  exactly steady and 2/4 processors match serial to 6e-15. k_f = k min(a_P, a_N)
  (a = material fraction), so void neither stores nor conducts heat; how
  temperature crosses the free surface with the material is T.2.
- T.2 (2026-11-26): the solid's remap carries T as a "passenger"
  (`eulerianSolid::addPassenger`, the same mass-weighted path as the cap's p_c), so
  sum m c T is conserved to round-off; outlets carry out m c T of the cell, open
  inlets bring the inlet patch value of `0/T`; heatTransfer counts both from
  `remapMass()` and prints them, so the balance is an independent check of the
  remap. Default limiter `superbee` (as the history fields, decision 0017): with
  `MC` a front carried 50 cells at the small material Courant number of the solid
  step spreads over 5 cells, with superbee 2. Faces between material and void take
  the material side's T in the gradient. Immersed STL walls and bodies are not yet
  supported: the small-cell mixing (decision 0047) would need T too; split as T.2b.
  A plain Mohr-Coulomb stream landing in a box parts serial and parallel runs of the
  solid itself after about 0.2 s (amplified round-off, decision 0134); T matches to
  2e-15 until then and to 2e-14 across processor cuts in a carried block.
- T.3 (2026-11-26): walls in `thermalProperties` (`walls { patch { h; Tinf; } }` or
  `{ q; }`), as the solid's `walls` live in materialProperties, rather than new
  `0/T` patch types: no new boundary-condition classes to register, and the film
  1/h is put in series with the half cell 1/(k delta) so the slab matches the
  Carslaw-Jaeger series within 0.02 % at Bi 10. A wall at a fixed temperature with
  a contact resistance (the moving-bed air gap of T.6) is the same with h = 1/R.
  `zehnerSchlunder` takes the packing of the filled part per cell and faces use
  the harmonic mean (series conduction); `constant` keeps k exactly.
- T.4 (2026-11-26): with the pore gas, `equilibrium` adds the gas heat capacity
  phi_g rho_g c_g V to the grains' and carries T with the gas; `twoTemperature`
  reads `0/Tg`. The gas mass flux carried is the gas solve's (relative to the
  grains) plus the gas the grains carry, (dt phiL - dm/rhoP) at the donor density,
  with the same face rules as the gas solve (decisions 0128, 0129). Transport is
  upwind in advective form, C_P dT_P = dt c_g sum_in |F| (T_up - T_P), so a
  uniform temperature stays uniform under any flow (1e-14), sub-cycled to Courant
  0.5. The exchange is Gunn's h on a_s = 6 phi_s/d, applied exactly per cell
  (the difference decays as exp(-H (1/Cs + 1/Cg) dt), Cs Ts + Cg Tg kept), after
  the transport (first-order splitting; dt is far below the exchange time in
  practice). Gas conduction uses the gas part of Zehner-Schluender,
  (1 - sqrt(phi_s)) kGas. The Schumann test was planned at 3 % but 20 % in h moves
  theta by only 3 %, so the test asks 1 % (got 0.38 %).
- T.5 (2026-11-26): no new code. In the void of a `gasFlow` case the gas cells
  have no grains (C_s 0, no exchange), the gas conducts with kGas and is carried
  by the gas solve's mass flux; walls hold Tg through `0/Tg` (fixedValue). Laminar
  air between plates at a fixed wall temperature gives the developed Nu 7.457
  against 7.5407 (-1.1 %, from the upwind transport on 0.3 mm cells and the
  axial conduction at Peclet 50). gasFlow is not stable at a 2e-4 s step on
  these cells (the solid sets the step, decision 0131); the test uses 4e-5 s.
