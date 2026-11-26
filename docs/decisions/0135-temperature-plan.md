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
