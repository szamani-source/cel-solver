# 0126 Gas fractions and the celFluid skeleton (WP F.1)

Date: 2026-11-24. Status: accepted (follows decision 0125).

## Context

Version 2 adds gas in the pores (decision 0125). F.1 sets up the library and the
fractions every later package uses, with no pressure solve yet. Sia's idea for the
solids fraction is phi_s = alpha rho_bulk / rho_particle. Sia also wants several drag
laws for users, with Ergun-Wen-Yu for every test (decision 0125 item 7).

## Decision

1. **Library `src/celFluid`, class `poreGas`**, built by Allwmake after celCoupling and
   linked by celFoam. `poreGas::New(solid)` returns null when `constant/gasProperties`
   is absent or says `model none`; celFoam then does exactly what v1 did (the only
   extra work is looking for the file). With `model poreGas`, celFoam calls
   `correctFractions()` after the remap and after moving walls, every step.
2. **Fractions from v1 fields.** eulerianSolid's `rho` is the cell density (mass over
   the open cell volume), which is alpha times the material density rho/alpha. So
   Sia's phi_s = alpha rho_m / rho_p is `rho / rhoP` in code, and phi_g = 1 - phi_s.
   Both are fractions of the open volume. Cells inside a body (open fraction below
   1e-6) get phi_s = phi_g = 0.
3. **Permeability of the filled part** (decision 0125 item 4): e_m = 1 - rho/(alpha
   rhoP) and k = mu e_m^2 / K(e_m, slip 0) of the Ergun-Wen-Yu law: Ergun's
   e^3 d^2 / (150 (1-e)^2) below e_m = 0.8, Wen and Yu's Stokes limit
   d^2 e^4.65 / (18 (1-e)) above. Void cells (alpha below 1e-6) and cells packed past
   rhoP get k = 0; the void will hold the boundary pressure (F.4).
4. **Drag laws.** `drag` is read from gasProperties (default `ergunWenYu`); any other
   name stops with an error until F.2 adds the run-time `dragModel` selector and F.7
   the other laws. The zero-slip k above is what that selector will give for
   `ergunWenYu`, so F.2 replaces the function without changing F.1's results.
5. **Fields** `phiS`, `phiG`, `k` are written with the case's fields and never read: a
   restart recomputes them from rho and alpha.
6. **gasProperties** keys: `model` (none or poreGas, default poreGas), `rhoP` [kg/m^3],
   `d` [m], `T` [K, 293.15], `p0` [Pa, 101325], `drag` [ergunWenYu]. Air is ideal
   (R 287.05) with Sutherland viscosity, as `tools/gas_reference.py`. Documented in
   the manual (`docs/manual/case-setup.html#gas`, which `tests/test_manual.sh` requires
   for every keyword the code reads); the first gas tutorial is F.9's `gasHopper`.

## Consequences

- `tests/test_gas_fractions.sh` (about 80 s): settled column phi_s 0.048 % from the set
  packing; k equals `tools/gas_reference.py` to 7e-16; solid fields byte for byte the
  same with no gasProperties, `model none` and `model poreGas`; 2 processors match
  serial to 1e-13.
- The parallel check uses plain Mohr-Coulomb. With the T0 rule set the solid itself
  parts between serial and 2 processors by about 1e-5 in this column (flowImpedance and
  shear-cap switches flip on round-off, the known issue of WP 6.1, decision 0080).
  Without those rules it matches to 1e-14. F.1 does not change this; F.9 must allow
  for it when it compares parallel runs of F.5 and F.8.
