# 0125 Version 2: gas in the pores and the void (plan, F series)

Date: 2026-11-24. Status: accepted (Sia chose "Go: Stage A then B", 2026-11-24 14:31 UTC).
Plan doc: https://claude.ai/artifact/MzB8dAtSPDiU77jKdDifCi. References and values:
`docs/spec/gas-references.md`, `tools/gas_reference.py`.

## Context

Sia (2026-11-24): version 2 couples a fluid to the granular solver; plan first and
report each step. Earlier talk: solids fraction phi_s = alpha rho / rho_p, Ergun or
Di Felice drag, volume change from dilatancy or the cap so venting is realistic,
care at partly filled cells. Sia accepted the recommendations: gas only (air), no
temperature, both stages.

## Decision

1. **Own library `src/celFluid`**, called from the celFoam driver beside
   `celEulerian`; OpenFOAM's `fvMatrix` only for the implicit solve (decision 0001).
2. **Fractions from v1 fields**: phi_s = alpha rho / rho_p, phi_g = 1 - phi_s. No new
   transported fraction; gas mass per cell phi_g rho_g is remapped with the solid.
3. **Gas**: ideal, isothermal at the case temperature (default 293.15 K), Sutherland
   viscosity.
4. **Stage A (F.1-F.5), pore pressure**: gas mass balance in the bulk with Darcy-Ergun
   flux relative to the grains, solved implicitly each step for p:
   d(phi_g rho_g)/dt + div(phi_g rho_g u_s) + div(rho_g q) = 0, q = -k/mu grad p
   (Ergun's inertial term as a p-dependent correction to k). Linearised about p0 this is
   a diffusion equation with c = k p0/(mu e) for a rigid bed and
   c = k/(mu (e/p0 + 1/M)) for a skeleton with constrained modulus M; those closed
   forms are the F.2 and F.3 tests. Void cells hold the boundary pressure
   (atmosphere); partly filled cells take k from the filled part; free surface and
   outlet faces carry the void pressure.
5. **Coupling to the solid**: effective stress. The constitutive models get
   sigma' = sigma + p I (tension positive); the momentum equation sees the total
   stress sigma' - p I. The solid volume change enters the gas balance through u_s.
6. **Stage B (F.6-F.8), gas momentum**: a gas velocity in void and pores, low-Mach
   with an implicit pressure (step set by the solid CFL, not by the 340 m/s sound
   speed), drag treated implicitly. Stage A must be reproduced within 1 %.
7. **Drag laws, selectable** (Sia, 2026-11-24: "I want multiple drag laws to be
   available to the user, but for testing, we should use ergun-wen-yu"): a
   run-time-selected `dragModel` in `celFluid` (`drag ergunWenYu;` in
   `constant/gasProperties`). F.2 adds the selector with `ergunWenYu` (Gidaspow:
   Ergun below gas fraction 0.8, Wen and Yu above), the default; F.7 adds `ergun`,
   `wenYu`, `diFelice`, `syamlalOBrien` and `beetstra`, each checked against its
   closed form. Every validation test (F.2 to F.9) runs with `ergunWenYu`. In Stage A
   the pore flux uses the same law (k = phi_g^2 / K).
8. **Parallel and fast** (Sia, 2026-11-24): every gas field and solve runs under
   MPI. Every F test from F.2 on also runs on 2 and 4 processors and must match serial
   to 1e-10 (pressure solve converged to 1e-12) with gas mass conserved to 1e-12. F.9
   adds a scaling test: gas on at most 1.5 times the cost of gas off per step, and at
   least 2.5 times faster on 4 processors than serial.
9. **Tests**: each package one closed form or published experiment; with gas off,
   every v1 result is unchanged bit for bit.

## Consequences

- One or two more fields per cell (p, later U_g), more restart data, one linear solve
  per step (MPI through OpenFOAM).
- Dilation (A.4) and the cap (K series) now change discharge through suction; the
  3.6g dilatancy option may be needed for fine powders.
- The Couto wheat (4 mm) relaxes pore pressure in about 1 ms, so it is a control case:
  gas on should leave its rate nearly unchanged.
- Found while reading (F.0): the Crewdson et al. (1977) record has no abstract online
  and the Roche et al. (2010) page refuses automated reading; F.5 and F.8 read them
  first, with Muite et al. (2004) and Rondon et al. (2011) as the fall-backs.
