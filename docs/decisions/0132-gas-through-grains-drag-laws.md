# 0132 Gas through and around the grains; more drag laws (WP F.7)

Date: 2026-11-24. Status: accepted (follows decisions 0125 items 6-9, 0127, 0131).

## Context

F.7 adds the gas momentum inside the material (F.6 solved it only in the void; the pores
kept Darcy-Ergun), the drag on the grains, the user drag laws beside the default
(decision 0125 item 7) and validation pack test 2, the onset of fluidisation. Sia's
rules: ergunWenYu in every test; 2 and 4 processors match serial to 1e-10 with gas mass
to 1e-12; report the cost with gas on against off (F.9 target 1.5x).

## Decision

1. **Five more drag laws**, run-time selected like `ergunWenYu` (`drag <name>;` in
   `constant/gasProperties`): `ergun` (Ergun at every gas fraction), `wenYu` (Wen and Yu
   at every gas fraction), `diFelice` (one Dallavalle sphere at the superficial slip
   times e^-beta, Re = e rho u d/mu), `syamlalOBrien` (MFiX form, Re = rho u d/mu) and
   `beetstra` (Beetstra, van der Hoef and Kuipers 2007, K = 18 mu e (1-e) F/d^2,
   Re = e rho u d/mu). Each is written in a Cd Re form that is finite at zero slip. The
   forms and the conventions are in each class's header and in
   `tools/gas_reference.py`, which `tests/test_drag_laws.sh` checks them against
   (celDragTable prints K). Sources: El Geitani et al. (2022) eq. 30 for Di Felice (they
   write Re without e; Di Felice's own definition uses the superficial slip, kept here);
   Syamlal-O'Brien and Beetstra in their standard published forms (papers not read; the
   El Geitani paper lists but does not give them). Beetstra's packed limit is
   Carman-Kozeny and its dilute limit one Stokes sphere (both checked). The generic slip
   root finder now widens its bracket: Syamlal-O'Brien's K falls with the slip, and the
   old bracket returned the Darcy speed (an 18 % error in a permeameter).
2. **Pore gas momentum** (`model gasFlow`, `poreInertia yes` by default): in cells with
   grains the gas obeys e rho_g (w - w_o)/dt = -e grad p - K w (w the slip, q = e w the
   superficial flux relative to the grains), the drag implicit. That is the Stage A flux
   law with K + e rho_g/dt for K and the gradient less rho_g q_o/(e dt): the same Newton
   solve, a wrapped drag law (`inertialGasDrag`) and an explicit face mass flux from the
   old relative flux (as the void's H/A flux in F.6). It applies on faces between gas
   cells, between gas and void (scaled by the gas side's share of the centre distance)
   and on fixed-pressure boundaries; inflow patches (fixed flux) keep Darcy. A face that
   was in the void flow last step starts from e (u_g - u_s).n. As dt/tau grows
   (tau = e rho_g/K, 7 microseconds for 0.2 mm sand) it is Stage A; in coarse or dilute
   beds it caps the mobility at e dt/rho_g instead of e^2/K.
3. **Drag on the grains** stays the mixture force of F.3: each cell gets -sum(p_f S_f),
   which is the buoyancy -phi_s grad p plus the drag K w = -e grad p - e rho_g Dw/Dt
   less the pore gas inertia. The reaction of that inertia, -V (rho_g/dt)(q - q_o), is
   built (`inertiaForce yes`) but off by default: applied explicitly one step late it
   fed back on the lifting bed (the relative flux falls as the grains rise) and the
   onset run blew up within 0.05 s of lift; it is a 1e-3 effect (rho_g against
   rho_s), the order of the gas weight and the grain acceleration term already left out
   (decision 0131 item 5). The force on the grains therefore lacks the gas's own
   inertia; F.8/F.9 can revisit it with the solid's acceleration in the gas equation.
4. **Inflow as a Function1 of time** (`inflow { bottom table ((0 0) (0.28 0.048)); }`),
   for ramps; a number still works.
5. **Onset and settling cases**: sand d 0.2 mm packed at e_mf 0.38 (Wen and Yu's
   constants 1/(phi e^3) = 14 and (1-e)/(phi^2 e^3) = 11 stand for packings near 0.38;
   U_mf of an Ergun-type law moves about 9 % per 0.01 of e_mf), plain Mohr-Coulomb (the
   T0 rules part serial from parallel, decision 0130), air ramped through the floor
   over 0.28 s. U_mf is where the floor's pressure drop (u + u^2 fit before lift)
   reaches the bed weight. Settling: a slab at e 0.9 (50 micrometre grains) falling in a
   closed column.

## Results

`tests/test_drag_laws.sh` (about 5 s): each law equals its closed form to 2e-16 over
gas fraction 0.35 to 0.999 and Re 1e-3 to 1e3 and at zero slip (limit 1e-10);
ergunWenYu equals ergun below 0.8 and wenYu above exactly; dilute and packed limits as
above; a permeameter with each law within 1.3e-7 of the law integrated up the column.

`tests/test_fluidisation.sh` (about 5 min):

- Pore gas momentum, pressure step across a frozen 5 mm bed: the flux follows
  e rho_g dq/dt = e^2 G - K q within 1.1 % from 0.4 tau (backward Euler's lag at
  dt/tau = 0.025), 0.1 % at 4 tau.
- Onset: U_mf 0.03372 m/s, -2.2 % of Wen and Yu (0.03446), +0.5 % of the law's own
  U_mf at e 0.38 (0.03355). The bed rises from 0.0335 m/s (extrapolated). The pressure
  drop after lift is the weight +3.6 % (the emptying floor cell). Without the pore gas
  momentum (F.6) the onset is the same (+0.5 %) but the pressure profile bends in the
  top three cells (the void cell above holds 4.4 Pa); with it the gradient is uniform
  to the surface.
- Settling: 0.1155 m/s, +7.6 % of Richardson and Zaki (0.1073; Dallavalle U_t, Rowe's
  n), +1.9 % of the drag law's own hindered speed (Wen-Yu gives e^4.65 in the Stokes
  limit; at Re_t 0.57 it is 5.6 % above Richardson-Zaki).
- Stage A: the venting column and the outlet box with the pore gas momentum within
  0.04 % and 0.1 % of model poreGas; `model poreGas` runs are untouched (the momentum is
  gasFlow only); `tests/test_gas_flow.sh` still passes (0.08 %, 0.12 %).
- Parallel, onset run 2 x 1 and 1 x 4 (cuts through the bed and the free surface): up to
  the onset p, Ug and U agree to 7e-11 of their largest values (per write time up to
  1e-9 for U while the bed is nearly still, at 1e-5 m/s), gas mass 1e-15. After lift
  1 x 4 parts by up to 1e-1 (2 x 1 stayed at 3e-11 in this run): the lifting bed sits
  at the Mohr-Coulomb apex everywhere and amplifies round-off. With other write
  intervals (2.5 ms) a run parted by 1e-3 before lift (0.0175 s); from a restart at
  0.015 s it did not, cause not found. Without the pore gas momentum the same case parts
  by 1e-6 before lift; linear elastic grains part too. F.9.
- Cost per step, gas on against off: outlet box (480 cells) 2.5x with the pore gas
  momentum against 2.0x poreGas (F.6's gasFlow 2.3x); fluidised bed (160 cells, the dry
  bed only sits) 3.8x. The inertial law keeps Ergun's closed-form slip (it is affine in
  u); on the Wen-Yu branch and the other laws every face calls the root finder. F.9.

Validation pack test 2: `tests/validation/gas/fluidOnset` (see its results.md).

## Consequences

- gasFlow cases now carry the pore gas momentum; `poreInertia no` gives F.6.
- Not restarted yet: the old face fluxes (a restart starts the pore gas from rest
  relative to the grains), F.9.
- The lifted bed and settling slab are coarse: plug lift, no bubbles; F.8 shows a
  collapse.
