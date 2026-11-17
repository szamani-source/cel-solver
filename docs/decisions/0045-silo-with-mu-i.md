# 0045 Silo discharge with mu(I) and the flow impedance (WP 3.6d)

Date: 2026-11-17. Status: accepted (WP 3.6d closed, done-when not met; finding recorded).

## Context

WP 3.6d runs Sia's target T0 (tank W 0.4 m, fill 0.6 m, outlet D 0.06 m, walls and floor
20 deg, settle 0.3 s then open) with the mu(I) law (decision 0037 placeholders) and
`flowImpedance yes` (decisions 0040, 0044), which made the inclined layer independent of
E. Done-when: `tests/manual/test_discharge_stiffness.sh` (rate at 1 and 10 MPa within
10 % of 0.1 MPa) at 10 mm and 5 mm, C within decision 0019.

## Findings

1. **Crash at the first step** of a stress-free fill: `flowViscosity` took
   sqrt(rhoS/P) with P about 1e-320 Pa (floating-point trap). Fixed with the guard
   `FrictionLaw::inertialNumber` already has (no viscosity unless P > 1e-300 rhoS).
2. **The nearly-full snap of decision 0044 broke the stiff silo.** Above the outlet the
   converging flow opens about 1e-6 of void per step; with a 1e-4 tolerance those cells
   were set full every step, which adds material volume each time. At 10 MPa (28 000
   steps) the stream above the outlet inflated to half its density (967 kg/m3) and the
   rate fell to 12 kg/s/m. Because the snap is per step, its effect grows with the step
   count, so with E. Bisection at 1 MPa, 10 mm (each flow-impedance part switched off in
   turn): only the snap mattered (49.8 -> 61.6 kg/s/m without it). Without the snap,
   dropping the wall ghost, the plastic predictor or the phi_m cap changed the rate by
   under 6 %; dropping the flow shear jump made it far slower (18.9).
3. **Fix (behind `flowImpedance yes`):** a cell is snapped full only while its density
   stays within 1e-3 of the density its stress implies, rho0 (1 + p/K). Noise void in a
   dense cell is still removed; material that is loosening is not. Tried and dropped:
   snapping only next to a free surface (the inclined layer failed: 11 % and 86 % profile
   errors), a bound at rho0 alone (1 MPa unchanged).
4. **Run-to-run scatter.** The same 10 MPa case through celFoam (the tutorial; its time
   steps shorten to hit write times) follows the test run to 0.3 % until t = 0.9 s, then
   drifts: 57.5 against 67.7 kg/s/m. Plastic flow parts by round-off growth (as the
   collapse, decision 0032), so a single 10 mm rate is uncertain by about 15 %.

## Results (rate in kg/s per m depth, 97 to 88 % window; C = rate / (rho sqrt(g) D^1.5))

| Run | 0.1 MPa | 1 MPa | 10 MPa |
| --- | --- | --- | --- |
| Mohr-Coulomb only, 10 mm (decision 0034) | 74.1 | 65.7 (-11 %) | 70.9 (-4 %) |
| mu(I), no flow impedance, 10 mm | 70.8 | | 68.8 (-3 %) |
| mu(I) + flow impedance as in decision 0044, 10 mm | 62.3 | 49.8 (-20 %) | 12.4 (-80 %) |
| no snap at all, 10 mm | 62.1 | 61.6 (-1 %) | 63.2 (+2 %) |
| **this package, 10 mm** | **65.5 (C 0.79)** | **53.8 (-18 %)** | **67.7 (+3 %)** |
| **this package, 5 mm, 2 processors** | **77.6 (C 0.94)** | **61.2 (-21 %)** | not run (about 70 min) |

`tests/manual/test_inclined_layer.sh` still passes with the fix (profile error 4.1 / 0.9
/ 0.8 % at 0.1 / 1 / 10 MPa; was 3.6 / 0.9 / 0.8).

## Decision

1. WP 3.6d closes with its done-when not met: 1 MPa is about 20 % slow on both meshes
   (10 MPa is within 3 % at 10 mm). Every snap variant gave 49 to 54 at 1 MPa, 10 mm;
   only no snap at all gave 61.6, but the layer needs the snap.
2. C is 0.69 to 0.94, below the decision 0019 band (1.0 to 1.9). Expected in part: the
   continuum references (Staron 2014, Zhou 2017) use mu_s 0.32 to 0.38; T0 has 0.70, and
   a rougher material discharges slower. Rising with refinement (0.79 -> 0.94 at
   0.1 MPa).
3. The T0 tank is now a celFoam tutorial, `tutorials/siloDischarge` (Allrun, plot.py),
   the first case for Sia to look at. It needs no test driver: celFoam reads 0/alpha and
   0/rho from setFields. celFoam's start-up print of the wave speed divided by the
   smallest density (0 with void); it now uses rho0.
4. Next (a later package): why 1 MPa is slow (the snap at 1 MPa, or the normal Riemann
   jump, decision 0028), measured as a scatter band over a few perturbed runs, not one run.

## Consequences

- Runs with `flowImpedance yes` and a stress-free start no longer crash, and the stiff
  silo is no longer throttled by the snap. Runs without `flowImpedance` are unchanged.
- `tests/manual/test_discharge_stiffness.sh` now runs with mu(I) (`--mui`);
  `CEL_NO_MUI=1` gives the old plain Mohr-Coulomb sweep.
- Rate differences under about 15 % at 10 mm are within run-to-run scatter.
