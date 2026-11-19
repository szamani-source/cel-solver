# Conical hopper validation case (Sia's validation case, 2026-11-19)

Sia supplied two papers on one conical-hopper silo (2026-11-19) and asked us to run the
same case on the axisymmetric wedge (decision 0061) with Mohr-Coulomb and with
Drucker-Prager, using the papers' material values. This page holds everything a
package needs to set up and judge the runs. The numbers are in `hopper-data/`; the
choices behind them are in decision 0100. Packages H.1 to H.6 in `work-packages.md`.

## The papers

| # | Paper | What it gives us |
| --- | --- | --- |
| P1 | Q.J. Zheng, A.B. Yu, "Finite element investigation of the flow and stress patterns in conical hopper during discharge", Chem. Eng. Sci. 129 (2015) 49-57, doi 10.1016/j.ces.2015.02.022 | Abaqus 6.10 Eulerian (CEL) runs, Mohr-Coulomb, quarter model. Basic hopper: outlet-centre velocity over time (Fig. 3), mass/funnel flow (Figs. 6-9), flat-bottom rate against D0 (Fig. 10, Beverloo C = 0.59), wall pressure in the Couto silo at 0.5 s (Fig. 11) |
| P2 | A. Couto, A. Ruiz, P.J. Aguado, "Experimental study of the pressures exerted by wheat stored in slender cylindrical silos, varying the flow rate of material during discharge. Comparison with Eurocode 1 part 4", Powder Technol. 237 (2013) 450-467, doi 10.1016/j.powtec.2012.12.030 | **Experiments**: wheat in a 1 m steel silo with a 34.3 deg hopper; wall pressures at 8 load cells after filling, at rest and during discharge (Table 1, seven tests) |

Copies: Sia's laptop, `CEL_ref/` (not in the repository; publisher copyright).

## Case A: the Couto test silo (experiment, wall pressure)

Geometry (P2 Fig. 3), all in m, z measured up from the outlet plane:

| Item | Value |
| --- | --- |
| Cylinder diameter d_c | 1.000 |
| Cylinder height h_c | 2.000 (z = 0.48 to 2.48) |
| Hopper half angle from vertical beta | 34.3 deg |
| Hopper height h_t | 0.480 (apex 0.738 below the transition) |
| Outlet diameter d_h | 0.350, central |
| Fill | level to the top, z = 2.48 (fill weight 14.75 to 14.89 kN; geometry gives 1.756 m^3 x 8.397 kN/m^3 = 14.75 kN) |
| Load cells | GA1 to GA7 on the cylinder, 250 mm apart, GA1 at z = 0.575; GA0 on the hopper wall at z = 0.41 (`hopper-data/couto2013_cells.csv`) |

Material (wheat, P2 Section 2.3; P1 Table 1 and Section 3.4):

| Item | Value | Source |
| --- | --- | --- |
| Bulk density rho | 856 kg/m^3 (8.397 kN/m^3) | P2, P1 |
| Young's modulus E | 10.674 MPa | P2 (P1 basic case 10 MPa) |
| Poisson's ratio nu | 0.3 | P1 |
| Friction angle phi | 30.17 deg | P2 |
| Cohesion c | 0 | P1 (P2: slight moisture cohesion, not measured) |
| Dilation angle psi | 1 deg (Sia: 0.5 to 1); 5 deg as P1's sensitivity run | Sia; P1 |
| Wall friction mu | 0.20 (measured per test 0.17 to 0.28, Table 1) | P2 |
| Grain size d | about 4 mm (equivalent sphere of a wheat kernel, roughly 6.5 x 3.3 x 3.0 mm) | typical published kernel sizes; only for the cell size and mu(I) if used |

Reference results:

1. **At rest after filling** (P2 Table 1, static max/min rows, seven tests): cylinder
   wall pressure rises from 1.4 to 1.7 kPa at GA7 to 4.0 to 4.5 kPa at GA1, hopper
   wall GA0 8.6 to 9.6 kPa; vertical stress at the transition p_vt 11.7 to 13.4 kPa,
   K 0.31 to 0.37.
2. **During discharge** (Table 1, discharge max rows): GA0 18.7 to 24.8 kPa (the
   switch peak), GA1 7.1 to 8.2 kPa, GA2 5.8 to 6.5 kPa, upper cells close to the
   static values. The experiment's rate is set by a screw conveyor (0.6 to 4.9 kg/s);
   free discharge through 0.35 m would be about 110 kg/s. The peak barely depends on
   the rate (P2), so P1 compares it with free discharge at 0.5 s, and so do we.
3. **P1's FEM at 0.5 s** (Fig. 11, `hopper-data/zheng2015_fig11_fem.csv`): peak 27.5
   kPa at the transition; Nedderman's closed form (P1 Eqs. 9-10: K_A 0.338, Q0 13.0
   kPa, C_p 2.23, m 3.73) gives 4.4 kPa at the bottom of the cylinder and 29 kPa at the
   top of the hopper.

## Case B: P1's basic hopper (discharge)

| Item | Value |
| --- | --- |
| Silo diameter D | 1.0 m |
| Total height H | 2.48 m |
| Hopper half angle alpha | 30 deg (hopper height 0.779 m from the outlet edge, so the cylinder is 1.701 m; P1 gives only D, H, alpha, D0) |
| Outlet D0 | 0.10 m |
| Fill | to the top (assumed; P1 says "up to an expected level") |
| Material | rho 856, E 10 MPa, nu 0.3, phi 30 deg, c 0, psi 5 deg in P1 (we run 1 deg and 5 deg), wall mu 0.3 |

Reference results:

1. Vertical velocity at the outlet centre (P1 Fig. 3, `hopper-data/zheng2015_fig3.csv`):
   rises to a plateau of 1.37 m/s by about 0.45 s.
2. Flow mode: alpha 30 deg, mu 0.3, phi 30 deg is **mass flow** (P1 Fig. 6, the middle
   inset; Fig. 7 gives MFI > 0.3, the Jenike line). Mass flow index MFI = mean speed at
   the wall / mean speed on the axis.
3. Mass rate: P1 gives none for the hopper. We record ours against Beverloo with a
   hopper factor (Rose-Tanaka), as in WP 4.5, without a pass band.

## Case C: flat-bottomed silo (P1 Fig. 10)

D = 1 m, alpha = 90 deg, D0 from 0.05 to 0.5 m, same material as case B. P1's points
(`hopper-data/zheng2015_fig10.csv`) follow Beverloo with C = 0.59, z = 0:
W = 0.59 rho sqrt(g) D0^2.5 (5.0 kg/s at 0.1 m, 78 kg/s at 0.3 m).

## Models and run settings

- Both Mohr-Coulomb and Drucker-Prager, rate independent (no mu(I)), as in P1. P1
  uses Abaqus' Mohr-Coulomb with a hyperbolic flow potential (eccentricity 0.1); ours
  uses the plain cone potential with psi; the difference is at low pressure only.
- Drucker-Prager `fit compression` (triaxial compression, the stress state on the
  axis of a hopper) as the main run; `plane_strain` as a sensitivity run only.
- Wall friction mu on the silo walls (Coulomb, no adhesion).
- Cell size: Sia asked that cells resolve the shear band, 10 to 20 grain diameters
  (40 to 80 mm for wheat). With 3 to 4 cells across the band, cells are 10 to 20 mm.
  Base runs at 10 mm (50 cells across the radius, 5 across case B's outlet
  half-width; decision 0062 needs 2 on a wedge); 5 mm as the mesh check.
- Keep runs lean: case A to 0.5 s of discharge, case B to about 1 s (the plateau).

## Pass bands (set before the runs; decision 0100)

| Check | Band |
| --- | --- |
| A, at rest: cylinder cells GA1 to GA3 | within the experiment's min-max over seven tests widened by 15 % |
| A, at rest: cylinder cells GA4 to GA7 | within 15 % of the experiment range or of Janssen with K_A 0.338, mu 0.2 (the cells read 11 to 31 % above Janssen, and P1's FEM sits on Janssen there; see note) |
| A, at rest: GA0 | 8.6 to 9.6 kPa widened by 15 % |
| A, discharge 0.5 s: GA0 peak | 18.7 to 24.8 kPa (experiment) widened by 15 %; also report against P1's 27.5 kPa |
| A, discharge: GA1, GA2 | experiment range widened by 15 % |
| B: outlet-centre velocity plateau | 1.37 m/s within 20 % (Sia: about 20 % accepted for flow rates) |
| B: flow mode | mass flow (MFI > 0.3) |
| C: Beverloo C | 0.59 within 20 % at D0 = 0.1, 0.2, 0.3 m |

Note on the upper cylinder: P2's static pressures at GA4 to GA7 are 11 to 31 % above
Janssen (and P1's FEM), while GA1 to GA3 match it. Filling stresses are the likely
cause; a rate-independent continuum filled in one step cannot be expected to show them,
so either reference is accepted there. `tests/test_hopper_spec.sh` prints the gaps.

Both models must pass on their own; MC against DP differences are reported, not banded.
