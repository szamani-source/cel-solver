# 0112 Cohesion before the Linux install: simple cases (packages C.1-C.3)

Date: 2026-11-20. Status: accepted (plan; results added as the packages close).

## Context

Sia (2026-11-20 16:04 UTC): cohesion is the second item to address before the
Linux installation, after mesh dependence (decision 0111) and before
compaction. Sia (18:19 UTC): develop on simple test cases, not the silo.
WP 3.5 (decision 0054) already has cohesion in Mohr-Coulomb and
Drucker-Prager (c, tensile strength sigmaT, else c cot(phi), at contact faces
and in void opening) and checks the silo arch against Jenike. Open:

- the nonlocal law of decision 0111 read P alone, so with c > 0 it let
  material flow at tau = mu_s P, below the cohesive yield c + mu_s P;
- no simple case checks a cohesive body standing or failing against a
  closed-form limit, and none on more than one mesh.

## Decision

Packages, each with one done-when test:

| WP | Title | Done when |
| --- | --- | --- |
| C.1 | Cohesion in the nonlocal law | The law acts on P* = P + c cot(phi) (mu = tau/P*, I with P*); `tests/test_nonlocal.sh` kernel check repeats homogeneous shear with c = 500 Pa: mu* = mu_s + b I* exactly |
| C.2 | Vertical cut (cohesive column) | The a = 0.5 column with c, free right face: the critical c at which it stands lies between Terzaghi's plane-wedge limits H_c = 2.67 to 4 (c/gamma) tan(45 + phi/2) (with and without a tension crack), on 10 / 5 / 2.5 mm cells, critical c within 10 % between 5 and 2.5 mm; MC and DP |
| C.3 | Cohesion with the nonlocal law in the sheared layer | The sheared layer of N.2 with c: band and profile against the 1D nonlocal solution with P* (within 5 %) |

## Consequences

- With `nonlocalA` and c > 0, results change (the law now respects the
  cohesive yield); cases without nonlocalA are unchanged.
- The face test of decision 0111 (2 eta dc < z_s) uses P* too.

## C.1 result (2026-11-20)

`celNonlocal` reads P* = P + c cot(phi) (c from the model's `c`): mu, g_loc,
the root solve and the face test. `tests/test_nonlocal.sh` kernel check now
also runs c = 500 Pa (MC and DP): mu* = mu_s + b I* exactly (24 runs).

## C.2 result (2026-11-20)

`tests/test_vertical_cut.sh` (about 4 minutes), kinetic-energy peak over
0.5 s (stands below 1e-5 J/m):

| Model, cells | 230 Pa | 250 Pa | 275 Pa | 300 Pa | 344 Pa |
| --- | --- | --- | --- | --- | --- |
| MC 10 mm | | stands | | stands | |
| MC 5 mm | | fails | stands | stands | |
| MC 2.5 mm | | fails | stands | stands | |
| DP (plane strain) 5 mm | fails | | | | stands |

Critical c in (250, 275) Pa on 5 and 2.5 mm, inside Terzaghi's 230 to
344 Pa, converged to 10 %; 10 mm is weaker (200 to 250 Pa, earlier scan).
The test checks (250, 300) on both meshes and reports 275 Pa.

## C.3 result (2026-11-20): the face viscosity is the law's tangent

First C.3 runs (c = 50 Pa) carried 7 to 13 % more shear than the 1D
solution: with P* the secant test 2 (P*/g) dc < z_s flagged no cell as
flowing, so every face kept the elastic jump z_s. Changed (both with and
without c): a cell flows where 2 (P*/g) dc < 10 z_s, and its face
viscosity is the law's tangent d tau/d gamma_dot = b d sqrt(rho_s P*) (as
the mu(I) flow impedance), or the secant if smaller. Margins 10 and 100
agree to 0.2 %. `tests/test_nonlocal_shear.sh` (about 6 minutes) now:

| Case | 2 d | 1 d | 0.5 d |
| --- | --- | --- | --- |
| c = 0: band (1D) | 3.83 (3.39) d, 1.4 % of V | 3.28 (3.25) d, 0.6 % | 3.23 (3.21) d, 0.2 % |
| c = 50 Pa: band (1D) | | 4.17 (4.15) d, 0.3 % | 4.08 (4.11) d, 0.4 % |

The shear stress equals the 1D value to 0.2 % (it was 4 to 9 % high with
the secant viscosity, decision 0111 N.2 table).
