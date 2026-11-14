# Silo reference values (WP 3.0)

Status: written 2026-11-14. The numbers the silo packages (3.1p, 3.2p, 3.3p, then the
C++ 3.1 to 3.4 and 3.5) are judged against. Each value names its source. "Opened" means
the source was fetched and read on that date; anything not found in an opened source is
marked **not verified**. "Derived" means worked out here and checked by
`tools/silo_reference.py --check`.

Pressures in this page are **compression positive** (the solver convention is tension
positive). T0 inputs are from `docs/spec/target-problems.md`: W = 0.40 m, H = 0.60 m,
rho = 1800 kg/m3, phi = 35 deg, c = 0, nu = 0.3, wall friction delta_w = 20 deg
(mu_w = tan 20 = 0.364), outlet D = 0.04, 0.06, 0.08 m.

Regenerate every table in section 5 with `python3 tools/silo_reference.py` (and
`--H 1.2` for the tall fill). `tests/test_silo_reference.sh` fails if this page and the
script disagree.

## 1. Sources

| # | Reference | Link | Opened | Used for |
| --- | --- | --- | --- | --- |
| S1 | M. Sperl, "Experiments on corn pressure in silo cells: translation and comment of Janssen's paper from 1895", Granular Matter 8, 59 (2006) | https://arxiv.org/abs/cond-mat/0511618 | 2026-11-14 | Janssen's equation, his measured K and wall friction, the assumptions |
| S2 | L. Staron, P.-Y. Lagree, S. Popinet, "The granular silo as a continuum plastic flow: the hour-glass vs the clepsydra" (2012; journal not shown in the opened text) | https://arxiv.org/abs/1211.5916 | 2026-11-14 | Continuum 2D discharge constant, fill-height independence, constant friction vs mu(I) |
| S3 | L. Staron, P.-Y. Lagree, S. Popinet, "Continuum simulation of the discharge of the granular silo: a validation test for the mu(I) visco-plastic flow law" (arXiv 2013; journal not shown in the opened text) | https://arxiv.org/abs/1306.6279 | 2026-11-14 (sections I to III only) | Continuum vs discrete (contact dynamics) constants, silo width effect |
| S4 | S. Dunatunga, K. Kamrin, "Continuum modelling and simulation of granular flows through their many phases", J. Fluid Mech. (2015) | https://arxiv.org/abs/1411.5447 | 2026-11-14 | MPM plane-strain silo constant; continuum column-collapse run-out; table of published collapse laws |
| S5 | A. Janda, I. Zuriguel, D. Maza, "Flow rate of particles through apertures obtained from self-similar density and velocity profiles", Phys. Rev. Lett. 108, 248001 (2012) | https://www.unav.edu/documents/14434042/14611288/2012_PRL108_248001_Janda.pdf/efa70593-a794-c11c-78ad-5982051a42a4 | 2026-11-14 | 2D experiments: D^(3/2) law, outlet velocity sqrt(2 g R), packing at the outlet |
| S6 | C. Mankoc et al., "The flow rate of granular materials through an orifice", Granular Matter 9, 407 (2007) | https://arxiv.org/abs/0707.4550 | 2026-11-14 | Beverloo form and 3D constant range, fill-height and width conditions |
| S7 | A. Bhateja, "Connection between the mass flow rate and the base and bulk normal stresses in silo discharge" (arXiv 2017) | https://arxiv.org/abs/1706.09582 | 2026-11-14 | 2D DEM constant with bulk density; quotes the continuum constants of S2/S3 |
| S8 | E. Lajeunesse, J. B. Monnier, G. M. Homsy, "Granular slumping on a horizontal surface", Phys. Fluids 17, 103302 (2005) | https://www.ipgp.fr/~lajeunes/download/2005-pof-slumping.pdf | 2026-11-14 (prefactors not readable in the extracted text) | Column collapse exponents and regime limits; failure angle |
| S9 | L. Staron, E. J. Hinch, "Study of the collapse of granular columns using DEM numerical simulation" (2005) | https://arxiv.org/abs/physics/0501022 | 2026-11-14 | Quotes Lube et al.'s quasi-2D run-out law with prefactors |
| S10 | Wikipedia, "Lateral earth pressure" | https://en.wikipedia.org/wiki/Lateral_earth_pressure | 2026-11-14 (formulas did not survive extraction) | Jaky's references only |

Not opened (cited by the sources above, listed so a later session can fetch them):
Beverloo, Leniger, van de Velde (1961); Lube, Huppert, Sparks, Freundt, Phys. Rev. E 72,
041301 (2005); Lagree, Staron, Popinet, J. Fluid Mech. 686, 378 (2011); Jaky (1944,
1948); Nedderman, *Statics and kinematics of granular materials* (1992).

## 2. Settled state: Janssen in a plane (2D) bin

### Equation

Janssen (S1) balances a horizontal slice of a square cell of side s:
`dp/dx = gamma - (4 K / s) p`, so `p = (s gamma / 4K)(1 - exp(-4 K x / s))` (S1, eqs. 2
and 3). His K includes the wall friction: `K_Janssen = f p_s / p`.

For our plane bin (two walls, width W, unit depth) the wall length per unit area is 2/W
instead of 4/s. With K = sigma_h / sigma_v (the stress ratio alone, wall friction
mu_w = tan delta_w kept separate) (derived):

    d(sigma_v)/dz = rho g - (2 K mu_w / W) sigma_v
    sigma_v(z)    = rho g lambda (1 - exp(-z / lambda)),   lambda = W / (2 K mu_w)
    sigma_h,wall  = K sigma_v,   tau_wall = mu_w K sigma_v
    saturation sigma_v,sat = rho g lambda;  90 percent of it at z90 = lambda ln 10

sigma_v is the mean over the width. Total vertical balance: weight = base force + wall
friction on both walls (derived; checked to 1e-12 by the script).

### Assumptions (S1 commentary)

Janssen assumes (a) wall friction fully mobilised, (b) a constant stress ratio,
(c) stress uniform across a slice. Sperl notes (a) is not true in general but holds well
when the wall has slipped relative to the fill; (c) contradicts symmetry (it implies
shear at the centreline). Refined continuum theories change the prefactors, not the
form. Janssen's own data for wheat in wooden cells: wall friction f = 0.302 to 0.346, his
K (which includes f) 0.20 to 0.235, so the stress ratio alone was about 0.7 (S1).

Consequence for the tests: a simulation fills a tank in a way the theory does not
describe (here: gravity switched on in a full tank). The wall friction is only fully
mobilised where the material has slid down the wall. Expect Janssen to hold best away
from the top and bottom, and expect the fitted K to depend on how the fill settled.

### Which K

| K | Formula | Value for T0 | When it applies | Source |
| --- | --- | --- | --- | --- |
| Elastic, uniaxial strain | nu / (1 - nu) | 0.4286 | Material loaded under gravity between smooth rigid walls with no lateral strain and no yield. Also the plane-strain out-of-plane stress is nu (sigma_h + sigma_v) = 0.4286 sigma_v, so the state is axisymmetric about vertical and far from yield | Derived (linear elasticity) |
| Jaky, at rest | 1 - sin phi | 0.4264 | Normally consolidated sand or clay in the field | Formula **not verified** (S10 names Jaky and his 1944 and 1948 papers, but its formulas were not readable) |
| Wall yield, active | (1 + sin^2 phi - 2 sqrt(sin^2 phi - mu_w^2 cos^2 phi)) / (4 mu_w^2 + cos^2 phi) | 0.2906 | Material at yield next to a slipping wall, major principal stress near vertical | Derived (Mohr circle through the wall point; checked against a direct construction; equals Rankine active when mu_w = 0). That this is the classical "Walker" K is **not verified** |
| Rankine active | (1 - sin phi)/(1 + sin phi) | 0.2710 | Smooth wall, material at active failure | Derived |

For cohesionless Mohr-Coulomb with phi = 35 deg, any ratio between 0.271 and 3.69
(Rankine passive) is admissible without yield. The elastic value 0.4286 is inside that
range, so the smooth-wall settled state of WP 3.1p is elastic, and the expected ratio
is **0.4286** (nu/(1 - nu)), not Jaky (they differ by 0.5 percent here, by
coincidence). With frictional walls the ratio near the walls lies between the wall-yield
value 0.29 and the elastic value 0.43; which one depends on how much the material slid.

### Numbers for T0

Section 5, tables J1 to J4. Highlights:

- T0 fill (H = 0.60 m = 1.5 W): lambda = 1.28 to 1.89 m, so the bed is shallow in
  Janssen terms. Mean vertical stress at the bottom is 80 to 86 percent of hydrostatic
  (8.5 to 9.1 kPa against rho g H = 10.6 kPa). Saturation is not visible at this depth.
- Tall fill for WP 3.1p (H = 1.2 m = 3 W): bottom stress 65 to 75 percent of
  hydrostatic (13.8 to 15.7 kPa against 21.2 kPa). The spread between the two K values
  (14 percent) is larger than the 10 percent test tolerance, so the fit does determine K.
- 90 percent of saturation needs z90 = 3.0 m (K = 0.43) to 4.4 m (K = 0.29), that is 7
  to 11 tank widths. Seeing true saturation in a test needs H of about 8 W.

## 3. Discharge rate: Beverloo/Hagen law in 2D

### Form

Beverloo's law with grain size d (S5, S6; S3 eq. 4):
`Q = C sqrt(g) (D - k d)^(N - 1/2)` in N dimensions, so `D^(3/2)` in 2D. As a mass rate
per unit depth with bulk density rho_b (S7 eq. 2): `Wdot = C rho_b sqrt(g) (D - k d)^(3/2)`.

A continuum has no grain size, so k = 0 (S2: "no volume of exclusion"; S4: k = 0 "as
expected for a local relation"). For T0, with rho = 1800 kg/m3 taken as the bulk density:

    Wdot = C rho sqrt(g) D^(3/2)   [kg/s per metre depth]

Physical picture (S5, 2D experiments): grains cross the outlet at about the free-fall
speed from an arch of height about R = D/2, `v(x) = sqrt(2 g R) sqrt(1 - (x/R)^2)` (S5 eq.
3; measured centre speed 1.07 sqrt(2 g R)); the packing at the outlet is lower than in
the bulk (area fraction at the centre 0.83 for large outlets against above 0.90 in the
bulk, S5). A continuum code that keeps the material at bulk density through the outlet
will therefore tend to a higher C than experiments, unless it allows dilation.

### The constant C

| C (k) | Method | Constitutive law, geometry | Source |
| --- | --- | --- | --- |
| 1.4 (k = 0); free fit k = 0.85 d | Continuum, Navier-Stokes + VOF (Gerris) | mu(I): mu_s = 0.32, Delta mu = 0.28, I0 = 0.40; flat-bottomed 2D silo, no-slip walls (free slip gave the same result); outlet 0.06 to 0.28 W; volume rate, so mass rate = rho_bulk Q | S2 sec. III, eq. 2 |
| 1.55 (k = 0.56) | Same continuum code | Same law; W = 90 d and 180 d, H = 90 d; outlet 5.6 d to 22.5 d | S3 sec. III |
| 1.22 (k = 2.17) | Discrete, contact dynamics | Rigid discs, contact friction 0.5, same geometry | S3 eq. 5 |
| 1.42 (k = 0); free exponent: C = 1.31, exponent 1.47 | Continuum, MPM, plane strain | mu(I) with a yield criterion and elasticity (E = 1 GPa, nu = 0.3, mu_s = 0.382, mu_2 = 0.644), separation allowed; silo 1 m wide, 1 m fill, rough bottom, frictionless sides; outlet 0.06 to 0.20 m; rate fitted while 900 to 1200 kg remain | S4 sec. 4.2 |
| 1.46 with bulk density (1.23 without), asymptotic in fill height; k varies with fill height | Discrete (soft discs), 2D, W = 40 d | packing 0.84 | S7 eq. 2 and 3 |
| 0.55 to 0.65 (3D, k = 1 to 2) | Experiments, literature range | 3D | S6 (quoting Beverloo) |

Reference for T0: **C = 1.45** (centre of the three continuum values 1.40, 1.42 and
1.55); continuum codes with mu(I) span 1.40 to 1.55, about 25 percent above discrete
simulations of the same geometry (1.55/1.22 = 1.27, S3). Rates in table B1 (section 5).

Important caveat for our material: every continuum value above used mu(I) with a
static friction mu_s of 0.32 to 0.38 (phi about 18 to 21 deg). T0 uses rate-independent
Mohr-Coulomb with phi = 35 deg (tan phi = 0.70). S2 reports that a lower friction gives a
faster discharge and that at mu_s = 0.10 the rate no longer follows D^(3/2); it also
reports that constant friction (Delta mu = 0) still gives the Beverloo scaling, but
gives no number. How C changes with phi for Mohr-Coulomb is **not verified**: expect a
somewhat smaller C than 1.45, not a different exponent.

### When the rate does not depend on the head

- S2: with outlet 0.19 W, the rate is "essentially constant and independent of the
  filling height" for H = 1.5 W to 4.5 W (1.8 percent higher at 2.5 W, 3.9 percent higher
  at 4.5 W, against 1.5 W), and volume against time stays linear even from H = 0.9 W.
- S6: base pressure saturates and the rate is about constant once the fill is deeper
  than about 1.2 silo diameters (3D); experiments kept the fill above about 2 diameters.
- S3: silo width 90 d or 180 d had no effect on the rate: discharge is controlled
  locally at the outlet, not by the stress at the walls.
- S6 (quoting Nedderman) and S7: the silo must be wider than about 2.5 D (and D + 30 d)
  for the width not to matter. T0: W / D = 5 to 10, satisfied.

For T0 the fill (1.5 W) passes 1 W after 0.2 m of drop, that is after 144 kg/m has left.
At C = 1.45 that is 2.2 s for D = 0.04 m, 1.2 s for D = 0.06 m and 0.8 s for D = 0.08 m.
The steady-rate window is short for the large outlet; WP 3.3p should measure the rate
over the mass window from 95 percent to 70 percent of the initial mass (after the
start-up transient), or use H = 0.9 m as well.

### Flow pattern

- Bottom corners stay static; a fast flow forms above the outlet (S2 intro and sec. III).
- In mu(I) continuum runs, a low-pressure zone sits above the outlet between two
  high-pressure "domes"; it needs friction above about 0.3 and is absent for a Newtonian
  fluid (S2 sec. III and VI).
- The continuum leaves more material behind than the discrete simulation: about 80
  percent evacuated against 95 percent (S3 sec. III). Discrete and continuum velocities
  agree in the fast region; slow creep in the stagnant zones is not captured (S3
  abstract).
- Stagnant-zone angle: **not verified** (no opened source gives a number). Estimate
  (derived, not a test value): if the residual heaps slope at phi = 35 deg from the
  outlet edges, about 36 kg/m (8 percent of the initial 432 kg/m) remains.

### Where rate-independent Mohr-Coulomb is known to differ

- Rate-independent Drucker-Prager does not give a unique steady flow over a range of
  inclination angles (S4, discussing inclined chute flow); the mu(I) rate dependence is
  what selects it.
- S2 notes that the I-dependence raises the friction near the outlet, which matters for a
  quantitative comparison with discrete silos, while constant friction already gives
  the height independence and the D^(3/2) shape.
- With no grain size, a local continuum law cannot produce k > 0 or clogging (S4). Real
  T0 material with 1 mm grains would have k d of about 2 mm, a 5 percent effect at
  D = 0.04 m.

## 4. Granular column collapse in 2D (for WP 3.2p)

Definitions: a = H_i / L_i (initial height over initial length of the column, collapsing
one way from a wall); run-out (L_f - L_i)/L_i; H_f the final height.

| Law | Range | Source |
| --- | --- | --- |
| (L_f - L_i)/L_i = 1.2 a | a below about 2.3 | Lube et al. (quasi-2D experiments), as quoted in S9 sec. III |
| (L_f - L_i)/L_i = 1.9 a^(2/3) | a above about 2.3 | same |
| Same prefactors 1.2 and 1.9, with limits a < 1.8 and a > 2.8 | | S4 table 3, attributed to Lube et al. (2005); the table text was garbled, so the limits are uncertain |
| (L_f - L_i)/L_i proportional to a, then to a^(2/3) | change at a = 3 | S8 sec. III.B (glass beads, channel 45 mm wide); prefactors **not verified** (not readable) |
| H_f / L_i proportional to a (H_f = H_i: the top survives), then to a^(1/3) | change at a = 0.7 | S8; prefactors **not verified** |
| Failure plane at 50 to 55 deg, consistent with 45 + phi/2 | | S8 |
| Continuum mu(I), Navier-Stokes: run-out 2.2 a | a < 7 | Lagree et al. (2011), as listed in S4 table 3 (definition of a not checked) |
| Continuum MPM mu(I): 2.28 a, then 3.24 a^0.66 | change at a = 2.83 | S4 table 3 (a uses the half-width of a symmetric heap) |
| DEM: 2.5 a for a below about 2 | | Staron and Hinch, as listed in S4 table 3 |

Material dependence: S8 states the exponents do not depend on the friction angle and
that a Rankine-type coefficient changes only the constant; S8 (quoting Balmforth and
Kerswell) says the constants do depend on the internal friction angle. The experiments
used glass beads (repose 22 deg, avalanche 27 deg, S8 table I). T0 material has
phi = 35 deg, so its run-out should be at or below the glass-bead fit at the same a.

Expected values for the 3.2p cases (derived from the rows above):

| a | Experimental fit (1.2 a) | Continuum mu(I) (about 2.2 a) | Top survives? |
| --- | --- | --- | --- |
| 0.5 | 0.60 | 1.1 | Yes: a < 0.7 (S8) and a failure plane at 45 + phi/2 = 62.5 deg from the toe leaves the top intact for a < tan 62.5 = 1.92 (derived) |
| 2.0 | 2.4 | 4.4 | No |

## 5. Tables (generated by `tools/silo_reference.py`)

Table J1 to J2 and B1, T0 fill (`python3 tools/silo_reference.py`):

<!-- BEGIN silo_reference.py output -->
Inputs: rho = 1800 kg/m3, g = 9.81 m/s2, W = 0.4 m, H = 0.6 m, phi = 35 deg, delta_w = 20 deg (mu_w = 0.3640), nu = 0.3

| Stress ratio K | K | lambda (m) | sigma_v,sat (Pa) | z90 (m) | mean sigma_v at H (Pa) | sigma_v(H) / rho g H |
| --- | --- | --- | --- | --- | --- | --- |
| Elastic uniaxial strain, nu/(1-nu) | 0.4286 | 1.2822 | 22640 | 2.952 | 8461 | 0.7986 |
| Jaky at rest, 1 - sin(phi) | 0.4264 | 1.2886 | 22754 | 2.967 | 8470 | 0.7995 |
| Wall yield (active, wall slipping) | 0.2906 | 1.8908 | 33389 | 4.354 | 9078 | 0.8569 |
| Rankine active, (1-sin)/(1+sin) | 0.2710 | 2.0277 | 35806 | 4.669 | 9171 | 0.8656 |

Profile with K = 0.4286 (elastic), wall normal = K sigma_v, wall shear = mu_w K sigma_v:

| z (m) | rho g z (Pa) | sigma_v (Pa) | sigma_h wall (Pa) | tau_w (Pa) |
| --- | --- | --- | --- | --- |
| 0.000 | 0 | 0 | 0 | 0 |
| 0.100 | 1766 | 1699 | 728 | 265 |
| 0.200 | 3532 | 3270 | 1401 | 510 |
| 0.300 | 5297 | 4723 | 2024 | 737 |
| 0.400 | 7063 | 6068 | 2600 | 946 |
| 0.500 | 8829 | 7311 | 3133 | 1140 |
| 0.600 | 10595 | 8461 | 3626 | 1320 |

| D (m) | D/W | Wdot, C = 1.40 (kg/s/m) | Wdot, C = 1.45 (kg/s/m) | Wdot, C = 1.55 (kg/s/m) | emptying time at C = 1.45 (s) |
| --- | --- | --- | --- | --- | --- |
| 0.040 | 0.100 | 63.1 | 65.4 | 69.9 | 6.61 |
| 0.060 | 0.150 | 116.0 | 120.1 | 128.4 | 3.60 |
| 0.080 | 0.200 | 178.6 | 185.0 | 197.7 | 2.34 |
<!-- END silo_reference.py output -->

Tables J3 to J4 and B2, tall fill for WP 3.1p (`python3 tools/silo_reference.py --H 1.2`):

<!-- BEGIN silo_reference.py output -->
Inputs: rho = 1800 kg/m3, g = 9.81 m/s2, W = 0.4 m, H = 1.2 m, phi = 35 deg, delta_w = 20 deg (mu_w = 0.3640), nu = 0.3

| Stress ratio K | K | lambda (m) | sigma_v,sat (Pa) | z90 (m) | mean sigma_v at H (Pa) | sigma_v(H) / rho g H |
| --- | --- | --- | --- | --- | --- | --- |
| Elastic uniaxial strain, nu/(1-nu) | 0.4286 | 1.2822 | 22640 | 2.952 | 13760 | 0.6494 |
| Jaky at rest, 1 - sin(phi) | 0.4264 | 1.2886 | 22754 | 2.967 | 13788 | 0.6507 |
| Wall yield (active, wall slipping) | 0.2906 | 1.8908 | 33389 | 4.354 | 15688 | 0.7404 |
| Rankine active, (1-sin)/(1+sin) | 0.2710 | 2.0277 | 35806 | 4.669 | 15993 | 0.7548 |

Profile with K = 0.4286 (elastic), wall normal = K sigma_v, wall shear = mu_w K sigma_v:

| z (m) | rho g z (Pa) | sigma_v (Pa) | sigma_h wall (Pa) | tau_w (Pa) |
| --- | --- | --- | --- | --- |
| 0.000 | 0 | 0 | 0 | 0 |
| 0.200 | 3532 | 3270 | 1401 | 510 |
| 0.400 | 7063 | 6068 | 2600 | 946 |
| 0.600 | 10595 | 8461 | 3626 | 1320 |
| 0.800 | 14126 | 10509 | 4504 | 1639 |
| 1.000 | 17658 | 12261 | 5255 | 1913 |
| 1.200 | 21190 | 13760 | 5897 | 2146 |

| D (m) | D/W | Wdot, C = 1.40 (kg/s/m) | Wdot, C = 1.45 (kg/s/m) | Wdot, C = 1.55 (kg/s/m) | emptying time at C = 1.45 (s) |
| --- | --- | --- | --- | --- | --- |
| 0.040 | 0.100 | 63.1 | 65.4 | 69.9 | 13.21 |
| 0.060 | 0.150 | 116.0 | 120.1 | 128.4 | 7.19 |
| 0.080 | 0.200 | 178.6 | 185.0 | 197.7 | 4.67 |
<!-- END silo_reference.py output -->

## 6. Tolerances for the tests

Reasoning is in `docs/decisions/0019-silo-reference-tolerances.md`.

| Package | Check | Reference | Tolerance |
| --- | --- | --- | --- |
| 3.1p | Smooth walls: bottom vertical stress | rho g H | 2 percent (unchanged) |
| 3.1p | Smooth walls: horizontal / vertical stress | nu/(1 - nu) = 0.4286 (elastic, no yield) | 5 percent (unchanged) |
| 3.1p | Frictional walls (20 deg), H = 3 W: mean sigma_v against depth below one width from the top | Janssen with the fitted K | 10 percent (unchanged) |
| 3.1p | Fitted K | Between the wall-yield value 0.29 and the elastic value 0.43 | 0.27 to 0.45; record the value |
| 3.1p | Wall shear / normal where the wall slips | tan 20 = 0.364 | 5 percent (unchanged) |
| 3.1p | Vertical force balance | weight = base + walls | 1 percent (unchanged) |
| 3.2p | Column a = 0.5: final height | H_f = H_i (top survives) | within 1 cell |
| 3.2p | Column a = 0.5: run-out (L_f - L_i)/L_i | 1.2 a = 0.60 | 0.5 to 1.5 times the fit: 0.30 to 0.90 |
| 3.2p | Column a = 2: run-out | 1.2 a = 2.4 | 0.5 to 1.5 times the fit: 1.2 to 3.6 |
| 3.2p | Final deposit slope | below the friction angle | max surface slope at most phi + 2 deg = 37 deg (one cell of slack) |
| 3.3p | Rate exponent over D = 0.04, 0.06, 0.08 m | 1.5 | 1.3 to 1.7 (unchanged) |
| 3.3p | Rate constant C (exponent fixed at 1.5, k = 0) | 1.45 (continuum mu(I), 1.40 to 1.55) | 1.0 to 1.9 (30 percent); if outside, record and do not tune |
| 3.3p | Head independence, H = 0.6 and 0.9 m | S2: within 4 percent | 10 percent (unchanged) |
| 3.3p | Steady rate while the bed is deeper than W | S2 | variation below 10 percent (unchanged) |
| 3.3p | Stagnant corners | S2 | velocity below 1 percent of the outlet velocity (unchanged) |
| 3.5 | Cohesive arching | **Not collected yet**: the critical outlet width for a cohesive arch (Jenike) is to be sourced in the 3.5 brief | |
