# Gas references for version 2 (WP F.0)

Sources opened 2026-11-24 unless marked "not read". Values in the generated block come
from `tools/gas_reference.py`; `tests/test_gas_reference.sh` checks that this page matches
the script and that the closed forms pass their checks. Plan: decision 0125.

## Closed forms used by the F-series tests

| Quantity | Formula | Used in | Source |
| --- | --- | --- | --- |
| Air density | ideal gas, R = 287.05 J/(kg K), 20 C, 101325 Pa | all | standard |
| Air viscosity | Sutherland, mu0 1.716e-5 Pa s, T0 273.15 K, S 110.4 K | all | standard constants (book not opened); checked against the tabulated 1.81e-5 Pa s at 20 C |
| Packed-bed drag | Ergun: -dp/dx = 150 mu (1-e)^2/(e^3 d^2) U + 1.75 rho_g (1-e)/(e^3 d) U^2 | F.2, F.3, F.5 | Ergun (1952), constants as written in [El Geitani, Golshan and Blais (2022)](https://arxiv.org/pdf/2206.02842), eq. 48 |
| Pore pressure diffusivity | rigid bed c = k p0/(mu e); deformable c = k/(mu (e/p0 + 1/M)) | F.2, F.3 | linearised isothermal gas mass balance with Darcy flow (derivation in decision 0125) |
| 1D consolidation | Terzaghi series U(T), T = c t/H^2 | F.3 | closed form, checked against finite differences |
| Step at a face | p - p0 = dp erfc(x/(2 sqrt(c t))) | F.2, F.4 | closed form, checked against finite differences |
| Minimum fluidisation | Re_mf = sqrt(33.7^2 + 0.0408 Ar) - 33.7 | F.7 | Wen and Yu (1966), AIChE J. 12(3) 610-612 (not read); constants as written in [El Geitani et al. (2022)](https://arxiv.org/pdf/2206.02842), eq. 45 |
| Drag, default and for every test (Sia, 2026-11-24) | Ergun-Wen-Yu (Gidaspow): Ergun below gas fraction 0.8, Wen and Yu with Schiller-Naumann Cd and e^-2.65 above | F.2 to F.9 | Gidaspow (1994), Multiphase Flow and Fluidization (not read; standard form); checked: dense branch equals Ergun to 1e-12, dilute limit one Schiller-Naumann sphere |
| Drag, packed to dilute (user option) | Di Felice: F = F_single e^(-beta), beta = 3.7 - 0.65 exp(-(1.5 - log10 Re)^2/2), Cd = (0.63 + 4.8/sqrt(Re))^2 | F.7, F.8 | [Di Felice (1994), Int. J. Multiphase Flow 20(1) 153-159](https://www.scilit.net/publications/ce60d2e801ff2d196451c46a19d06d52) (paper not read); formula as written in [El Geitani et al. (2022)](https://arxiv.org/pdf/2206.02842), eq. 30 |
| Hindered settling | Richardson and Zaki, U_s = U_t e^n, n from Rowe (1987) | F.7 | Richardson and Zaki (1954), Rowe (1987) (not read; standard form) |

## Experiments for validation

| Package | Case | What the source gives | Source |
| --- | --- | --- | --- |
| F.5 | Air-impeded hopper discharge | Record only: Powder Technology 16(2) 197-207 (1977), no abstract online; values to be read in F.5 if the paper can be reached | [Crewdson, Ormond and Nedderman (1977)](https://www.scilit.net/publications/44d4c06ce61dcd19ac55a6bb28eae6d3), doi 10.1016/0032-5910(77)87007-1 |
| F.5 | Hourglass with interstitial fluid | Abstract: the fluid matters when N = (single-grain terminal velocity) / (hopper discharge velocity) is below 10; rate depends on the Archimedes number and D/d; flow stops for d above 25 % of the effective outlet | [Muite, Hunt and Joseph (2004), Phys. Fluids 16(9) 3415-3425](https://authors.library.caltech.edu/1793) |
| F.8 | Initially fluidised dam break in gas | Page refused automated reading; to be read in F.8 | [Roche et al. (2010)](https://hal.archives-ouvertes.fr/hal-00381402) (not read) |
| F.8 backup | Collapse of a column in a liquid, loose against dense | Glass beads 225 microns, 2500 kg/m^3; viscosity 12 and 23 cP; phi_i 0.55 to 0.62, transition at 0.58; loose: long thin deposit, positive pore pressure; dense: run-out half as long, slow, negative pore pressure; extremes about 0.2 Delta rho g H_i | [Rondon, Aussillous and Pouliquen (2011)](https://core.ac.uk/download/15523693.pdf) |

Rondon is in a liquid, which v2 does not cover (gas only). It stays the backup for F.8
because the same pore-pressure equation holds with the liquid's compressibility.

## What the numbers say

The table shows why gas matters only for fine grains. Pore pressure in a 0.5 m bed of
0.1 mm sand takes about 1.5 s to relax, as long as a discharge; for 1 mm grains it takes
0.015 s, and for the 4 mm Couto wheat about 1 ms. Muite's criterion N < 10 should be
checked for each case before a run (F.5 records it for the Couto hopper too).

<!-- BEGIN gas_reference.py output -->
Air at 20 C, 101325 Pa: rho_g = 1.2041 kg/m^3, mu = 1.8133e-05 Pa s

Quartz grains (rho_p = 2650 kg/m^3), packed bed e = 0.40, drainage length 0.5 m:

| d [mm] | k [m^2] | c rigid [m^2/s] | H^2/c [s] | U_mf Wen-Yu [m/s] | U_mf Ergun [m/s] | U_t [m/s] | n R-Z |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.1 | 1.185e-11 | 1.656e-01 | 1.510e+00 | 8.667e-03 | 1.018e-02 | 5.333e-01 | 3.966 |
| 0.5 | 2.963e-10 | 4.139e+00 | 6.040e-02 | 1.976e-01 | 2.227e-01 | 3.535e+00 | 2.668 |
| 1.0 | 1.185e-09 | 1.656e+01 | 1.510e-02 | 5.593e-01 | 5.819e-01 | 6.188e+00 | 2.482 |
| 4.0 | 1.896e-08 | 2.649e+02 | 9.437e-04 | 1.754e+00 | 1.683e+00 | 1.521e+01 | 2.369 |

Ergun-Wen-Yu jump at e = 0.8 (known in Gidaspow's blend; d = 0.5 mm, slip 0.1 m/s): K = 6.2828e+02 below, 6.1010e+02 above [kg/(m^3 s)]

Terzaghi degree of consolidation U(T):

| T | 0.01 | 0.05 | 0.1 | 0.2 | 0.5 | 1.0 |
| --- | --- | --- | --- | --- | --- | --- |
| U | 0.1128 | 0.2523 | 0.3568 | 0.5041 | 0.7640 | 0.9313 |
<!-- END gas_reference.py output -->
