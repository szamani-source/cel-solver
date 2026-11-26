# Temperature references (WP T.0)

Sources opened 2026-11-26 unless marked "not read". Values in the generated block come
from `tools/thermal_reference.py`; `tests/test_thermal_reference.sh` checks that this
page matches the script and that every closed form passes its check against an
independent finite-difference solution or a published tabulated value. Plan:
decision 0135.

## Equations

Volume-averaged two-temperature model for packed and moving beds (Nield and Bejan,
*Convection in Porous Media*, not read, standard form; the same balances as the
two-fluid codes in Kuipers et al. 1992 and the papers below), e = gas fraction,
phi_s = 1 - e, q = Darcy flux of gas relative to the grains:

- grains: d(phi_s rho_s c_s T_s)/dt + div(phi_s rho_s c_s T_s u_s)
  = div(phi_s k_s,eff grad T_s) + h a_s (T_g - T_s) + beta sigma:d(eps_p)/dt
- gas: d(e rho_g c_g T_g)/dt + div(rho_g c_g T_g (e u_s + q))
  = div(e k_g,eff grad T_g) + h a_s (T_s - T_g)

The phase fraction appears once, in front of the effective conductivity (as written
in Armstrong et al. 2010, Table 1, and Mahmood et al., eqs. 7-8).

## Closed forms used by the T-series tests

| Quantity | Formula | Used in | Source |
| --- | --- | --- | --- |
| Air | ideal gas, cp 1006 J/(kg K), Sutherland viscosity (as `gas_reference.py`) and conductivity k0 0.0241 W/(m K), T0 273.15 K, S 194 K | all | Sutherland constants for k: White, *Viscous Fluid Flow* (not read; standard values); checked against Incropera Table A.4 at 300 K (k 0.0263, mu 1.846e-5, Pr 0.707) within 1 % |
| Step at a face | T - T0 = dT erfc(x / (2 sqrt(a t))) | T.1 | Carslaw and Jaeger (1959), *Conduction of Heat in Solids* (not read; standard); checked against finite differences |
| Constant heat flux at a face | T - T0 = 2 q/k sqrt(a t/pi) exp(-x^2/(4 a t)) - q x/k erfc(x/(2 sqrt(a t))) | T.1 | Carslaw and Jaeger (1959), section 2.9 (not read; standard); checked against finite differences |
| Slab cooled by convection | theta = sum C_n exp(-z_n^2 Fo) cos(z_n x/L), z_n tan z_n = Bi, C_n = 4 sin z_n/(2 z_n + sin 2 z_n) | T.3 | Carslaw and Jaeger (1959); first roots and C_1 checked against Incropera et al., *Fundamentals of Heat and Mass Transfer*, Table 5.1 (Bi 1: 0.8603, 1.1191; Bi 10: 1.4289, 1.2620), series against finite differences |
| Stagnant bed conductivity | Zehner and Schluender (1970), split between the phases as in Kuipers et al. (1992): e k_g,eff = (1 - sqrt(1-e)) k_g; (1-e) k_s,eff = sqrt(1-e) (omega A + (1-omega) C) k_g; A = k_p/k_g; B = 1.25 ((1-e)/e)^(10/9); omega = 7.26e-3; C = 2/(1-B/A) [(A-1) B/A/(1-B/A)^2 ln(A/B) - (B-1)/(1-B/A) - (B+1)/2] | T.3 | Zehner and Schluender (1970), Chem. Ing. Tech. 42, 933; Kuipers, Prins and van Swaaij (1992), AIChE J. 38, 1079-1091 (both not read). Formulas as written in [Armstrong, Gu and Luo (2010), Int. J. Multiphase Flow 36, 916-929](https://eprints.soton.ac.uk/165821/1/published_version.pdf), eqs. 4-7, and [Mahmood, Zhang and Koksal](https://yunus.hacettepe.edu.tr/~koksalm/FBC19_CFD_KOKSAL.pdf), eq. 18 (the two agree). Checks: k_p = k_g gives k_g; values inside the series and parallel bounds |
| Grain-gas heat transfer | Gunn (1978): Nu = (7 - 10e + 5e^2)(1 + 0.7 Re^0.2 Pr^(1/3)) + (1.33 - 2.4e + 1.2e^2) Re^0.7 Pr^(1/3); Re = e rho_g abs(u_g - u_s) d/mu; h = Nu k_g/d; a_s = 6 (1-e)/d | T.4 | Gunn (1978), Int. J. Heat Mass Transfer 21, 467-476 (not read); formula as written in [Yusuf, Melaaen and Mathiesen (2005)](https://www.sintef.no/project/CFD2008/CFD2005_papers/42_Yusuf_Melaaen_Mathiesen.pdf), eqs. 15-16. Checks: single sphere in still gas Nu = 2 |
| Packed-bed thermal front | Schumann (1929): theta_g = J(xi, tau), theta_s = 1 - J(tau, xi), J(x, y) = 1 - exp(-y) int_0^x exp(-s) I0(2 sqrt(y s)) ds; xi = h a_s x/(e rho_g c_g u), tau = h a_s (t - x/u)/((1-e) rho_s c_s), u interstitial | T.4 | Schumann (1929), J. Franklin Inst. 208, 405-416 (not read; the J-function form is standard). Checked against an independent finite-difference solution to 2e-6 |
| Local equilibrium front | w = e rho_g c_g u/(e rho_g c_g + (1-e) rho_s c_s) | T.4 | energy balance; checked as the far-field limit of Schumann (tau_50 = xi within 0.25 % at xi = 200) |
| Laminar flow between heated plates | fully developed Nu_Dh = 7.541 (wall temperature fixed), D_h = 4 b | T.5 | Shah and London (1978), *Laminar Flow Forced Convection in Ducts* (not read; standard value); computed in the script by shooting: 7.5407 |
| Moving bed past a heated wall | slug flow between plates, wall flux: Nu_Dh(z) = [1/12 - 1/2 sum exp(-16 n^2 pi^2/Gz)/(n pi)^2 + R_nw k_eff/(4b)]^-1, Gz = (D_h/z) U D_h/a; near-wall air gap R_nw = D_air/k_air | T.6 | Muzychka et al. and Kakac et al. (not read), as written in [Adapa et al. (2024), arXiv 2403.19892](https://arxiv.org/pdf/2403.19892), eqs. 1-10; checked against an independent marching solution to 1e-4; fully developed 12 |
| Near-wall gap | D_air = 0.085 d (Sullivan and Sabersky 1975); 0.09 d smooth and 0.114 d rough (Natarajan and Hunt); 0.1 d (Denloye and Botterill) | T.6 | as quoted by Adapa et al. (2024); Sullivan and Sabersky (1975), Int. J. Heat Mass Transfer 18, 97-107 (not read); [Sullivan's thesis (Caltech 1973)](https://thesis.caltech.edu/1146/) is a scan without text |
| Adiabatic plastic heating | dT = beta W_p/(rho c) | T.7 | energy balance |
| Sealed box of ideal gas | p/T constant | T.7 | ideal gas law |

## Experiments for validation

| Package | Case | What the source gives | Source |
| --- | --- | --- | --- |
| T.6 | Ceramic particles flowing at 5-15 mm/s down a 5 mm gap between heated plates (300 x 30 x 5 mm), 300-650 C | Measured bed k_eff and near-wall gap D_air for three particle sizes (275, 404, 956 microns, Tables 1-6), average wall heat transfer coefficients 225-350 W/(m^2 K) (figures only), and the slug-flow model above fitted to them | [Adapa et al. (2024), arXiv 2403.19892](https://arxiv.org/pdf/2403.19892) |
| T.6 | Glass beads and sands flowing past a heated surface | Semi-empirical correlation with the 0.085 d gap; thesis is a scan without text | Sullivan and Sabersky (1975) (not read) |

The Adapa et al. (2024) case is the T.6 validation: their measured k_eff and D_air are
inputs, their slug-flow formula is the analytical check, and their 225-350 W/(m^2 K)
range is the experimental check. It replaces the plan's Sullivan and Sabersky (1975)
data, which could not be read.

## What the numbers say

Conduction never limits the time step: the stable explicit step is about 20 s on 4 mm
cells against a solid step near 1e-5 s. Grain-gas exchange is the stiff part: for
0.1 mm grains the gas reaches the grain temperature in about 1e-5 s, the solid's own
step, so the exchange is solved implicitly and exactly per cell. For 4 mm grains the
gas takes 0.01 s and the grains 30 s, so slow flows of coarse grains are close to local
equilibrium only in the gas.

<!-- BEGIN thermal_reference.py output -->
Air at 20 C, 101325 Pa: rho_g = 1.2041 kg/m^3, cp = 1006 J/(kg K), k = 0.0257 W/(m K), mu = 1.8133e-05 Pa s, Pr = 0.710

Stagnant bed conductivity (Zehner-Schluender, Kuipers split), still air:

| grains | k_p [W/(m K)] | e | e k_g,eff | (1-e) k_s,eff | bed k_eff [W/(m K)] |
| --- | --- | --- | --- | --- | --- |
| glass beads | 1.0 | 0.36 | 0.0051 | 0.1876 | 0.1927 |
| glass beads | 1.0 | 0.40 | 0.0058 | 0.1654 | 0.1712 |
| glass beads | 1.0 | 0.45 | 0.0066 | 0.1412 | 0.1478 |
| quartz sand | 3.0 | 0.36 | 0.0051 | 0.2863 | 0.2915 |
| quartz sand | 3.0 | 0.40 | 0.0058 | 0.2485 | 0.2542 |
| quartz sand | 3.0 | 0.45 | 0.0066 | 0.2084 | 0.2150 |

Grain-gas heat transfer (Gunn), e = 0.40, slip 0.01 m/s, quartz sand, and the time scales that set the numerics (4 mm cells):

| d [mm] | Re | Nu | h [W/(m^2 K)] | h a_s [W/(m^3 K)] | gas relax [s] | solid relax [s] | conduction step [s] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.1 | 0.027 | 4.988 | 1282 | 4.614e+07 | 1.050e-05 | 2.860e-02 | 20.8 |
| 0.5 | 0.133 | 5.507 | 283 | 2.037e+06 | 2.378e-04 | 6.477e-01 | 20.8 |
| 1.0 | 0.266 | 5.819 | 149.5 | 5.382e+05 | 9.003e-04 | 2.452e+00 | 20.8 |
| 4.0 | 1.062 | 6.725 | 43.2 | 3.888e+04 | 1.246e-02 | 3.394e+01 | 20.8 |

Slab cooled on both faces, centre (x = 0) and face (x = L):

| Bi | z_1 | C_1 | theta centre Fo 0.2 | Fo 0.5 | Fo 1.0 | theta face Fo 0.5 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.8603 | 1.1191 | 0.95064 | 0.77253 | 0.53386 | 0.50452 |
| 10 | 1.4289 | 1.2620 | 0.82925 | 0.45464 | 0.16382 | 0.06433 |

Schumann packed bed (hot gas in, cold bed), theta = (T - T_0)/(T_in - T_0):

| xi | tau | theta_gas | theta_solid |
| --- | --- | --- | --- |
| 2 | 1 | 0.39430 | 0.18258 |
| 5 | 5 | 0.56392 | 0.43608 |
| 10 | 5 | 0.11979 | 0.07439 |
| 10 | 10 | 0.54489 | 0.45511 |
| 20 | 20 | 0.53164 | 0.46836 |

Channel flow between plates, Nu_Dh (D_h = 4 b):

| case | value |
| --- | --- |
| laminar, wall temperature fixed, fully developed (computed) | 7.5407 |
| slug flow, wall flux, no gap, Gz = 10 | 12.0000 |
| slug flow, wall flux, no gap, Gz = 100 | 13.7237 |
| slug flow, wall flux, no gap, Gz = 1000 | 31.5632 |
| slug flow, wall flux, fully developed, no gap | 12.0000 |
<!-- END thermal_reference.py output -->
