#!/usr/bin/env python3
"""Closed-form reference values for version 2, gas in the pores and the void (WP F.0).

Prints the numbers tabulated in docs/spec/gas-references.md and checks the closed
forms that the F-series tests compare against. Pure Python, no OpenFOAM needed.
Pressures here are gauge or absolute as named; COMPRESSION POSITIVE for pore
pressure p (the solver keeps stress tension positive and adds -p I).

Usage:
    python3 tools/gas_reference.py           # tables for docs/spec/gas-references.md
    python3 tools/gas_reference.py --check   # internal consistency checks

Formulas (sources in docs/spec/gas-references.md):
    Air, ideal gas:       rho_g = p / (R_air T),  R_air = 287.05 J/(kg K)
    Air viscosity:        Sutherland, mu = mu0 (T/T0)^1.5 (T0 + S)/(T + S),
                          mu0 = 1.716e-5 Pa s, T0 = 273.15 K, S = 110.4 K
    Ergun (1952):         -dp/dx = 150 mu (1-e)^2/(e^3 d^2) U + 1.75 rho_g (1-e)/(e^3 d) U^2
                          U = superficial gas velocity relative to the grains, e = gas fraction
    Darcy permeability:   k = e^3 d^2 / (150 (1-e)^2)   (viscous term of Ergun)
    Pore pressure diffusivity, rigid bed, isothermal gas (linearised about p0):
                          c = k p0 / (mu e)
    Deformable bed (constrained modulus M of the skeleton):
                          c = k / (mu (e/p0 + 1/M))
    Terzaghi (1D, drained top, sealed bottom, layer H):
                          U(T) = 1 - sum 8/((2m+1)^2 pi^2) exp(-(2m+1)^2 pi^2 T/4),  T = c t / H^2
    Step at a face of a half-space:  p - p0 = dp erfc(x / (2 sqrt(c t)))
    Wen and Yu (1966):    Re_mf = sqrt(33.7^2 + 0.0408 Ar) - 33.7,
                          Ar = d^3 rho_g (rho_p - rho_g) g / mu^2,  U_mf = Re_mf mu / (rho_g d)
    Ergun-Wen-Yu drag (Gidaspow 1994), the default and the law every F test uses.
    Momentum exchange K in  phi_g grad p = -K (u_g - u_s), e = phi_g, u = |u_g - u_s|:
        e <  0.8: K = 150 (1-e)^2 mu/(e d^2) + 1.75 (1-e) rho_g u/d          (Ergun)
        e >= 0.8: K = 0.75 Cd e (1-e) rho_g u/d e^-2.65,                     (Wen and Yu)
                  Cd = 24/Re (1 + 0.15 Re^0.687) (Re < 1000), 0.44 above, Re = e rho_g u d/mu
    Single sphere drag:   Cd = (0.63 + 4.8/sqrt(Re))^2  (Dallavalle, used by Di Felice)
    Di Felice (1994):     F = F_single(e) e^(-beta), beta = 3.7 - 0.65 exp(-(1.5 - log10 Re)^2 / 2)
    Richardson and Zaki (1954), exponent of Rowe (1987):
                          U_s = U_t e^n,  n = (4.7 + 0.41 Re_t^0.75)/(1 + 0.175 Re_t^0.75)
"""

import argparse
import math
import sys

G = 9.81
R_AIR = 287.05
P_ATM = 101325.0
T_REF = 293.15
MU0, T0_SUTH, S_SUTH = 1.716e-5, 273.15, 110.4


def air_density(p=P_ATM, T=T_REF):
    return p / (R_AIR * T)


def air_viscosity(T=T_REF):
    return MU0 * (T / T0_SUTH) ** 1.5 * (T0_SUTH + S_SUTH) / (T + S_SUTH)


def ergun_gradient(U, d, e, mu, rho_g):
    """Pressure drop per length [Pa/m] for superficial velocity U [m/s]."""
    a = 150.0 * mu * (1.0 - e) ** 2 / (e ** 3 * d ** 2)
    b = 1.75 * rho_g * (1.0 - e) / (e ** 3 * d)
    return a * U + b * U * abs(U)


def permeability(d, e):
    return e ** 3 * d ** 2 / (150.0 * (1.0 - e) ** 2)


def diffusivity_rigid(d, e, mu, p0=P_ATM):
    return permeability(d, e) * p0 / (mu * e)


def diffusivity_deformable(d, e, mu, M, p0=P_ATM):
    return permeability(d, e) / (mu * (e / p0 + 1.0 / M))


def terzaghi_degree(T, terms=200):
    s = 0.0
    for m in range(terms):
        a = (2 * m + 1) * math.pi / 2.0
        s += 2.0 / a ** 2 * math.exp(-a * a * T)
    return 1.0 - s


def terzaghi_pressure(z, T, terms=200):
    """Excess pressure / initial excess at depth z/H (0 = drained top, 1 = sealed base)."""
    s = 0.0
    for m in range(terms):
        a = (2 * m + 1) * math.pi / 2.0
        s += 2.0 / a * math.sin(a * z) * math.exp(-a * a * T)
    return s


def step_pressure(x, t, c):
    return math.erfc(x / (2.0 * math.sqrt(c * t)))


def archimedes(d, rho_p, rho_g, mu):
    return d ** 3 * rho_g * (rho_p - rho_g) * G / mu ** 2


def wen_yu_umf(d, rho_p, rho_g, mu):
    ar = archimedes(d, rho_p, rho_g, mu)
    re = math.sqrt(33.7 ** 2 + 0.0408 * ar) - 33.7
    return re * mu / (rho_g * d)


def cd_schiller_naumann(re):
    return 24.0 / re * (1.0 + 0.15 * re ** 0.687) if re < 1000.0 else 0.44


def drag_ergun_wen_yu(e, u, d, mu, rho_g):
    """Gidaspow momentum exchange coefficient K [kg/(m^3 s)], u = slip speed."""
    if e < 0.8:
        return 150.0 * (1 - e) ** 2 * mu / (e * d * d) + 1.75 * (1 - e) * rho_g * u / d
    re = max(e * rho_g * u * d / mu, 1e-300)
    return 0.75 * cd_schiller_naumann(re) * e * (1 - e) * rho_g * u / d * e ** -2.65


def cd_dallavalle(re):
    return (0.63 + 4.8 / math.sqrt(re)) ** 2


def di_felice_beta(re):
    return 3.7 - 0.65 * math.exp(-((1.5 - math.log10(re)) ** 2) / 2.0)


def terminal_velocity(d, rho_p, rho_g, mu):
    """Single sphere, Dallavalle drag, by bisection on U."""
    def f(u):
        re = rho_g * u * d / mu
        return 0.75 * cd_dallavalle(re) * rho_g * u * u / d - (rho_p - rho_g) * G
    lo, hi = 1e-12, 1e3
    for _ in range(300):
        mid = math.sqrt(lo * hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def rz_exponent(re_t):
    r = re_t ** 0.75
    return (4.7 + 0.41 * r) / (1.0 + 0.175 * r)


def ergun_umf(d, rho_p, rho_g, mu, e_mf):
    """Minimum fluidisation from Ergun = bed weight per volume (check of Wen-Yu)."""
    w = (1.0 - e_mf) * (rho_p - rho_g) * G
    a = 150.0 * mu * (1.0 - e_mf) ** 2 / (e_mf ** 3 * d ** 2)
    b = 1.75 * rho_g * (1.0 - e_mf) / (e_mf ** 3 * d)
    return (-a + math.sqrt(a * a + 4.0 * b * w)) / (2.0 * b)


# --------------------------------------------------------------------------- tables

GRAINS = [0.1e-3, 0.5e-3, 1.0e-3, 4.0e-3]   # m: fine sand to wheat (Couto)
RHO_P = 2650.0                              # quartz grain density
E_BED = 0.40                                # gas fraction of a packed bed
H_BED = 0.5                                 # m, drainage length for the time scale


def tables(out=sys.stdout):
    mu, rg = air_viscosity(), air_density()
    w = out.write
    w("<!-- BEGIN gas_reference.py output -->\n")
    w("Air at 20 C, 101325 Pa: rho_g = %.4f kg/m^3, mu = %.4e Pa s\n\n" % (rg, mu))
    w("Quartz grains (rho_p = %.0f kg/m^3), packed bed e = %.2f, drainage length %.1f m:\n\n"
      % (RHO_P, E_BED, H_BED))
    w("| d [mm] | k [m^2] | c rigid [m^2/s] | H^2/c [s] | U_mf Wen-Yu [m/s] | U_mf Ergun [m/s] | U_t [m/s] | n R-Z |\n")
    w("| --- | --- | --- | --- | --- | --- | --- | --- |\n")
    for d in GRAINS:
        k = permeability(d, E_BED)
        c = diffusivity_rigid(d, E_BED, mu)
        ut = terminal_velocity(d, RHO_P, rg, mu)
        w("| %.1f | %.3e | %.3e | %.3e | %.3e | %.3e | %.3e | %.3f |\n" % (
            d * 1e3, k, c, H_BED ** 2 / c, wen_yu_umf(d, RHO_P, rg, mu),
            ergun_umf(d, RHO_P, rg, mu, E_BED), ut, rz_exponent(rg * ut * d / mu)))
    w("\nErgun-Wen-Yu jump at e = 0.8 (known in Gidaspow's blend; d = 0.5 mm, slip 0.1 m/s):"
      " K = %.4e below, %.4e above [kg/(m^3 s)]\n" % (
          drag_ergun_wen_yu(0.8 - 1e-9, 0.1, 0.5e-3, mu, rg),
          drag_ergun_wen_yu(0.8, 0.1, 0.5e-3, mu, rg)))
    w("\nTerzaghi degree of consolidation U(T):\n\n")
    w("| T | 0.01 | 0.05 | 0.1 | 0.2 | 0.5 | 1.0 |\n| --- | --- | --- | --- | --- | --- | --- |\n")
    w("| U | " + " | ".join("%.4f" % terzaghi_degree(T) for T in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0))
      + " |\n")
    w("<!-- END gas_reference.py output -->\n")


# --------------------------------------------------------------------------- checks

def check():
    ok = True

    def report(name, err, tol):
        nonlocal ok
        good = err <= tol
        ok = ok and good
        print("%s  %-58s err %.2e (tol %.0e)" % ("ok " if good else "BAD", name, err, tol))

    mu, rg = air_viscosity(), air_density()
    # Air values against the common tabulated 20 C values (1.204 kg/m^3, 1.81e-5 Pa s).
    report("air density at 20 C vs 1.204 kg/m^3", abs(rg - 1.204) / 1.204, 1e-3)
    report("air viscosity at 20 C vs 1.81e-5 Pa s", abs(mu - 1.81e-5) / 1.81e-5, 5e-3)

    # Ergun -> Darcy at low velocity.
    d, e = 0.1e-3, 0.4
    U = 1e-6
    darcy = mu * U / permeability(d, e)
    report("Ergun equals Darcy at U = 1e-6 m/s", abs(ergun_gradient(U, d, e, mu, rg) - darcy) / darcy, 1e-4)

    # Terzaghi series: U(0+) = 0, U(inf) = 1, early time 2 sqrt(T/pi).
    report("Terzaghi U(T = 1e-4) vs 2 sqrt(T/pi)",
           abs(terzaghi_degree(1e-4, 4000) - 2 * math.sqrt(1e-4 / math.pi)) / (2 * math.sqrt(1e-4 / math.pi)), 1e-3)
    report("Terzaghi U(T = 3) vs 1", abs(terzaghi_degree(3.0) - 1.0), 1e-3)
    # Terzaghi degree = 1 - average of the pressure profile (midpoint rule).
    n = 2000
    avg = sum(terzaghi_pressure((i + 0.5) / n, 0.2) for i in range(n)) / n
    report("Terzaghi U(0.2) vs 1 - mean pressure profile", abs((1 - avg) - terzaghi_degree(0.2)), 1e-6)

    # Terzaghi against an explicit finite-difference solution of p_t = p_zz.
    m = 200
    dz = 1.0 / m
    p = [1.0] * (m + 1)
    p[0] = 0.0
    dt = 0.4 * dz * dz
    t = 0.0
    while t < 0.2 - 1e-12:
        h = min(dt, 0.2 - t)
        q = p[:]
        for i in range(1, m):
            q[i] = p[i] + h / dz ** 2 * (p[i + 1] - 2 * p[i] + p[i - 1])
        q[m] = p[m] + h / dz ** 2 * 2 * (p[m - 1] - p[m])
        q[0] = 0.0
        p = q
        t += h
    report("Terzaghi p(z = 1, T = 0.2) vs finite differences",
           abs(p[m] - terzaghi_pressure(1.0, 0.2)), 2e-4)

    # erfc step against the same scheme on a long domain.
    L, m = 1.0, 400
    dx = L / m
    c = 1.0
    p = [0.0] * (m + 1)
    p[0] = 1.0
    dt = 0.4 * dx * dx / c
    t, tend = 0.0, 0.01
    while t < tend - 1e-15:
        h = min(dt, tend - t)
        q = p[:]
        for i in range(1, m):
            q[i] = p[i] + c * h / dx ** 2 * (p[i + 1] - 2 * p[i] + p[i - 1])
        p = q
        t += h
    i = int(0.1 / dx)
    report("erfc step at x = 0.1, c t = 0.01 vs finite differences",
           abs(p[i] - step_pressure(0.1, tend, c)), 1e-3)

    # Wen-Yu is Ergun at the bed weight with their bed-averaged shape factors
    # 1/(phi_s e^3) = 14 and (1 - e)/(phi_s^2 e^3) = 11: Ar = 24.5 Re^2 + 1650 Re.
    # Solving that quadratic must give their constants (33.7, 0.0408) back.
    for d in (0.1e-3, 0.5e-3, 4e-3):
        ar = archimedes(d, RHO_P, rg, mu)
        re = (-1650.0 + math.sqrt(1650.0 ** 2 + 4 * 24.5 * ar)) / (2 * 24.5)
        wy = wen_yu_umf(d, RHO_P, rg, mu) * rg * d / mu
        report("Wen-Yu Re_mf vs Ergun with shape factors 14, 11, d = %.1f mm" % (d * 1e3),
               abs(wy - re) / re, 2e-3)
    # With spheres at e_mf = 0.40 (the table) the two U_mf differ by 4 to 15 %.
    for d in GRAINS:
        wy = wen_yu_umf(d, RHO_P, rg, mu)
        er = ergun_umf(d, RHO_P, rg, mu, E_BED)
        report("Wen-Yu vs Ergun (spheres, e 0.40) U_mf, d = %.1f mm" % (d * 1e3), abs(wy - er) / er, 0.2)

    # Ergun-Wen-Yu (Gidaspow): dense branch is Ergun exactly (grad p = K U / e^2
    # with superficial U = e u); dilute branch tends to one Schiller-Naumann sphere
    # per particle volume as e -> 1.
    d, e = 0.5e-3, 0.4
    for U in (1e-4, 0.1, 1.0):
        u = U / e
        dpdx = drag_ergun_wen_yu(e, u, d, mu, rg) * u / e
        report("Ergun-Wen-Yu dense branch vs Ergun, U = %g m/s" % U,
               abs(dpdx - ergun_gradient(U, d, e, mu, rg)) / ergun_gradient(U, d, e, mu, rg), 1e-12)
    e, u = 1.0 - 1e-9, 0.2
    re = e * rg * u * d / mu
    single = 0.75 * cd_schiller_naumann(re) * rg * u / d   # K / (1 - e) for one sphere
    report("Ergun-Wen-Yu dilute limit vs Schiller-Naumann sphere",
           abs(drag_ergun_wen_yu(e, u, d, mu, rg) / (1 - e) - single) / single, 1e-6)

    # Di Felice limits and Richardson-Zaki exponents.
    report("Di Felice beta at Re -> 0 vs 3.7", abs(di_felice_beta(1e-8) - 3.7), 1e-6)
    report("Di Felice beta at Re = 10^1.5 vs 3.05", abs(di_felice_beta(10 ** 1.5) - 3.05), 1e-12)
    report("Richardson-Zaki n at Re -> 0 vs 4.7", abs(rz_exponent(1e-12) - 4.7), 1e-6)
    report("Richardson-Zaki n at Re -> inf vs 2.343", abs(rz_exponent(1e12) - 0.41 / 0.175), 1e-3)

    # Terminal velocity: Stokes limit for a 20 micron grain.
    d = 20e-6
    stokes = (RHO_P - rg) * G * d * d / (18 * mu)
    report("terminal velocity, 20 micron, vs Stokes (Dallavalle ~ +)",
           abs(terminal_velocity(d, RHO_P, rg, mu) - stokes) / stokes, 0.2)

    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    if args.check:
        sys.exit(0 if check() else 1)
    tables()


if __name__ == "__main__":
    main()
