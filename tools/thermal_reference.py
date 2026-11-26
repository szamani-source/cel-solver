#!/usr/bin/env python3
"""Closed-form reference values for temperature in celFoam (WP T.0).

Prints the numbers tabulated in docs/spec/thermal-references.md and checks the closed
forms that the T-series tests compare against, each against an independent finite
difference solution or a published tabulated value. Pure Python, no OpenFOAM, no numpy.

Usage:
    python3 tools/thermal_reference.py           # tables for docs/spec/thermal-references.md
    python3 tools/thermal_reference.py --check   # internal consistency checks

Formulas (sources in docs/spec/thermal-references.md):
    Air:                  ideal gas R = 287.05 J/(kg K), cp = 1006 J/(kg K),
                          Sutherland viscosity (as tools/gas_reference.py) and
                          conductivity k = k0 (T/T0)^1.5 (T0 + S)/(T + S),
                          k0 = 0.0241 W/(m K), T0 = 273.15 K, S = 194 K
    Step at a face (T.1): T - T0 = dT erfc(x / (2 sqrt(a t))),  a = k/(rho c)
    Flux q into a face (T.1): T - T0 = 2 q/k sqrt(a t/pi) exp(-x^2/(4 a t))
                                       - q x/k erfc(x / (2 sqrt(a t)))
    Slab cooled by convection on both faces (T.3), half thickness L, Bi = h L/k:
                          theta = sum C_n exp(-z_n^2 Fo) cos(z_n x/L),  Fo = a t/L^2
                          z_n tan z_n = Bi,  C_n = 4 sin z_n / (2 z_n + sin 2 z_n)
    Zehner and Schluender (1970) stagnant bed, split between the phases as in
    Kuipers et al. (1992), e = gas fraction, A = k_particle/k_gas,
    B = 1.25 ((1-e)/e)^(10/9), omega = 7.26e-3:
                          e k_g,eff     = (1 - sqrt(1-e)) k_gas
                          (1-e) k_s,eff = sqrt(1-e) (omega A + (1-omega) C) k_gas
                          C = 2/(1-B/A) [ (A-1) B/A /(1-B/A)^2 ln(A/B)
                                          - (B-1)/(1-B/A) - (B+1)/2 ]
                          bed k_eff = e k_g,eff + (1-e) k_s,eff
    Gunn (1978) grain-gas heat transfer, Re = e rho_g |u_g - u_s| d / mu (superficial
    slip velocity), Pr = mu cp / k:
                          Nu = (7 - 10 e + 5 e^2)(1 + 0.7 Re^0.2 Pr^(1/3))
                             + (1.33 - 2.4 e + 1.2 e^2) Re^0.7 Pr^(1/3)
                          h = Nu k_gas / d,  a_s = 6 (1-e)/d
    Schumann (1929) packed bed, hot gas into a cold bed, no conduction (T.4):
                          xi  = h a_s x / (e rho_g cp u)            (u interstitial)
                          tau = h a_s (t - x/u) / ((1-e) rho_s c_s)
                          theta_g = J(xi, tau),  theta_s = 1 - J(tau, xi),
                          J(x, y) = 1 - exp(-y) int_0^x exp(-s) I0(2 sqrt(y s)) ds
    Local equilibrium front speed (T.4):
                          w = e rho_g cp u / (e rho_g cp + (1-e) rho_s c_s)
    Graetz, parallel plates, laminar, wall temperature fixed (T.5): fully developed
                          Nu_Dh = 7.541 (Shah and London 1978), D_h = 4 b, b = half gap
    Slug (plug) flow between plates with wall heat flux and a near-wall air gap
    (T.6, Muzychka et al.; as written in Adapa et al. 2024, eqs. 1-9):
                          Nu_Dh(z) = [1/12 - 1/2 sum 1/(n^2 pi^2) exp(-16 n^2 pi^2/Gz)
                                      + R_nw k_eff/(4 b)]^-1,
                          Gz = (D_h/z) Pe,  Pe = U D_h / a,  R_nw = D_air / k_air,
                          D_air = 0.085 d (Sullivan and Sabersky 1975, as quoted there)
    Adiabatic plastic heating (T.7):  dT = beta W_p / (rho c)
    Sealed rigid box of ideal gas (T.7):  p/T constant
"""

import argparse
import math
import sys

G = 9.81
R_AIR = 287.05
CP_AIR = 1006.0
P_ATM = 101325.0
T_REF = 293.15
MU0, T0_SUTH, S_SUTH = 1.716e-5, 273.15, 110.4
K0_AIR, S_K_AIR = 0.0241, 194.0
OMEGA_ZS = 7.26e-3
NU_GRAETZ_PLATES = 7.541   # Shah and London (1978), fully developed, wall T fixed


# ----------------------------------------------------------------------------- air

def air_density(p=P_ATM, T=T_REF):
    return p / (R_AIR * T)


def air_viscosity(T=T_REF):
    return MU0 * (T / T0_SUTH) ** 1.5 * (T0_SUTH + S_SUTH) / (T + S_SUTH)


def air_conductivity(T=T_REF):
    return K0_AIR * (T / T0_SUTH) ** 1.5 * (T0_SUTH + S_K_AIR) / (T + S_K_AIR)


def air_prandtl(T=T_REF):
    return air_viscosity(T) * CP_AIR / air_conductivity(T)


# ------------------------------------------------------------------- conduction

def erfc_step(x, t, a, T0=0.0, dT=1.0):
    return T0 + dT * math.erfc(x / (2.0 * math.sqrt(a * t)))


def flux_step(x, t, a, k, q, T0=0.0):
    """Half-space heated by a constant flux q at x = 0 (Carslaw and Jaeger 2.9)."""
    r = math.sqrt(a * t)
    return T0 + 2.0 * q / k * r / math.sqrt(math.pi) * math.exp(-x * x / (4.0 * r * r)) \
        - q * x / k * math.erfc(x / (2.0 * r))


def slab_roots(Bi, n=60):
    """First n roots of z tan z = Bi (one in each ((k) pi, (k + 1/2) pi))."""
    roots = []
    for k in range(n):
        lo, hi = k * math.pi + 1e-14, k * math.pi + 0.5 * math.pi - 1e-14
        f = lambda z: z * math.sin(z) - Bi * math.cos(z)
        flo = f(lo)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi = mid
        roots.append(0.5 * (lo + hi))
    return roots


def slab_theta(x_over_L, Fo, Bi, n=60):
    """(T - T_inf)/(T_i - T_inf) in a slab |x| <= L cooled on both faces."""
    s = 0.0
    for z in slab_roots(Bi, n):
        C = 4.0 * math.sin(z) / (2.0 * z + math.sin(2.0 * z))
        s += C * math.exp(-z * z * Fo) * math.cos(z * x_over_L)
    return s


def zs_phase_conductivities(e, k_gas, k_part):
    """(e k_g,eff, (1-e) k_s,eff) of Zehner-Schluender as split by Kuipers et al."""
    A = k_part / k_gas
    B = 1.25 * ((1.0 - e) / e) ** (10.0 / 9.0)
    r = B / A
    if abs(1.0 - r) < 1e-9:
        C = 1.0 if abs(A - 1.0) < 1e-12 else zs_C(A, B * (1 + 1e-6))
    else:
        C = zs_C(A, B)
    kg = (1.0 - math.sqrt(1.0 - e)) * k_gas
    ks = math.sqrt(1.0 - e) * (OMEGA_ZS * A + (1.0 - OMEGA_ZS) * C) * k_gas
    return kg, ks


def zs_C(A, B):
    r = B / A
    return 2.0 / (1.0 - r) * (
        (A - 1.0) * r / (1.0 - r) ** 2 * math.log(A / B)
        - (B - 1.0) / (1.0 - r)
        - 0.5 * (B + 1.0)
    )


def zs_bed(e, k_gas, k_part):
    kg, ks = zs_phase_conductivities(e, k_gas, k_part)
    return kg + ks


def gunn_nu(e, re, pr):
    return (7.0 - 10.0 * e + 5.0 * e * e) * (1.0 + 0.7 * re ** 0.2 * pr ** (1.0 / 3.0)) \
        + (1.33 - 2.4 * e + 1.2 * e * e) * re ** 0.7 * pr ** (1.0 / 3.0)


def gunn_h(e, slip, d, T=T_REF, p=P_ATM):
    mu, k, rho = air_viscosity(T), air_conductivity(T), air_density(p, T)
    re = e * rho * abs(slip) * d / mu
    return gunn_nu(e, re, air_prandtl(T)) * k / d


# ------------------------------------------------------------------- Schumann

def bessel_i0(z):
    term, s, k = 1.0, 1.0, 0
    q = 0.25 * z * z
    while True:
        k += 1
        term *= q / (k * k)
        s += term
        if term < 1e-17 * s:
            return s


def schumann_J(x, y, n=4000):
    """J(x, y) = 1 - exp(-y) int_0^x exp(-s) I0(2 sqrt(y s)) ds (Simpson)."""
    if x <= 0.0:
        return 1.0
    if n % 2:
        n += 1
    h = x / n
    s = 0.0
    for i in range(n + 1):
        si = i * h
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        s += w * math.exp(-si - y) * bessel_i0(2.0 * math.sqrt(y * si))
    return 1.0 - s * h / 3.0


def schumann_gas(xi, tau, n=4000):
    return schumann_J(xi, tau, n)


def schumann_solid(xi, tau, n=4000):
    return 1.0 - schumann_J(tau, xi, n)


def schumann_fd(xi_max, tau_max, n):
    """Independent check: trapezoidal (Crank-Nicolson) on the characteristic grid of
    d(theta_g)/d(xi) = theta_s - theta_g, d(theta_s)/d(tau) = theta_g - theta_s,
    theta_g(0, tau) = 1, theta_s(xi, 0) = 0. Returns (theta_g, theta_s) at the corner."""
    h, k = xi_max / n, tau_max / n
    # row j = tau level; march xi inside each tau level
    gs_prev = [0.0] * (n + 1)      # theta_s on tau = 0
    gg_prev = None
    for j in range(n + 1):
        gg = [0.0] * (n + 1)
        gs = [0.0] * (n + 1)
        gg[0] = 1.0
        if j == 0:
            gs = [0.0] * (n + 1)
            for i in range(1, n + 1):
                gg[i] = (gg[i - 1] * (1 - h / 2) + h / 2 * (gs[i] + gs[i - 1])) / (1 + h / 2)
        else:
            # node (i, j): unknowns gg[i], gs[i]
            # gs: (gs - gs_prev[i])/k = 0.5[(gg - gs) + (gg_prev[i] - gs_prev[i])]
            # gg: (gg - gg[i-1])/h   = 0.5[(gs - gg) + (gs[i-1] - gg[i-1])]
            gs[0] = (gs_prev[0] + k / 2 * (gg[0] + gg_prev[0] - gs_prev[0])) / (1 + k / 2)
            for i in range(1, n + 1):
                a11, a12 = 1 + k / 2, -k / 2
                b1 = gs_prev[i] + k / 2 * (gg_prev[i] - gs_prev[i])
                a21, a22 = -h / 2, 1 + h / 2
                b2 = gg[i - 1] + h / 2 * (gs[i - 1] - gg[i - 1])
                det = a11 * a22 - a12 * a21
                gs[i] = (b1 * a22 - a12 * b2) / det
                gg[i] = (a11 * b2 - a21 * b1) / det
        gg_prev, gs_prev = gg, gs
    return gg_prev[n], gs_prev[n]


def equilibrium_front_speed(e, u, rho_g, cp_g, rho_s, c_s):
    return e * rho_g * cp_g * u / (e * rho_g * cp_g + (1.0 - e) * rho_s * c_s)


# ---------------------------------------------------------- channel heat transfer

def graetz_plates_nu(n=4000):
    """Fully developed Nu_Dh, laminar flow between plates, wall temperature fixed:
    phi'' + mu (1 - eta^2) phi = 0, phi'(0) = 0, phi(1) = 0, Nu_Dh = 8 mu / 3."""
    def shoot(mu):
        h = 1.0 / n
        y, z, eta = 1.0, 0.0, 0.0
        f = lambda et, yy, zz: (zz, -mu * (1 - et * et) * yy)
        for _ in range(n):
            k1 = f(eta, y, z)
            k2 = f(eta + h / 2, y + h / 2 * k1[0], z + h / 2 * k1[1])
            k3 = f(eta + h / 2, y + h / 2 * k2[0], z + h / 2 * k2[1])
            k4 = f(eta + h, y + h * k3[0], z + h * k3[1])
            y += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            z += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            eta += h
        return y
    lo, hi = 1.0, 5.0
    flo = shoot(lo)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        fm = shoot(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 8.0 * 0.5 * (lo + hi) / 3.0


def slug_plates_nu_local(gz, rnw_keff_over_4b=0.0, terms=2000):
    s = sum(math.exp(-16.0 * (n * math.pi) ** 2 / gz) / (n * math.pi) ** 2
            for n in range(1, terms + 1))
    return 1.0 / (1.0 / 12.0 - 0.5 * s + rnw_keff_over_4b)


def slug_plates_nu_mean(dh_over_l, pe, rnw_keff_over_4b=0.0, terms=2000):
    """Mean Nu_Dh over a heated length L (Adapa et al. 2024, eq. 9): the wall
    minus mean temperature difference averaged over L, h = q / mean(dT)."""
    s = sum(dh_over_l * pe / (32.0 * (n * math.pi) ** 4)
            * (math.exp(-16.0 * (n * math.pi) ** 2 / (dh_over_l * pe)) - 1.0)
            for n in range(1, terms + 1))
    return 1.0 / (1.0 / 12.0 + s + rnw_keff_over_4b)


def slug_plates_h_mean(k_eff, d_air, k_air, rho_bulk, cp, u, b, length, terms=2000):
    """Mean wall heat transfer coefficient [W/(m^2 K)] of a moving bed between
    plates 2b apart, heated over length with a uniform flux (Adapa et al. 2024,
    eqs. 9-10): D_h = 4 b, Pe = U D_h rho c/k_eff, R_nw = d_air/k_air."""
    dh = 4.0 * b
    pe = u * dh * rho_bulk * cp / k_eff
    nu = slug_plates_nu_mean(dh / length, pe, d_air / k_air * k_eff / (4.0 * b), terms)
    return nu * k_eff / dh


def slug_plates_fd(zeta, n=200):
    """Independent check of the gap-free local Nu: plug flow between plates of half gap
    b = 1, theta_zeta = theta_yy (zeta = a z/(U b^2)), q = 1 into the bed at y = 1,
    symmetry at y = 0, theta = 0 at the inlet. Crank-Nicolson, cell-centred.
    Returns Nu_Dh = 4 b q / (k (T_wall - T_mean))."""
    h = 1.0 / n
    T = [0.0] * n
    nsteps = max(400, int(zeta / (0.25 * h)))
    dz = zeta / nsteps
    r = dz / (h * h)
    for _ in range(nsteps):
        # tridiagonal: (1 + r) T_i - r/2 (T_i-1 + T_i+1) = rhs; flux 1 at the wall
        a = [0.0] * n
        b = [0.0] * n
        c = [0.0] * n
        d = [0.0] * n
        for i in range(n):
            left = T[i - 1] if i > 0 else T[i]
            right = T[i + 1] if i < n - 1 else T[i]
            lap = left - 2 * T[i] + right
            d[i] = T[i] + 0.5 * r * lap
            if i == n - 1:
                d[i] += dz / h          # wall flux q = 1 (k = 1) into the last cell
            a[i] = -0.5 * r if i > 0 else 0.0
            c[i] = -0.5 * r if i < n - 1 else 0.0
            b[i] = 1 + 0.5 * r * ((1 if i > 0 else 0) + (1 if i < n - 1 else 0))
        for i in range(1, n):
            m = a[i] / b[i - 1]
            b[i] -= m * c[i - 1]
            d[i] -= m * d[i - 1]
        T[n - 1] = d[n - 1] / b[n - 1]
        for i in range(n - 2, -1, -1):
            T[i] = (d[i] - c[i] * T[i + 1]) / b[i]
    mean = sum(T) / n
    wall = T[n - 1] + 0.5 * h * 1.0   # half cell to the wall with gradient q/k = 1
    return 4.0 / (wall - mean)


# ------------------------------------------------------------------------ tables

GRAINS = [
    # name, rho_s [kg/m^3], c_s [J/(kg K)], k_particle [W/(m K)]
    ("glass beads", 2500.0, 840.0, 1.0),
    ("quartz sand", 2650.0, 830.0, 3.0),
]


def tables(out=sys.stdout):
    w = out.write
    w("<!-- BEGIN thermal_reference.py output -->\n")
    w("Air at 20 C, 101325 Pa: rho_g = %.4f kg/m^3, cp = %.0f J/(kg K), "
      "k = %.4f W/(m K), mu = %.4e Pa s, Pr = %.3f\n\n"
      % (air_density(), CP_AIR, air_conductivity(), air_viscosity(), air_prandtl()))
    w("Stagnant bed conductivity (Zehner-Schluender, Kuipers split), still air:\n\n")
    w("| grains | k_p [W/(m K)] | e | e k_g,eff | (1-e) k_s,eff | bed k_eff [W/(m K)] |\n")
    w("| --- | --- | --- | --- | --- | --- |\n")
    for name, rho, c, kp in GRAINS:
        for e in (0.36, 0.40, 0.45):
            kg, ks = zs_phase_conductivities(e, air_conductivity(), kp)
            w("| %s | %.1f | %.2f | %.4f | %.4f | %.4f |\n" % (name, kp, e, kg, ks, kg + ks))
    w("\nGrain-gas heat transfer (Gunn), e = 0.40, slip 0.01 m/s, quartz sand, and the "
      "time scales that set the numerics (4 mm cells):\n\n")
    w("| d [mm] | Re | Nu | h [W/(m^2 K)] | h a_s [W/(m^3 K)] | gas relax [s] | "
      "solid relax [s] | conduction step [s] |\n")
    w("| --- | --- | --- | --- | --- | --- | --- | --- |\n")
    e = 0.40
    rho_s, c_s, kp = GRAINS[1][1], GRAINS[1][2], GRAINS[1][3]
    keff = zs_bed(e, air_conductivity(), kp)
    dx = 0.004
    dt_cond = (1 - e) * rho_s * c_s * dx * dx / (2 * 2 * keff)
    for d in (0.1e-3, 0.5e-3, 1e-3, 4e-3):
        re = e * air_density() * 0.01 * d / air_viscosity()
        nu = gunn_nu(e, re, air_prandtl())
        h = nu * air_conductivity() / d
        ha = h * 6 * (1 - e) / d
        tg = e * air_density() * CP_AIR / ha
        ts = (1 - e) * rho_s * c_s / ha
        w("| %.1f | %.3f | %.3f | %.4g | %.4g | %.3e | %.3e | %.3g |\n"
          % (d * 1e3, re, nu, h, ha, tg, ts, dt_cond))
    w("\nSlab cooled on both faces, centre (x = 0) and face (x = L):\n\n")
    w("| Bi | z_1 | C_1 | theta centre Fo 0.2 | Fo 0.5 | Fo 1.0 | theta face Fo 0.5 |\n")
    w("| --- | --- | --- | --- | --- | --- | --- |\n")
    for Bi in (1.0, 10.0):
        z1 = slab_roots(Bi, 1)[0]
        C1 = 4 * math.sin(z1) / (2 * z1 + math.sin(2 * z1))
        w("| %.0f | %.4f | %.4f | %.5f | %.5f | %.5f | %.5f |\n"
          % (Bi, z1, C1, slab_theta(0, 0.2, Bi), slab_theta(0, 0.5, Bi),
             slab_theta(0, 1.0, Bi), slab_theta(1, 0.5, Bi)))
    w("\nSchumann packed bed (hot gas in, cold bed), theta = (T - T_0)/(T_in - T_0):\n\n")
    w("| xi | tau | theta_gas | theta_solid |\n| --- | --- | --- | --- |\n")
    for xi, tau in ((2, 1), (5, 5), (10, 5), (10, 10), (20, 20)):
        w("| %g | %g | %.5f | %.5f |\n" % (xi, tau, schumann_gas(xi, tau),
                                            schumann_solid(xi, tau)))
    w("\nChannel flow between plates, Nu_Dh (D_h = 4 b):\n\n")
    w("| case | value |\n| --- | --- |\n")
    w("| laminar, wall temperature fixed, fully developed (computed) | %.4f |\n"
      % graetz_plates_nu())
    for gz in (10.0, 100.0, 1000.0):
        w("| slug flow, wall flux, no gap, Gz = %g | %.4f |\n" % (gz, slug_plates_nu_local(gz)))
    w("| slug flow, wall flux, fully developed, no gap | %.4f |\n" % slug_plates_nu_local(1e-9))
    w("<!-- END thermal_reference.py output -->\n")


# ------------------------------------------------------------------------ checks

def check():
    ok = True

    def report(name, err, tol):
        nonlocal ok
        good = err <= tol
        ok &= good
        print("%s %-62s err %.2e (tol %.0e)" % ("ok " if good else "BAD", name, err, tol))

    # air properties against tabulated values at 300 K (Incropera, Table A.4:
    # k = 0.0263 W/(m K), mu = 184.6e-7 Pa s, Pr = 0.707)
    report("air k at 300 K vs 0.0263", abs(air_conductivity(300) / 0.0263 - 1), 1e-2)
    report("air mu at 300 K vs 1.846e-5", abs(air_viscosity(300) / 1.846e-5 - 1), 1e-2)
    report("air Pr at 300 K vs 0.707", abs(air_prandtl(300) / 0.707 - 1), 1.5e-2)

    # erfc step against explicit finite differences
    a, L, n = 1.0, 6.0, 600
    h = L / n
    dt = 0.4 * h * h / a
    T = [0.0] * n
    t = 0.0
    while t < 1.0 - 1e-12:
        step = min(dt, 1.0 - t)
        new = T[:]
        for i in range(n):
            left = 2.0 - T[0] if i == 0 else T[i - 1]   # face at x = 0 held at 1
            right = T[i] if i == n - 1 else T[i + 1]
            new[i] = T[i] + a * step / (h * h) * (left - 2 * T[i] + right)
        T, t = new, t + step
    err = max(abs(T[i] - erfc_step((i + 0.5) * h, 1.0, a)) for i in range(n))
    report("erfc step vs finite differences", err, 2e-3)

    # constant flux q = 1 (k = 1) into the half-space against finite differences
    T = [0.0] * n
    t = 0.0
    while t < 1.0 - 1e-12:
        step = min(dt, 1.0 - t)
        new = T[:]
        for i in range(n):
            left = T[0] + h * 1.0 if i == 0 else T[i - 1]   # -dT/dx = q/k at x = 0
            right = T[i] if i == n - 1 else T[i + 1]
            new[i] = T[i] + a * step / (h * h) * (left - 2 * T[i] + right)
        T, t = new, t + step
    err = max(abs(T[i] - flux_step((i + 0.5) * h, 1.0, a, 1.0, 1.0)) for i in range(n))
    report("constant-flux face vs finite differences", err / flux_step(0, 1.0, a, 1.0, 1.0),
           2e-3)

    # slab: roots against Incropera Table 5.1 and series against finite differences
    for Bi, z1, C1 in ((1.0, 0.8603, 1.1191), (10.0, 1.4289, 1.2620)):
        zz = slab_roots(Bi, 1)[0]
        CC = 4 * math.sin(zz) / (2 * zz + math.sin(2 * zz))
        report("slab Bi %g first root vs Incropera %.4f" % (Bi, z1), abs(zz - z1), 6e-5)
        report("slab Bi %g C_1 vs Incropera %.4f" % (Bi, C1), abs(CC - C1), 6e-4)
        n = 200
        h = 1.0 / n
        T = [1.0] * n
        dt = 0.2 * h * h
        t = 0.0
        while t < 0.3 - 1e-12:
            step = min(dt, 0.3 - t)
            new = T[:]
            for i in range(n):
                left = T[i] if i == 0 else T[i - 1]
                if i == n - 1:
                    # Robin face: -k dT/dx = h T at x = 1, ghost from the face value
                    Tf = T[i] / (1 + Bi * h / 2)
                    right = 2 * Tf - T[i]
                else:
                    right = T[i + 1]
                new[i] = T[i] + step / (h * h) * (left - 2 * T[i] + right)
            T, t = new, t + step
        err = max(abs(T[i] - slab_theta((i + 0.5) * h, 0.3, Bi)) for i in range(n))
        report("slab Bi %g series vs finite differences, Fo 0.3" % Bi, err, 2e-4)

    # Zehner-Schluender: equal conductivities give k_gas; formula continuous near A = B
    kg = air_conductivity()
    report("ZS bed with k_p = k_gas equals k_gas", abs(zs_bed(0.4, kg, kg) / kg - 1), 1e-12)
    for e in (0.36, 0.40, 0.45):
        report("ZS phase split adds to the bed value, e %.2f" % e,
               abs(sum(zs_phase_conductivities(e, kg, 3.0)) - zs_bed(e, kg, 3.0)), 1e-15)
    # bounds: between the series and parallel arrangements
    for e in (0.36, 0.40, 0.45):
        kb = zs_bed(e, kg, 3.0)
        ser = 1.0 / (e / kg + (1 - e) / 3.0)
        par = e * kg + (1 - e) * 3.0
        report("ZS bed inside series/parallel bounds, e %.2f" % e,
               0.0 if ser < kb < par else 1.0, 0.0)

    # Gunn: single sphere in still gas (e -> 1, Re -> 0) gives Nu = 2
    report("Gunn Nu, e = 1, Re = 0 equals 2", abs(gunn_nu(1.0, 0.0, 0.7) - 2.0), 1e-14)
    report("Gunn Nu, e = 0.4, Re = 0 equals 7-4+0.8", abs(gunn_nu(0.4, 0.0, 0.7) - 3.8), 1e-14)

    # Schumann: closed form against an independent finite-difference solution
    for xi, tau in ((2.0, 1.0), (5.0, 5.0), (10.0, 5.0)):
        g_fd, s_fd = schumann_fd(xi, tau, 400)
        report("Schumann gas xi %g tau %g vs finite differences" % (xi, tau),
               abs(schumann_gas(xi, tau) - g_fd), 1e-4)
        report("Schumann solid xi %g tau %g vs finite differences" % (xi, tau),
               abs(schumann_solid(xi, tau) - s_fd), 1e-4)
    # equilibrium limit: far into the bed the gas front (theta_g = 0.5) sits at
    # tau = xi, i.e. it moves at the local-equilibrium speed w
    lo, hi = 150.0, 250.0
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if schumann_gas(200.0, mid, 20000) > 0.5:
            hi = mid
        else:
            lo = mid
    report("Schumann front at xi = 200: tau_50/xi - 1 (equilibrium speed)",
           abs(0.5 * (lo + hi) / 200.0 - 1.0), 1e-2)

    # Graetz plates against Shah and London
    report("Graetz plates Nu_Dh vs Shah and London 7.541",
           abs(graetz_plates_nu() - NU_GRAETZ_PLATES), 1e-3)
    # slug flow with wall flux: fully developed 12; local series against marching
    report("slug flow, wall flux, fully developed Nu 12",
           abs(slug_plates_nu_local(1e-9) - 12.0), 1e-9)
    for zeta in (0.02, 0.1, 0.5):
        # Gz = (D_h / z) Pe = 16 b^2 U / (a z) = 16 / zeta
        nu_s = slug_plates_nu_local(16.0 / zeta)
        nu_f = slug_plates_fd(zeta)
        report("slug flow local Nu, a z/(U b^2) = %g, vs marching" % zeta,
               abs(nu_s / nu_f - 1), 5e-3)
    # an air gap in series lowers Nu exactly as the resistance sum says
    nu0 = slug_plates_nu_local(1e-9)
    nug = slug_plates_nu_local(1e-9, 0.1)
    report("slug flow with gap: 1/Nu adds the gap term", abs(1 / nug - 1 / nu0 - 0.1), 1e-12)

    # mean Nu (eq. 9) is the length average of 1/Nu_local (eq. 4), with the gap
    for pe, g in ((800.0, 0.0), (300.0, 0.15)):
        dh_l = 0.02
        m = 4000
        inv = 0.0
        for i in range(m):
            z_over_dh = (i + 0.5) / m / dh_l
            inv += 1.0 / slug_plates_nu_local(pe / z_over_dh, g)
        report("slug flow mean Nu (Pe %g, gap %g) is the mean of 1/Nu_local" % (pe, g),
               abs(slug_plates_nu_mean(dh_l, pe, g) * inv / m - 1), 1e-4)

    print("thermal_reference --check: %s" % ("PASS" if ok else "FAIL"))
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
