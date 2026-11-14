#!/usr/bin/env python3
"""Reference values for the 2D silo benchmark T0 (WP 3.0).

Prints the Janssen and Beverloo numbers tabulated in docs/spec/silo-reference.md.
Pressures here are COMPRESSION POSITIVE (the solver convention is tension positive).

Usage:
    python3 tools/silo_reference.py                  # T0 values, markdown tables
    python3 tools/silo_reference.py --W 0.4 --H 1.2  # any inputs
    python3 tools/silo_reference.py --check          # internal consistency checks

Formulas (derivations in docs/spec/silo-reference.md):
    Janssen, plane bin of width W between two walls (perimeter/area = 2/W):
        d(sv)/dz = rho g - (2 K mu_w / W) sv
        sv(z)   = rho g lam (1 - exp(-z/lam)),   lam = W / (2 K mu_w)
        sh = K sv (wall normal),  tau_w = mu_w K sv (wall shear)
        90 percent of saturation at z90 = lam ln 10
    Beverloo/Hagen, 2D slot, continuum (no grain size, k = 0):
        Wdot = C rho sqrt(g) D^(3/2)    [kg/s per metre depth]
"""

import argparse
import math
import sys


def k_jaky(phi):
    """At-rest coefficient, Jaky: 1 - sin(phi)."""
    return 1.0 - math.sin(phi)


def k_elastic(nu):
    """Lateral-to-vertical ratio for elastic uniaxial strain: nu/(1-nu)."""
    return nu / (1.0 - nu)


def k_rankine_active(phi):
    """Rankine active coefficient (smooth wall, c = 0)."""
    s = math.sin(phi)
    return (1.0 - s) / (1.0 + s)


def k_wall_yield(phi, delta_w):
    """Active wall-yield ratio sigma_h/sigma_v at a wall where both the material
    (Mohr-Coulomb, c = 0, friction phi) and the wall (friction delta_w) slip.

    Mohr circle tangent to the yield line, passing through the wall point
    (sigma_w, tau_w = mu_w sigma_w); sigma_v is the other end of the diameter
    through that point, taking the major principal stress near vertical
    (active case). Closed form (Walker/Nedderman type), derived in the doc.
    """
    s = math.sin(phi)
    c2 = math.cos(phi) ** 2
    mu = math.tan(delta_w)
    root = math.sqrt(s * s - mu * mu * c2)
    return (1.0 + s * s - 2.0 * root) / (4.0 * mu * mu + c2)


def k_wall_yield_numeric(phi, delta_w):
    """Same as k_wall_yield by direct Mohr-circle construction (check)."""
    s = math.sin(phi)
    mu = math.tan(delta_w)
    # circle centre p = 1, radius R = s. Wall point at angle 2w from the
    # sigma axis measured from the minor end: sigma = 1 - s cos(2w),
    # tau = s sin(2w). Solve tau = mu sigma for the smaller 2w by bisection.
    f = lambda a: s * math.sin(a) - mu * (1.0 - s * math.cos(a))
    lo, hi = 0.0, math.pi / 2.0 + math.asin(min(1.0, mu))
    # f(0) < 0; find first sign change
    n = 2000
    a_prev = lo
    for i in range(1, n + 1):
        a = lo + (hi - lo) * i / n
        if f(a) >= 0.0:
            lo, hi = a_prev, a
            break
        a_prev = a
    else:
        raise RuntimeError("wall friction exceeds material friction")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    a = 0.5 * (lo + hi)
    sigma_h = 1.0 - s * math.cos(a)
    sigma_v = 1.0 + s * math.cos(a)
    return sigma_h / sigma_v


def janssen(z, rho, g, W, K, mu_w):
    lam = W / (2.0 * K * mu_w)
    sv = rho * g * lam * (1.0 - math.exp(-z / lam))
    return sv, K * sv, mu_w * K * sv


def janssen_rows(rho, g, W, K, mu_w):
    lam = W / (2.0 * K * mu_w)
    return {
        "lambda": lam,
        "s_inf": rho * g * lam,
        "z90": lam * math.log(10.0),
    }


def beverloo(C, rho, g, D):
    return C * rho * math.sqrt(g) * D ** 1.5


def k_table(phi, nu, delta_w):
    return [
        ("Elastic uniaxial strain, nu/(1-nu)", k_elastic(nu)),
        ("Jaky at rest, 1 - sin(phi)", k_jaky(phi)),
        ("Wall yield (active, wall slipping)", k_wall_yield(phi, delta_w)),
        ("Rankine active, (1-sin)/(1+sin)", k_rankine_active(phi)),
    ]


def print_tables(a):
    phi = math.radians(a.phi)
    dw = math.radians(a.delta_w)
    mu_w = math.tan(dw)
    g = a.g
    print("<!-- BEGIN silo_reference.py output -->")
    print(f"Inputs: rho = {a.rho:g} kg/m3, g = {g:g} m/s2, W = {a.W:g} m, "
          f"H = {a.H:g} m, phi = {a.phi:g} deg, delta_w = {a.delta_w:g} deg "
          f"(mu_w = {mu_w:.4f}), nu = {a.nu:g}")
    print()
    print("| Stress ratio K | K | lambda (m) | sigma_v,sat (Pa) | z90 (m) "
          "| mean sigma_v at H (Pa) | sigma_v(H) / rho g H |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for name, K in k_table(phi, a.nu, dw):
        r = janssen_rows(a.rho, g, a.W, K, mu_w)
        sv, _, _ = janssen(a.H, a.rho, g, a.W, K, mu_w)
        print(f"| {name} | {K:.4f} | {r['lambda']:.4f} | {r['s_inf']:.0f} "
              f"| {r['z90']:.3f} | {sv:.0f} | {sv / (a.rho * g * a.H):.4f} |")
    print()
    K = k_elastic(a.nu)
    print(f"Profile with K = {K:.4f} (elastic), wall normal = K sigma_v, "
          f"wall shear = mu_w K sigma_v:")
    print()
    print("| z (m) | rho g z (Pa) | sigma_v (Pa) | sigma_h wall (Pa) | tau_w (Pa) |")
    print("| --- | --- | --- | --- | --- |")
    nz = 6
    for i in range(nz + 1):
        z = a.H * i / nz
        sv, sh, tw = janssen(z, a.rho, g, a.W, K, mu_w)
        print(f"| {z:.3f} | {a.rho * g * z:.0f} | {sv:.0f} | {sh:.0f} | {tw:.0f} |")
    print()
    print("| D (m) | D/W | Wdot, C = 1.40 (kg/s/m) | Wdot, C = 1.45 (kg/s/m) "
          "| Wdot, C = 1.55 (kg/s/m) | emptying time at C = 1.45 (s) |")
    print("| --- | --- | --- | --- | --- | --- |")
    m0 = a.rho * a.W * a.H
    for D in a.D:
        w = [beverloo(C, a.rho, g, D) for C in (1.40, 1.45, 1.55)]
        print(f"| {D:.3f} | {D / a.W:.3f} | {w[0]:.1f} | {w[1]:.1f} | {w[2]:.1f} "
              f"| {m0 / w[1]:.2f} |")
    print("<!-- END silo_reference.py output -->")


def check(a):
    """Self-consistency checks; return number of failures."""
    fails = 0
    phi = math.radians(a.phi)
    for dw_deg in (5.0, 10.0, 20.0, 30.0):
        dw = math.radians(dw_deg)
        k1 = k_wall_yield(phi, dw)
        k2 = k_wall_yield_numeric(phi, dw)
        ok = abs(k1 - k2) < 1e-9
        fails += not ok
        print(f"{'ok ' if ok else 'BAD'} wall-yield K closed form {k1:.10f} "
              f"vs Mohr circle {k2:.10f} (delta_w {dw_deg:g})")
    # delta_w -> 0 gives Rankine active
    k0 = k_wall_yield(phi, 1e-9)
    ok = abs(k0 - k_rankine_active(phi)) < 1e-9
    fails += not ok
    print(f"{'ok ' if ok else 'BAD'} wall-yield K at delta_w = 0 equals Rankine "
          f"active ({k0:.10f})")
    # Janssen ODE by RK4 against the closed form
    dw = math.radians(a.delta_w)
    mu_w = math.tan(dw)
    K = k_elastic(a.nu)
    rho, g, W = a.rho, a.g, a.W
    f = lambda sv: rho * g - 2.0 * K * mu_w / W * sv
    n = 4000
    h = 3.0 * a.H / n
    sv = 0.0
    for _ in range(n):
        k1 = f(sv)
        k2 = f(sv + 0.5 * h * k1)
        k3 = f(sv + 0.5 * h * k2)
        k4 = f(sv + h * k3)
        sv += h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    ref, _, _ = janssen(3.0 * a.H, rho, g, W, K, mu_w)
    ok = abs(sv - ref) / ref < 1e-10
    fails += not ok
    print(f"{'ok ' if ok else 'BAD'} Janssen ODE (RK4) {sv:.6f} vs closed form "
          f"{ref:.6f}")
    # force balance on the column of height Z: weight = base force + 2 walls
    Z = 3.0 * a.H
    lam = W / (2.0 * K * mu_w)
    base = janssen(Z, rho, g, W, K, mu_w)[0] * W
    # integral of tau_w over 0..Z, both walls
    wall = 2.0 * mu_w * K * rho * g * lam * (Z - lam * (1.0 - math.exp(-Z / lam)))
    weight = rho * g * W * Z
    ok = abs(weight - base - wall) / weight < 1e-12
    fails += not ok
    print(f"{'ok ' if ok else 'BAD'} vertical force balance: weight {weight:.6f} "
          f"= base {base:.6f} + walls {wall:.6f}")
    # Beverloo exponent: ratio of rates for D and 2D is 2^1.5
    r = beverloo(1.45, rho, g, 0.08) / beverloo(1.45, rho, g, 0.04)
    ok = abs(r - 2.0 ** 1.5) < 1e-12
    fails += not ok
    print(f"{'ok ' if ok else 'BAD'} Beverloo D^1.5 scaling ratio {r:.12f}")
    return fails


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rho", type=float, default=1800.0, help="bulk density kg/m3")
    p.add_argument("--g", type=float, default=9.81)
    p.add_argument("--W", type=float, default=0.40, help="tank width m")
    p.add_argument("--H", type=float, default=0.60, help="fill height m")
    p.add_argument("--phi", type=float, default=35.0, help="friction angle deg")
    p.add_argument("--delta_w", type=float, default=20.0, help="wall friction deg")
    p.add_argument("--nu", type=float, default=0.3)
    p.add_argument("--D", type=float, nargs="+", default=[0.04, 0.06, 0.08])
    p.add_argument("--check", action="store_true")
    a = p.parse_args(argv)
    if a.check:
        return 1 if check(a) else 0
    print_tables(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
