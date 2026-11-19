"""WP H.2 done-when checks: the conical hopper wall on the axisymmetric wedge
(decision 0101).

The Couto et al. (2013) test silo (docs/spec/conical-hopper.md, case A) on a 5 deg
wedge, 25 mm cells: cylinder wall = mesh patch, 34.3 deg hopper = STL cut-cell
wall extruded in the wedge's face planes, outlet closed. Wheat, Mohr-Coulomb
(phi 30.17, psi 1, E 10.674 MPa), wall friction 0.2, filled level to the top,
gravity ramped over 0.1 s with damping 100 1/s, run to 0.5 s.

1. Mass: the initial mass is the silo volume times rho (0.2 %): the cut-cell
   hopper holds the right open volume.
2. Closed outlet holds: nothing discharged, mass constant to 1e-10.
3. Settled: largest speed below 0.05 m/s at 0.5 s.
4. Cylinder wall pressure at the load cells GA1 to GA7 within 10 % of circular
   Janssen with K = nu/(1 - nu) (the at-rest ratio an elastic fill under gravity
   reaches; the interior ratio is checked within 0.03 of it).
5. Mean vertical stress at the transition within the experiment's static range
   (11.65 to 13.35 kPa, P2 Table 1) widened by 15 %.
Reported, not checked (H.3): the hopper wall pressure against GA0.

Run by tests/test_hopper_wedge.sh after celVoidCheck is built:

    python3 hopper_wedge_checks.py WORKDIR
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hopper_case as hc  # noqa: E402
import void_checks as vc  # noqa: E402

CELL = float(os.environ.get("CEL_CELL", "0.025"))
T_END, DAMPING, T_RAMP = 0.5, 100.0, 0.1
G = 9.81
CELLS = {"GA7": 2.075, "GA6": 1.825, "GA5": 1.575, "GA4": 1.325, "GA3": 1.075,
         "GA2": 0.825, "GA1": 0.575}
GA0_STATIC = (8.61, 9.61)
PVT_STATIC = (11.65, 13.35)


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    fails = []

    def check(ok, text):
        print(f"  {text}: {'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(text)

    p, c = hc.couto_case(work, "couto_mc", "mc", CELL, dict(tEnd=T_END, maxCo=0.4))
    open(f"{p}/system/controlDict", "a").write(
        f"gravityRamp {T_RAMP!r};\ndamping {DAMPING!r};\n")
    wall = vc.run(p)
    h = np.loadtxt(f"{p}/history.dat")
    b, ht = hc.hopper_geometry(c)
    H = ht + c["Hc"]
    R, Ro = c["R"], c["Rout"]
    t = math.radians(hc.THETA)
    scale = math.pi / (math.sin(t) * math.cos(t))      # wedge to full circle
    V = math.pi * R**2 * c["Hc"] + math.pi * ht / 3 * (R**2 + R * Ro + Ro**2)
    m0 = h[0, 2] * scale
    print(f"[1] Couto silo at rest on the wedge, {CELL * 1e3:.1f} mm cells "
          f"({wall:.0f} s)")
    check(abs(m0 / (c["rho"] * V) - 1) < 2e-3,
          f"mass {m0:.1f} kg against {c['rho'] * V:.1f} kg")
    check(abs(h[:, 14]).max() == 0 and abs(h[:, 2] - h[0, 2]).max() / h[0, 2] < 1e-10,
          f"closed outlet: discharged {abs(h[:, 14]).max():.1e} kg, mass drift "
          f"{abs(h[:, 2] - h[0, 2]).max() / h[0, 2]:.1e}")
    check(h[-1, 9] < 0.05, f"settled: max |U| {h[-1, 9]:.4f} m/s at {h[-1, 1]:.2f} s")

    a = np.loadtxt(f"{p}/cells.dat")
    x, y, sxx, syy = a[:, 0], a[:, 1], a[:, 6], a[:, 7]
    K = c["nu"] / (1 - c["nu"])
    inner = (y > 1.0) & (y < 1.5) & (x < 0.3)
    Ki = float(np.mean(sxx[inner] / syy[inner]))
    check(abs(Ki - K) < 0.03, f"interior stress ratio {Ki:.3f} against nu/(1-nu) {K:.3f}")
    cyl, hop = hc.wall_pressures(p, c, CELL)
    mu = c["mu"]
    for name, z in CELLS.items():
        i = int(np.argmin(abs(cyl[:, 0] - z)))
        zc, pc = cyl[i]
        pj = c["rho"] * G * 2 * R / (4 * mu) * (1 - math.exp(-4 * mu * K * (H - zc) / (2 * R))) / 1e3
        check(abs(pc / pj - 1) <= 0.10,
              f"{name} z {zc:.3f} m: {pc:.2f} kPa, Janssen {pj:.2f} ({100 * (pc / pj - 1):+.1f} %)")
    sel = np.abs(y - (ht + 0.5 * CELL)) < 1e-3 * CELL + 0.25 * CELL
    pvt = -(syy[sel] * x[sel]).sum() / x[sel].sum() / 1e3
    check(PVT_STATIC[0] / 1.15 <= pvt <= PVT_STATIC[1] * 1.15,
          f"mean vertical stress at the transition {pvt:.2f} kPa (experiment "
          f"{PVT_STATIC[0]} to {PVT_STATIC[1]})")
    near = hop[np.abs(hop[:, 0] - 0.41) < 1.5 * CELL]
    print(f"  report: hopper wall at z 0.41 m {near[:, 1].mean():.2f} kPa (GA0 at rest "
          f"{GA0_STATIC[0]} to {GA0_STATIC[1]}; H.3); hopper wall "
          f"{hop[:, 1].min():.2f} to {hop[:, 1].max():.2f} kPa")
    print(f"{len(fails)} failures ({time.perf_counter() - t0:.0f} s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
