"""WP G.4 done-when checks: the Couto silo at rest on the wall-aligned wedge
mesh (fanned trapezoid hopper block, decision 0103) with the free surface.

25 mm aligned mesh (hopper_case.aligned_case: hopper wall = mesh patch, the
hopper columns fan from vertical at the axis to the wall slope), MC with the
T0 rules and fullBlend, outlet closed, filled level to the top, gravity
ramped over 0.1 s with damping 100 1/s, run to 0.5 s.

1. Mass constant to 1e-10 and nothing through the closed outlet.
2. Settled: largest speed below 0.05 m/s at 0.5 s.
3. Interior stress ratio within 0.03 of nu/(1 - nu).
4. Cylinder wall pressure at GA1 to GA7 within 5 % of circular Janssen with
   K = nu/(1 - nu) (as tests/hopper/hopper_wedge_checks.py, there within
   10 % on the cut-cell hopper).
Reported: mean vertical stress at the transition, hopper wall at GA0.

    python3 wedge_aligned_checks.py WORKDIR
"""
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hopper"))
import hopper_case as hc  # noqa: E402
import hopper_pressure_checks as hp  # noqa: E402

vc = hp.vc
CELL, T_END, G = 0.025, 0.5, 9.81
CELLS = {"GA7": 2.075, "GA6": 1.825, "GA5": 1.575, "GA4": 1.325, "GA3": 1.075,
         "GA2": 0.825, "GA1": 0.575}


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    fails = []

    def check(ok, text):
        print(f"  {text}: {'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            fails.append(text)

    os.environ.setdefault("CEL_HOPPER_RULES", hp.MU_I + hp.RULES + "fullBlend yes;\n")
    p, c = hc.aligned_case(work, "couto_aligned_rest", "mc", CELL,
                           dict(tEnd=T_END, maxCo=0.4), outlet_ratio=1 / 3)
    open(f"{p}/system/controlDict", "a").write("gravityRamp 0.1;\ndamping 100;\n")
    wall = vc.run(p)
    h = vc.history(p)
    print(f"[1] Couto silo at rest on the aligned wedge, {CELL * 1e3:.0f} mm ({wall:.0f} s)")
    drift = abs(h[:, vc.MASS] - h[0, vc.MASS]).max() / h[0, vc.MASS]
    check(drift < 1e-10 and abs(h[:, 14]).max() == 0,
          f"mass drift {drift:.1e}, through the outlet {abs(h[:, 14]).max():.1e} kg")
    check(h[-1, vc.UMAX] < 0.05, f"settled: max |U| {h[-1, vc.UMAX]:.4f} m/s at {h[-1, vc.T]:.2f} s")
    a = np.loadtxt(f"{p}/cells.dat")
    x, y, sxx, syy = a[:, 0], a[:, 1], a[:, 6], a[:, 7]
    K = c["nu"] / (1 - c["nu"])
    inner = (y > 1.0) & (y < 1.5) & (x < 0.3) & (a[:, 2] > 0.99)
    Ki = float(np.mean(sxx[inner] / syy[inner]))
    check(abs(Ki - K) < 0.03, f"interior stress ratio {Ki:.3f} against nu/(1-nu) {K:.3f}")
    b, ht = hc.hopper_geometry(c)
    H, R, mu = ht + c["Hc"], c["R"], c["mu"]
    cyl, hop = hc.wall_pressures(p, c, CELL)
    for name, z in CELLS.items():
        i = int(np.argmin(abs(cyl[:, 0] - z)))
        zc, pc = cyl[i]
        pj = c["rho"] * G * 2 * R / (4 * mu) * (1 - math.exp(-4 * mu * K * (H - zc) / (2 * R))) / 1e3
        check(abs(pc / pj - 1) <= 0.05,
              f"{name} z {zc:.3f} m: {pc:.2f} kPa, Janssen {pj:.2f} ({100 * (pc / pj - 1):+.1f} %)")
    near = hop[np.abs(hop[:, 0] - 0.41) < 1.5 * CELL]
    print(f"  report: hopper wall at z 0.41 m {near[:, 1].mean():.2f} kPa (GA0 8.61 to 9.61)")
    print(f"{len(fails)} failures ({time.perf_counter() - t0:.0f} s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
