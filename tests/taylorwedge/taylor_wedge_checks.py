"""WP 2.6 done-when: the Taylor bar on the axisymmetric wedge (decisions 0061,
0063), the Eulerian counterpart of the Lagrangian WP 5.2 (decision 0071).

Copper rod of the Abaqus Benchmarks ("Impact of a copper rod", after
Ferencz 1989): 32.4 mm long, radius 3.2 mm, 227 m/s onto a rigid
frictionless floor, E 110 GPa, nu 0.3, rho 8970, von Mises perfectly
plastic at 314 MPa, 80 microseconds, no gravity. Reference: the axisymmetric
CAX4R run of that benchmark, shortening 13.11 mm, widening (radial
displacement of the impact face edge) 6.006 mm (other elements 12.71 to
13.18 / 5.52 to 6.06).

Eulerian wedge (2.5 deg half angle) 12 x 36 mm, void around the rod,
`freeSurface yes`, maxCo 0.4, cells 0.4 mm and 0.2 mm.
Final length: material in the axis column (sum of alpha times the cell
height). Final face radius: material in the bottom cell row, area weighted
(r^2 = sum alpha (r_out^2 - r_in^2)).

Checks:
1. Mass kept to 1e-10 at both sizes.
2. 0.2 mm: shortening within 5 % of 13.11 mm.
3. 0.2 mm: widening within 20 % of 6.006 mm (the thin lip at the face edge
   is smeared over a cell or two: decision 0063).
4. Both errors smaller at 0.2 mm than at 0.4 mm.
4a. The widening extrapolated to zero cell size from the two sizes (first
   order, decision 0098) within 10 % of 6.006 mm.
5. Kinetic energy at 80 microseconds below 5 % of the initial.

Run by tests/test_taylor_wedge.sh after celVoidCheck is built:

    python3 taylor_wedge_checks.py WORKDIR
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "wedgevoid"))
import wedge_void_checks as wv  # noqa: E402

vc = wv.vc
L0, R0, V0, RHO = 0.0324, 0.0032, 227.0, 8970.0
SHORT_REF, WIDE_REF = 13.11e-3, 6.006e-3
R, H, T_END, MAX_CO = 0.012, 0.036, 80e-6, 0.4
SIZES = (0.0004, 0.0002)
MATERIAL = "type vonMises;\nE 110e9;\nnu 0.3;\nsigma_y 314e6;\n"
MASS, KE = 2, 12        # history columns


def taylor(work, h):
    wv.CELL = vc.CELL = h
    check = {"fill": [0, 0, R0 * wv.COS_T, L0], "U0": [0, -V0, 0],
             "tEnd": T_END, "maxCo": MAX_CO}
    p = wv.wedge_void_case(f"{work}/taylor_{h * 1e3:g}mm", R, H, MATERIAL, check,
                           gravity=False)
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read().replace(f"rho {vc.RHO};", f"rho {RHO};")
    open(mp, "w").write(txt)
    sec = vc.run(p)
    c = np.loadtxt(f"{p}/cells.dat")
    r, y, a = c[:, 0] / wv.COS_T, c[:, 1], c[:, 2]
    axis = np.abs(r - r.min()) < 1e-9
    length = np.sum(a[axis]) * h
    bot = np.abs(y - y.min()) < 1e-9
    rf = math.sqrt(np.sum(a[bot] * ((r[bot] + h / 2)**2 - (r[bot] - h / 2)**2)))
    hist = np.loadtxt(f"{p}/history.dat")
    return {"sec": sec, "short": L0 - length, "wide": rf - R0,
            "mass": abs(hist[-1, MASS] - hist[0, MASS]) / hist[0, MASS],
            "ke": hist[-1, KE] / hist[0, KE]}


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    fails = []

    def check(ok, text):
        print(("PASS " if ok else "FAIL ") + text)
        if not ok:
            fails.append(text)

    res = {}
    for h in SIZES:
        s = res[h] = taylor(work, h)
        print(f"{h * 1e3:g} mm cells: {s['sec']:.0f} s; shortening {s['short'] * 1e3:.2f} mm "
              f"(13.11), widening {s['wide'] * 1e3:.2f} mm (6.006); kinetic energy left "
              f"{s['ke']:.3f}")
    for h, s in res.items():
        check(s["mass"] < 1e-10, f"{h * 1e3:g} mm: mass kept ({s['mass']:.1e})")
    f = res[SIZES[-1]]
    es = (f["short"] - SHORT_REF) / SHORT_REF
    ew = (f["wide"] - WIDE_REF) / WIDE_REF
    check(abs(es) <= 0.05, f"shortening within 5 % ({100 * es:+.1f} %)")
    check(abs(ew) <= 0.20, f"widening within 20 % ({100 * ew:+.1f} %)")
    c = res[SIZES[0]]
    check(abs(f["short"] - SHORT_REF) < abs(c["short"] - SHORT_REF)
          and abs(f["wide"] - WIDE_REF) < abs(c["wide"] - WIDE_REF),
          "both errors smaller at the finer size")
    check(f["ke"] < 0.05, f"kinetic energy left below 5 % ({f['ke']:.3f})")
    # A.7 (decision 0098): the lip error is first order in the cell size, so the
    # two sizes extrapolate (Richardson, h -> 0) to the converged widening
    wx = 2 * f["wide"] - c["wide"]
    sx = 2 * f["short"] - c["short"]
    ewx = (wx - WIDE_REF) / WIDE_REF
    print(f"extrapolated to h -> 0: widening {wx * 1e3:.2f} mm ({100 * ewx:+.1f} %), "
          f"shortening {sx * 1e3:.2f} mm ({100 * (sx - SHORT_REF) / SHORT_REF:+.1f} %)")
    check(abs(ewx) <= 0.10, f"extrapolated widening within 10 % ({100 * ewx:+.1f} %)")
    print(f"WP 2.6 checks: {'ALL PASS' if not fails else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
