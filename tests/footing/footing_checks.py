"""WP 4.4 done-when check: bearing capacity of a strip footing on undrained
clay (Prandtl, Nc = 2 + pi = 5.14), the footing a moving STL body (decision
0053). Run by tests/test_footing.sh after celVoidCheck is built:

    python3 footing_checks.py WORKDIR

Weightless Tresca clay (Mohr-Coulomb phi 0, c = cu = 10 kPa, E 3 MPa, nu 0.3,
rho 1800), 0.6 m wide and 0.3 m deep in a box with frictionless mesh walls,
void above, 10 mm cells. A smooth rigid strip footing B = 0.1 m (an STL box,
sides cutting cells at a quarter, bottom 2.5 mm below the surface) is pushed
down, its speed ramped from 0 to 0.025 m/s over 0.05 s (quasi-static: 0.05
m/s gives 4 % more), to a settlement of about 0.19 B. The footing pressure
q = Fy / (B depth) against settlement s must reach a limit: q / cu averaged
over s = 0.10 B to 0.18 B within 10 % of 5.14, and rising by less than 5 %
across that range.

    python3 footing_checks.py WORKDIR --fine   # 5 mm cells (manual, 9 min) Mass kept to 1e-10, no
mass lost in hand-overs.

Prints PASS/FAIL, exits non-zero on failure.
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("silo", "void", "lagrange", "stlwalls"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import silo_checks as sc  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_push_checks as sp  # noqa: E402

CELL = 0.01
B = 0.1
CU = 10e3
E_CLAY = 3e6
NX, NY = 60, 34                # 0.6 m wide, 0.3 m of clay plus 4 void rows
DEPTH = 0.3
X0 = 0.25 + 0.0025             # footing x0 .. x0 + B (quarter-cell cuts)
EMBED = 0.0025
SPEED, T_RAMP, T_END = 0.025, 0.05, 0.775
NC = 2 + math.pi
TOL = 0.10
WINDOW = (0.10, 0.18)          # s / B


def clay():
    return f"type mohrCoulomb;\nE {E_CLAY!r};\nnu 0.3;\nphi 0;\nc {CU!r};\npsi 0;\n"


def run_footing(work):
    p = vc.make_case(f"{work}/footing", NX, NY, clay(),
                     {"fill": [-1.0, -1.0, 1.0, DEPTH], "U0": [0, 0, 0],
                      "tEnd": T_END, "maxCo": 0.4}, gravity=False)
    v = SPEED
    table = (f"table ((0 (0 0 0)) ({T_RAMP!r} (0 {-v!r} 0)) ({T_END + 1!r} (0 {-v!r} 0)))")
    y0 = DEPTH - EMBED
    sp.add_moving_body(p, "footing", (X0, y0, -0.01), (X0 + B, y0 + 0.1, 0.02), table, 0.0)
    t = vc.run(p)
    return p, t


def settlement(t):
    """Footing displacement at times t (ramp then constant speed)."""
    t = np.asarray(t)
    return np.where(t < T_RAMP, 0.5 * SPEED * t ** 2 / T_RAMP,
                    0.5 * SPEED * T_RAMP + SPEED * (t - T_RAMP))


def check_footing(work):
    print("[1] Smooth strip footing on weightless Tresca clay (cu 10 kPa, E/cu 300),")
    print(f"    B = 0.1 m = {B / CELL:.0f} cells, pushed to s = 0.19 B")
    p, sec = run_footing(work)
    h = np.loadtxt(f"{p}/history.dat")
    s = settlement(h[:, sp.T]) / B
    q = -h[:, sp.FYB] / (B * CELL) / CU      # force of the body on the material is down
    w = (s >= WINDOW[0]) & (s <= WINDOW[1])
    qm = q[w].mean()
    a = np.polyfit(s[w], q[w], 1)
    rise = a[0] * (WINDOW[1] - WINDOW[0]) / qm
    mass = h[0, sp.MASS]
    dm = abs(h[-1, sp.MASS] / mass - 1)
    lost = h[:, sp.LOST].max() / mass
    err = qm / NC - 1
    ok = abs(err) <= TOL and abs(rise) <= 0.05 and dm <= 1e-10 and lost == 0
    det = [f"q / cu over s = {WINDOW[0]} B to {WINDOW[1]} B: {qm:.3f} (Prandtl 5.142, "
           f"{100 * err:+.1f} %, limit {100 * TOL:.0f} %), change across it {100 * rise:+.1f} % "
           f"(<= 5 %)"]
    for sb in (0.01, 0.02, 0.05, 0.1, 0.15, 0.18):
        i = int(np.argmin(np.abs(s - sb)))
        det.append(f"  s = {s[i]:.3f} B: q / cu = {q[i]:.3f}")
    det.append(f"mass change {dm:.0e} (<= 1e-10), lost in hand-overs {lost:.0e}; "
               f"{int(h[-1, sp.STEP])} steps, {sec:.0f} s")
    np.savetxt(f"{p}/q.dat", np.c_[s, q], header="s/B q/cu")
    return sc.report("Strip footing bearing capacity", ok, det)


def main():
    work = sys.argv[1]
    if "--fine" in sys.argv[2:]:
        global CELL, NX, NY
        vc.CELL = CELL = 0.005
        NX, NY = 120, 68
    t0 = time.perf_counter()
    ok = check_footing(work)
    print(f"\nWP 4.4 checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
