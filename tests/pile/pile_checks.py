"""WP 4.4b done-when check: limiting lateral pressure on a pile in undrained
clay (Randolph and Houlsby 1984, "The limiting pressure on a circular pile
loaded laterally in cohesive soil", Geotechnique 34(4), 613-623): a plane
section of a deep pile moving sideways through Tresca clay, full flow
around it, P / (cu D) = 9.14 for a smooth pile (11.94 fully rough). The pile
is a moving STL body (decision 0053). Run by tests/test_pile.sh after
celVoidCheck is built:

    python3 pile_checks.py WORKDIR [--fine]

Weightless Tresca clay (Mohr-Coulomb phi 0, cu 10 kPa, E 3 MPa, nu 0.3,
rho 1800) fills a 0.6 x 0.5 m box with frictionless mesh walls, under a
uniform confining pressure of 10 cu (so the clay closes behind the pile, as
in the full-flow mechanism of the reference, instead of leaving a gap), with
`fullBlend yes` (decision 0049): without it, cells left a trace of void by
the moving pile take the partly filled rules and the force is noisy and far
too low (decision 0056). The
pile, D = 0.1 m (10 cells) as a 64-sided prism, smooth, moves along x, its
speed ramped from 0 to 0.025 m/s over 0.05 s, to a travel of 0.3 D. The
force per length P = Fx / depth must reach a limit: P / (cu D) averaged over
a travel of 0.15 D to 0.3 D within 10 % of 9.14 and changing by less than 5 %
across that range. Mass kept to 1e-10, no mass lost in hand-overs.

Prints PASS/FAIL, exits non-zero on failure.
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("silo", "void", "lagrange", "stlwalls", "outlet"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import silo_checks as sc  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_push_checks as sp  # noqa: E402
import stl_silo_checks as ss  # noqa: E402

CELL = 0.01
D = 0.1
CU = 10e3
E_CLAY = 3e6
P0 = 10 * CU
NX, NY = 60, 50
CENTRE = (0.2, 0.2525)
NSIDES = 64
SPEED, T_RAMP = 0.025, 0.05
TRAVEL = 0.3 * D
NP = 9.14                       # smooth pile, Randolph and Houlsby (1984)
TOL = 0.10
WINDOW = (0.15, 0.3)            # travel / D


def t_end():
    return T_RAMP + TRAVEL / SPEED - 0.5 * T_RAMP


def travel(t):
    t = np.asarray(t)
    return np.where(t < T_RAMP, 0.5 * SPEED * t ** 2 / T_RAMP,
                    0.5 * SPEED * T_RAMP + SPEED * (t - T_RAMP))


def add_pile(path):
    cx, cy = CENTRE
    # circumscribed polygon of the same area as the circle
    r = 0.5 * D * math.sqrt((math.pi / NSIDES) / math.sin(math.pi / NSIDES)
                            / math.cos(math.pi / NSIDES))
    poly = [(cx + r * math.cos(2 * math.pi * (k + 0.5) / NSIDES),
             cy + r * math.sin(2 * math.pi * (k + 0.5) / NSIDES)) for k in range(NSIDES)]
    cs.write(f"{path}/constant/triSurface/pile.stl", ss.prism_stl("pile", poly))
    v, te = SPEED, t_end()
    table = f"table ((0 (0 0 0)) ({T_RAMP!r} ({v!r} 0 0)) ({te + 1!r} ({v!r} 0 0)))"
    body = (f"    pile\n    {{\n        surface pile.stl;\n        friction 0;\n"
            f"        velocity {table};\n    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def check_pile(work):
    print(f"[1] Smooth pile section D = 0.1 m ({D / CELL:.0f} cells) moving sideways through")
    print("    weightless Tresca clay (cu 10 kPa, E/cu 300, confined at 10 cu), travel 0.3 D")
    mat = (f"type mohrCoulomb;\nE {E_CLAY!r};\nnu 0.3;\nphi 0;\nc {CU!r};\npsi 0;\n"
           "fullBlend yes;\n")
    p = vc.make_case(f"{work}/pile", NX, NY, mat,
                     {"fill": [-1.0, -1.0, 1.0, 1.0], "U0": [0, 0, 0],
                      "sigma0": [-P0, 0, 0, -P0, 0, -P0],
                      "tEnd": t_end(), "maxCo": 0.4}, gravity=False)
    add_pile(p)
    sec = vc.run(p)
    h = np.loadtxt(f"{p}/history.dat")
    u = travel(h[:, sp.T]) / D
    q = h[:, sp.FXB] / CELL / (CU * D)      # force of the pile on the clay, along +x
    w = (u >= WINDOW[0]) & (u <= WINDOW[1])
    qm = q[w].mean()
    a = np.polyfit(u[w], q[w], 1)
    rise = a[0] * (WINDOW[1] - WINDOW[0]) / qm
    mass = h[0, sp.MASS]
    dm = abs(h[-1, sp.MASS] / mass - 1)
    lost = h[:, sp.LOST].max() / mass
    err = qm / NP - 1
    ok = abs(err) <= TOL and abs(rise) <= 0.05 and dm <= 1e-10 and lost == 0
    det = [f"P / (cu D) over a travel of {WINDOW[0]} D to {WINDOW[1]} D: {qm:.3f} "
           f"(Randolph and Houlsby smooth 9.14, {100 * err:+.1f} %, limit {100 * TOL:.0f} %), "
           f"change across it {100 * rise:+.1f} % (<= 5 %)"]
    for ub in (0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3):
        i = int(np.argmin(np.abs(u - ub)))
        det.append(f"  travel {u[i]:.3f} D: P / (cu D) = {q[i]:.3f}")
    det.append(f"mass change {dm:.0e} (<= 1e-10), lost in hand-overs {lost:.0e}; "
               f"{int(h[-1, sp.STEP])} steps, {sec:.0f} s")
    np.savetxt(f"{p}/p.dat", np.c_[u, q], header="travel/D P/(cu D)")
    return sc.report("Pile section moving through clay", ok, det)


def main():
    global CELL, NX, NY
    work = sys.argv[1]
    if "--fine" in sys.argv[2:]:
        vc.CELL = CELL = 0.005
        NX, NY = 120, 100
    t0 = time.perf_counter()
    ok = check_pile(work)
    print(f"\nWP 4.4b checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
