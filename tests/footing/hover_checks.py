"""Body hovering above a free surface (decision 0091): a rigid footing held
2.5 mm above a sand layer at rest must leave the layer at rest.

Layer: cohesionless Mohr-Coulomb phi 30, E 10 MPa, nu 0.3, rho 1800, at its
K0 stress, 0.6 m wide and 0.3 m deep, 10 mm cells (as tests/footing/
sand_checks.py). Footing B = 0.1 m, base 2.5 mm above the surface, held still,
0.2 s. Not passing yet (2026-11-19): with psi 0 the top cells at the footing
edges reach 0.1 m/s with stresses about four times the geostatic value at the
footing's depth; with psi 30 (associated) the run blows up near 0.12 s. The
same layer without the footing stays at rest to 1e-17 m/s.

Checks, for psi 0 and psi 30:
1. the run reaches 0.2 s;
2. material cells (alpha > 0.5) move slower than 1e-3 m/s at the end.

    python3 hover_checks.py WORKDIR
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("silo", "void", "lagrange", "stlwalls"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import void_checks as vc  # noqa: E402
import stl_push_checks as sp  # noqa: E402

DEPTH, GAP, B, X0 = 0.3, 0.0025, 0.1, 0.2525
T_END = 0.2
U_TOL = 1e-3


def case(work, psi):
    mat = (f"type mohrCoulomb;\nE 10e6;\nnu 0.3;\nphi 30;\nc 0;\npsi {psi!r};\n")
    p = vc.make_case(f"{work}/psi{psi:g}", 60, 34, mat,
                     {"fill": [-1.0, -1.0, 1.0, DEPTH], "U0": [0, 0, 0],
                      "initK0": 0.5, "tEnd": T_END, "maxCo": 0.4})
    y0 = DEPTH + GAP
    sp.add_moving_body(p, "footing", (X0, y0, -0.01), (X0 + B, y0 + 0.1, 0.02),
                       "table ((0 (0 0 0)) (10 (0 0 0)))", 0.0)
    try:
        vc.run(p)
    except Exception:  # noqa: BLE001  (a blown-up run is a failed check)
        pass
    h = np.genfromtxt(f"{p}/history.dat", usecols=range(10), invalid_raise=False)
    t_end = float(h[-1, 1])
    u = float("inf")
    if os.path.exists(f"{p}/cells.dat"):
        c = np.loadtxt(f"{p}/cells.dat")
        full = c[:, 2] > 0.5
        u = float(np.hypot(c[full, 4], c[full, 5]).max())
    return t_end, u


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    fails = 0
    for psi in (0.0, 30.0):
        t_end, u = case(work, psi)
        ok = t_end >= T_END - 1e-9 and u < U_TOL
        fails += not ok
        print(f"{'PASS' if ok else 'FAIL'} psi {psi:g}: reached t = {t_end:.3f} s of {T_END}; "
              f"max material speed {u:.2e} m/s (< {U_TOL:g})")
    print(f"Hover checks: {'ALL PASS' if not fails else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
