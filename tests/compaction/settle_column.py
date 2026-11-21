"""WP K.2 done-when checks: the compaction cap in the Eulerian solver
(decision 0113). A column settling under its own weight. About 1 minute.

A 1 m column (one 20 mm cell wide, 50 cells, frictionless walls, 5 void
cells above) of Mohr-Coulomb (phi 30, E 10 MPa, nu 0.3, rho0 1800) with the
cap lambda* = 0.01, capP0 100 Pa, starts stress-free; gravity ramps up over
2 s with damping 100 1/s and the run ends at 3 s. fullBlend yes: with the
plain partly filled rule the stress-free surface cell jolts the column each
time the settling surface crosses a cell face, and the cap keeps every
overshoot (cells 5 to 18 % over the line near the top and the base).

[1] Every material cell lies on the compression line it was loaded along:
    ln(rho/rho0) = P0/K + lambda* ln(P/P0) for P > P0 (elastic to P0, then
    on the line) within 2 % of the strain, with P from the cell's stress;
    and p_c = P in loaded cells (normally consolidated, at rest).
[2] Mass constant to 1e-10; the column at rest (largest speed in full cells below
    0.1 mm/s); the top settles (more than 20 mm).
[3] Reported: the same column without the cap (elastic compression only).

    python3 settle_column.py WORKDIR
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import void_checks as vc  # noqa: E402

LAM, P0, E, NU = 0.01, 100.0, 10.0e6, 0.3
K = E / (3 * (1 - 2 * NU))
H, CELL, NFILL, NVOID = 1.0, 0.02, 50, 5
fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def case(work, cap):
    vc.CELL = CELL
    name = "cap" if cap else "nocap"
    mat = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi 30;\nc 0;\npsi 0;\n")
    if cap:
        mat += f"capLambda {LAM!r};\ncapP0 {P0!r};\n"
    mat += "fullBlend yes;\n"
    check = dict(fill=[0.0, 0.0, CELL, H], U0=[0.0, 0.0, 0.0], tEnd=3.0, maxCo=0.4)
    p = vc.make_case(f"{work}/{name}", 1, NFILL + NVOID, mat, check)
    cd = f"{p}/system/controlDict"
    open(cd, "a").write("gravityRamp 2;\ndamping 100;\n")
    vc.run(p)
    c = np.loadtxt(f"{p}/cells.dat")
    h = np.loadtxt(f"{p}/history.dat")
    return c, h


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    c, h = case(work, True)
    al, rho = c[:, 2], c[:, 3]
    P = -(c[:, 6] + c[:, 7] + c[:, 9]) / 3
    full = al > 0.999
    load = full & (P > 2 * P0)
    ev = np.log(rho[load] / vc.RHO)
    ref = P0 / K + LAM * np.log(P[load] / P0)
    err = np.abs(ev - ref).max() / ref.max()
    print("[1] Cells on the compression line", flush=True)
    report("ln(rho/rho0) = P0/K + lambda* ln(P/P0)", err <= 0.02 and load.sum() >= 40,
           f"{load.sum()} cells, P {P[load].min():.0f} to {P[load].max():.0f} Pa, "
           f"largest error {err * 100:.2f} % of {ref.max():.4f}")
    print("[2] Mass and settlement", flush=True)
    dm = abs(h[-1, 2] / h[0, 2] - 1)
    settle = H - al.sum() * CELL
    umax = np.abs(c[full, 4:6]).max()
    report("mass", dm <= 1e-10, f"relative change {dm:.1e}")
    report("at rest, settled", umax <= 1e-4 and settle > 0.02,
           f"largest speed {umax:.1e} m/s; top down {settle * 1000:.1f} mm")
    c2, _ = case(work, False)
    top2 = c2[:, 2].sum() * CELL
    print(f"[3] Reported: without the cap the top settles {(H - top2) * 1000:.1f} mm (elastic only)",
          flush=True)
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
