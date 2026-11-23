"""WP K.4 done-when checks: the compaction cap with immersed (STL) walls
(decision 0122). The WP K.2 settling column (tests/compaction/settle_column.py)
held between two frictionless STL walls instead of mesh walls. About 1 minute.

Mesh: 3 columns of 20 mm, x from -15 to 45 mm. STL walls with faces at x = 0
and x = 20 mm: the left column is open over a quarter (a small cell, mixed
into its neighbour), the middle over three quarters, the right column is
covered. The material is the K.2 column (MC phi 30, E 10 MPa, cap
lambda* 0.01, P0 100 Pa, fullBlend), gravity ramped over 2 s with damping.

[1] Every loaded open cell lies on ln(rho/rho0) = P0/K + lambda* ln(P/P0)
    within 2 % (as K.2), and p_c = P there within 2 %.
[2] Mass constant to 1e-10; at rest; the top settles as the mesh-wall column
    (K.2 case) within 1 %.

    python3 settle_column_stl.py WORKDIR
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "stlwalls"))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import settle_column as k2  # noqa: E402
import stl_void_checks as sv  # noqa: E402
import void_checks as vc  # noqa: E402

fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


# wall faces (m); CEL_K4_SHIFT moves both (study: 0.005 puts them on mesh faces)
SHIFT = float(os.environ.get("CEL_K4_SHIFT", "0"))


def stl_case(work):
    vc.CELL = k2.CELL
    mat = (f"type mohrCoulomb;\nE {k2.E!r};\nnu {k2.NU!r};\nphi 30;\nc 0;\npsi 0;\n"
           + ("" if os.environ.get("CEL_K4_NOCAP") else f"capLambda {k2.LAM!r};\ncapP0 {k2.P0!r};\n")
           + "fullBlend yes;\n")
    check = dict(fill=[-1.0, 0.0, 1.0, k2.H], U0=[0.0, 0.0, 0.0], tEnd=3.0, maxCo=0.4)
    p = vc.make_case(f"{work}/stl", 3, k2.NFILL + k2.NVOID, mat, check)
    open(f"{p}/system/controlDict", "a").write("gravityRamp 2;\ndamping 100;\n")
    sv.move_mesh(p, -0.015, 0.0)
    top = (k2.NFILL + k2.NVOID + 1) * k2.CELL
    sv.add_bodies(p, {
        "leftWall": ((-0.2, -0.1, -0.01), (SHIFT, top, 0.03), 0.0, 0.0),
        "rightWall": ((0.02 + SHIFT, -0.1, -0.01), (0.2, top, 0.03), 0.0, 0.0)})
    vc.run(p)
    return np.loadtxt(f"{p}/cells.dat"), np.loadtxt(f"{p}/history.dat"), p


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    c, h, p = stl_case(work)
    # columns by cell centre x: left -5 mm (open 1/4), middle 15 mm (3/4)
    x = c[:, 0]
    mid = np.abs(x - 0.015) < 1e-6
    left = np.abs(x + 0.005) < 1e-6
    al, rho = c[:, 2], c[:, 3]
    P = -(c[:, 6] + c[:, 7] + c[:, 9]) / 3
    pc = c[:, -1]
    full = mid & (al > 0.999)
    load = full & (P > 2 * k2.P0)
    ev = np.log(rho[load] / vc.RHO)
    ref = k2.P0 / k2.K + k2.LAM * np.log(P[load] / k2.P0)
    err = np.abs(ev - ref).max() / ref.max()
    epc = np.abs(pc[load] / P[load] - 1).max()
    print("[1] Cells between STL walls on the compression line", flush=True)
    report("ln(rho/rho0) = P0/K + lambda* ln(P/P0)", err <= 0.02 and load.sum() >= 40,
           f"{load.sum()} cells, P {P[load].min():.0f} to {P[load].max():.0f} Pa, "
           f"largest error {err * 100:.2f} % of {ref.max():.4f}")
    report("p_c = P in loaded cells", epc <= 0.02, f"largest |p_c/P - 1| {epc * 100:.2f} %")
    lm = left & (al > 0.999) & (P > 2 * k2.P0)
    if lm.any():
        rm = {round(y, 6): r for y, r in zip(c[mid, 1], rho[mid])}
        d = max(abs(r - rm[round(y, 6)]) for y, r in zip(c[lm, 1], rho[lm])) / vc.RHO
        print(f"  reported: small (1/4) cells mixed with their neighbour, density "
              f"difference {d:.1e} of rho0", flush=True)
    print("[2] Mass, rest and settlement against the mesh-wall column", flush=True)
    dm = abs(h[-1, 2] / h[0, 2] - 1)
    umax = np.abs(c[full, 4:6]).max()
    settle = k2.H - al[mid].sum() * k2.CELL
    c2, _ = k2.case(work, True)
    settle2 = k2.H - c2[:, 2].sum() * k2.CELL
    rel = abs(settle / settle2 - 1)
    report("mass", dm <= 1e-10, f"relative change {dm:.1e}")
    report("at rest", umax <= 1e-4, f"largest speed {umax:.1e} m/s")
    report("settlement as the mesh-wall column", rel <= 0.01,
           f"STL {settle * 1000:.2f} mm, mesh walls {settle2 * 1000:.2f} mm, {rel * 100:.2f} %")
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
