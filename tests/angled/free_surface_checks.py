"""WP G.2 done-when checks: free surface on angled (leaning) cells
(decision 0104). About 2 minutes.

Meshes: a box of vertical walls meshed as two blocks split by an inclined
line, so every column of cells leans and the cells differ in shape; the
same box as one Cartesian block is the reference.

[1] Lake at rest: T0 Mohr-Coulomb (void_checks.MC) between smooth walls in a
    0.5 x 1.0 m box filled to 0.6 m, 25 mm rows, gravity ramp 0.1 s, damping
    100 1/s, 0.5 s: mass to 1e-10, surface flat (every column's top cell at
    the same height), s_yy hydrostatic within 5 % of rho g H at most and
    0.5 % on average below 0.5 m.
[2] Granular column collapse a = 0.5, E = 1 MPa (tests/test_collapse.sh
    set-up: shear cap, layer spreading), on cells leaning both ways: mass to
    1e-10, stops (kinetic energy below 0.1 of its peak at 1 s), front within
    3 cells and largest height within 1 cell of the Cartesian run, surface
    profile within 1.5 cells (RMS over 2 cm bands).

    python3 free_surface_checks.py WORKDIR
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("void", "lagrange", "collapse"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import void_checks as vc  # noqa: E402
import collapse_checks as cc  # noqa: E402

fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def split_mesh(p, W, H, nx, ny, X0, X1, depth):
    """Box W x H as two blocks split by the line (X0, 0) to (X1, H)."""
    v = [(0, 0, 0), (X0, 0, 0), (X1, H, 0), (0, H, 0), (W, 0, 0), (W, H, 0)]
    v += [(x, y, depth) for x, y, _ in v]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    n1 = int(round(nx * (X0 + X1) / 2 / W))
    cs.write(f"{p}/system/blockMeshDict", cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n"
             f"    hex (0 1 2 3 6 7 8 9) ({n1} {ny} 1) simpleGrading (1 1 1)\n"
             f"    hex (1 4 5 2 7 10 11 8) ({nx - n1} {ny} 1) simpleGrading (1 1 1)\n);\n"
             "boundary\n(\n"
             "    left { type patch; faces ((0 6 9 3)); }\n"
             "    right { type patch; faces ((4 5 11 10)); }\n"
             "    bottom { type patch; faces ((0 1 7 6) (1 4 10 7)); }\n"
             "    top { type patch; faces ((3 9 8 2) (2 8 11 5)); }\n"
             "    frontAndBack { type empty; faces ((0 3 2 1) (1 2 5 4) (6 7 8 9) "
             "(7 10 11 8)); }\n);\n")


def check_lake(work):
    print("[1] Lake at rest on leaning cells", flush=True)
    W, H, nx, ny, fill = 0.5, 1.0, 20, 40, 0.6
    p = vc.make_case(f"{work}/lake", nx, ny, vc.MC,
                     dict(fill=[-1, -1, 1, fill], U0=[0, 0, 0], tEnd=0.5, maxCo=0.4))
    split_mesh(p, W, H, nx, ny, 0.1, 0.4, W / nx)
    open(f"{p}/system/controlDict", "a").write("gravityRamp 0.1;\ndamping 100;\n")
    vc.run(p)
    h = vc.history(p)
    c = np.loadtxt(f"{p}/cells.dat")
    x, y, a, syy = c[:, 0], c[:, 1], c[:, 2], c[:, 7]
    dm = abs(h[-1, vc.MASS] / h[0, vc.MASS] - 1)
    report("mass", dm <= 1e-10, f"relative change {dm:.1e}")
    tops = y[a > 1e-6]
    top = tops.max()
    flat = np.all(np.abs(y[(a > 0.5) & (y > top - 0.5 * H / ny)] - top) < 1e-3) and \
        ((a > 1e-6) & (y > fill + H / ny)).sum() == 0
    report("surface flat", bool(flat), f"top row at {top:.4f} m, nothing above {fill + H / ny:.3f} m")
    band = (a > 0.99) & (y > 0.1) & (y < 0.5)
    e = np.abs(syy[band] + vc.RHO * vc.G_ACC * (fill - y[band])) / (vc.RHO * vc.G_ACC * fill)
    report("s_yy hydrostatic", e.max() <= 0.05 and e.mean() <= 0.005,
           f"largest {e.max() * 100:.2f} %, mean {e.mean() * 100:.2f} % of rho g H")


def profile(cells):
    c = np.loadtxt(cells)
    x, y, a = c[:, 0], c[:, 1], c[:, 2]
    m = a >= 0.5
    bands = np.arange(0.0, 0.36, 0.02)
    prof = np.array([y[m & (x >= b) & (x < b + 0.02)].max(initial=0.0) for b in bands])
    return x[m].max(), y[m].max(), prof


def check_collapse(work):
    print("[2] Column collapse on leaning cells against Cartesian", flush=True)
    L_i, H_i, nx, ny = cc.geometry(0.5)
    W, H, d = nx * vc.CELL, ny * vc.CELL, vc.CELL
    res = {}
    for name, X0, X1 in (("cartesian", W / 2, W / 2), ("lean right", 0.15, 0.25),
                         ("lean left", 0.25, 0.15)):
        p = cc.setup(work, 0.5, 1.0, 1.0, name.replace(" ", "_"))
        split_mesh(p, W, H, nx, ny, X0, X1, d)
        vc.run(p)
        h = vc.history(p)
        res[name] = profile(f"{p}/cells.dat") + (h,)
    f0, t0, p0, _ = res["cartesian"]
    for name in ("lean right", "lean left"):
        f, t, pr, h = res[name]
        dm = abs(h[-1, vc.MASS] / h[0, vc.MASS] - 1)
        ke = h[:, cc.KE]
        rms = np.sqrt(np.mean((pr - p0) ** 2)) / d
        ok = (dm <= 1e-10 and ke[-1] <= 0.1 * ke.max() and abs(f - f0) <= 3 * d
              and abs(t - t0) <= d and rms <= 1.5)
        report(name, ok, f"mass {dm:.1e}, KE end/peak {ke[-1] / ke.max():.1e}, front "
               f"{f:.3f} m (Cartesian {f0:.3f}), height {t:.3f} m ({t0:.3f}), "
               f"profile RMS {rms:.2f} cells")


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    check_lake(work)
    check_collapse(work)
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
