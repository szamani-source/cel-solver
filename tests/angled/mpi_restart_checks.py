"""WP G.6 done-when checks: MPI and restart on angled (leaning) cells.
About 1 minute.

[1] Parallel: the G.2 granular column collapse (a = 0.5, 1 MPa) on cells
    leaning right, to t = 0.1 s, serial against 4 processors (2 x 2 split,
    so processor boundaries cross the leaning columns and the moving
    surface), without the shear cap. Scaled differences <= 1e-10, as for
    the Cartesian collapse (tests/test_collapse.sh): velocity by 1 m/s,
    stress by rho g H_i. With the shear cap on, the cap's switch flips on
    round-off within a few steps on leaning cells and the runs part at
    about 1e-8 (decision 0110); that run is reported, not checked.
[2] Restart: the siloDischarge tutorial tank (T0 rules, E 0.1 MPa) with its
    three blocks leaning (the outlet block widens upward, the side blocks
    lean outward), outlet open at 0.3 s, continuous to 0.5 s against a
    restart from 0.4 s while it discharges: fields, state and mass counters
    identical bit for bit (restart_checks.check). Then the same on 2
    processors.

    python3 mpi_restart_checks.py WORKDIR [parallel] [restart] [restart2]
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("void", "lagrange", "collapse", "restart"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import void_checks as vc  # noqa: E402
import collapse_checks as cc  # noqa: E402
import restart_checks as rc  # noqa: E402
from free_surface_checks import split_mesh  # noqa: E402

fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def parallel_pair(work, tag, cap):
    L_i, H_i, nx, ny = cc.geometry(0.5)
    W, H, d = nx * vc.CELL, ny * vc.CELL, vc.CELL
    ps = cc.setup(work, 0.5, 1.0, 0.1, tag + "Ser")
    pp = cc.setup(work, 0.5, 1.0, 0.1, tag + "Par", nProcs=4)
    for p in (ps, pp):
        split_mesh(p, W, H, nx, ny, 0.15, 0.25, d)
        if not cap:
            mp = f"{p}/constant/materialProperties"
            txt = open(mp).read()
            open(mp, "w").write(txt.replace("shearCap yes;", "shearCap no;"))
    dp = f"{pp}/system/decomposeParDict"
    txt = open(dp).read()
    open(dp, "w").write(txt.replace("n (4 1 1)", "n (2 2 1)"))
    vc.run(ps)
    vc.run(pp, 4)
    a, b = vc.cells(ps), vc.cells(pp, 4)
    if a.shape != b.shape:
        return None, 0.0, ""
    sg = cc.RHO * 9.81 * 0.1
    scale = {"alpha": (2, 1.0), "rho": (3, cc.RHO), "Ux": (4, 1.0), "Uy": (5, 1.0),
             "s_xx": (6, sg), "s_yy": (7, sg), "s_xy": (8, sg)}
    dd = {k: np.abs(a[:, i] - b[:, i]).max() / sc for k, (i, sc) in scale.items()}
    return max(dd.values()), np.abs(a[:, 4]).max(), \
        ", ".join(f"{k} {v:.0e}" for k, v in dd.items())


def check_parallel(work):
    print("[1] Column collapse on leaning cells: 4 processors against serial", flush=True)
    worst, moved, det = parallel_pair(work, "lean", cap=False)
    report("4 processors, no shear cap", worst is not None and worst <= 1e-10 and moved > 0.01,
           f"largest |Ux| {moved:.3f} m/s; largest scaled difference {det} (<= 1e-10)")
    worst, moved, det = parallel_pair(work, "leanCap", cap=True)
    print(f"  report: with the shear cap, largest scaled difference {det}", flush=True)


def lean_silo(path, t_restart, t_end):
    """siloDischarge with leaning blocks: top ends of the inner block lines
    moved out by 50 mm."""
    rc.make_silo(path, t_restart, t_end)
    bm = f"{path}/system/blockMeshDict"
    s = open(bm).read()
    for z in ("0", "0.01"):
        s = s.replace(f"(0.17 0.62 {z})", f"(0.12 0.62 {z})")
        s = s.replace(f"(0.23 0.62 {z})", f"(0.28 0.62 {z})")
    s = s.replace("(17 62 1)", "(15 62 1)").replace("(6 62 1)", "(10 62 1)")
    open(bm, "w").write(s)
    return path


def check_restart(work, nProcs):
    tag = "restart" if nProcs == 1 else f"restart{nProcs}"
    ok = rc.check(work, f"lean_{tag}", lean_silo, 0.4, 0.5, nProcs=nProcs)
    if not ok:
        fails.append(tag)


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    only = sys.argv[2:] or ["parallel", "restart", "restart2"]
    t0 = time.perf_counter()
    if "parallel" in only:
        check_parallel(work)
    if "restart" in only:
        check_restart(work, 1)
    if "restart2" in only:
        check_restart(work, 2)
    print(f"{len(fails)} failures ({time.perf_counter() - t0:.0f} s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
