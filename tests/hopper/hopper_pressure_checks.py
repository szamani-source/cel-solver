"""WP H.3 done-when checks: Couto silo wall pressure at rest and in discharge,
Mohr-Coulomb or Drucker-Prager (docs/spec/conical-hopper.md case A; decisions
0100 to 0102). Manual: about 12 minutes per model at 25 mm cells.

The silo of tests/hopper/hopper_case.py with the T0 rule set (mu(I) with a
small rate term, flow impedance, shear and normal caps, layer spreading;
decision 0102): settles for 0.5 s (gravity ramp 0.1 s, damping 100 1/s), the
outlet opens at 0.5 s, run to 1.5 s. Pressures at rest are read at 0.5 s;
discharge pressures are averaged over the writes from 0.25 to 1.0 s after
opening (every 0.05 s), because the cell values swing with the flow.

Bands (spec page): at rest GA1 to GA3 and GA0 within the experiment's range
widened by 15 %; GA4 to GA7 within that or 15 % of Janssen (K_A 0.338, mu 0.2).
In discharge GA0, GA1, GA2 within the experiment's discharge range widened by
15 %. Also reported: every cell, the discharge rate, the largest speed.

WP H.7 (decisions 0103, 0104): by default the silo is on the wall-aligned
mesh (hopper_case.aligned_case: hopper wall a mesh patch, outlet rows 3x
finer, Sia 2026-11-19) with `fullBlend yes` added to the rule set; `cut`
runs the cut-cell (STL) silo of WP H.3 with the H.3 rules.

    python3 hopper_pressure_checks.py WORKDIR [mc|dp] [aligned|cut]
"""

import csv
import glob
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
NP = int(os.environ.get("CEL_NPROCS", "1"))
T_OPEN, T_END = 0.5, 1.5
WRITES = [round(T_OPEN + 0.05 * i, 2) for i in range(0, 21)]
AVG = (0.75, 1.5)
MU_I = "mu2 0.8412;\nI0 0.279;\nd 0.004;\nrhoS 1380;\n"   # decision 0102
RULES = "shearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n"
GA = {"GA7": 2.075, "GA6": 1.825, "GA5": 1.575, "GA4": 1.325, "GA3": 1.075,
      "GA2": 0.825, "GA1": 0.575}
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))


def table():
    f = os.path.join(ROOT, "docs", "spec", "hopper-data", "couto2013_table1.csv")
    with open(f) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def exp_range(rows, phase, cell):
    v = [float(r[cell]) for r in rows if r["phase"] == phase]
    return min(v), max(v)


def read_cells(path, name):
    files = sorted(glob.glob(f"{path}/{name}.dat.proc*")) or [f"{path}/{name}.dat"]
    return np.vstack([np.loadtxt(f) for f in files])


def pressures(path, c, name):
    np.savetxt(f"{path}/cells.dat", read_cells(path, name),
               header="x y alpha rho Ux Uy sxx syy sxy szz open")
    cyl, hop = hc.wall_pressures(path, c, CELL)
    out = {n: cyl[int(np.argmin(abs(cyl[:, 0] - z))), 1] for n, z in GA.items()}
    out["GA0"] = hop[np.abs(hop[:, 0] - 0.41) < 1.5 * CELL][:, 1].mean()
    return out


def leak_volume(p, a, out):
    """Full-silo material volume of the flagged cells: alpha x open volume,
    the cell volume from the wedge (r dr dy over the wedge angle)."""
    if not out.any():
        return 0.0
    # cell sizes from neighbouring centres are not stored: use the open
    # fraction and a local size estimate from the nearest distinct centres
    x, y = a[:, 0], a[:, 1]
    vol = 0.0
    for i in np.where(out)[0]:
        dx = np.min(np.abs(x[np.abs(y - y[i]) < 1e-9] - x[i])[np.abs(x[np.abs(y - y[i]) < 1e-9] - x[i]) > 1e-12], initial=0.007)
        dy = np.min(np.abs(y[np.abs(x - x[i]) < 1e-9] - y[i])[np.abs(y[np.abs(x - x[i]) < 1e-9] - y[i]) > 1e-12], initial=0.007)
        vol += a[i, 2] * a[i, 10] * 2 * math.pi * x[i] * dx * dy
    return vol


def main():
    global CELL
    work = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else "mc"
    t0 = time.perf_counter()
    mesh = sys.argv[3] if len(sys.argv) > 3 else "aligned"
    blend = "fullBlend yes;\n" if mesh != "cut" else ""
    os.environ.setdefault("CEL_HOPPER_RULES", MU_I + RULES + blend)
    check = dict(tEnd=T_END, maxCo=0.4, writeTimes=WRITES, tOpen=T_OPEN)
    if mesh == "aligned":
        # grading (Sia): CEL_OUTLET_RATIO (outlet row over top row, 1/3),
        # CEL_WALL_RATIO (wall column over axis column, 1), CEL_NCOLS
        p, c = hc.aligned_case(
            work, f"couto_{model}_{mesh}", model, CELL, check, nProcs=NP,
            outlet_ratio=float(os.environ.get("CEL_OUTLET_RATIO", 1 / 3)),
            wall_ratio=float(os.environ.get("CEL_WALL_RATIO", 1.0)),
            ncols=int(os.environ["CEL_NCOLS"]) if "CEL_NCOLS" in os.environ else None)
    elif mesh == "aligned_buffered":
        # STL shell inside a wall-aligned mesh (Sia, 2026-11-20; WP G.5)
        import buffered_case as bc
        CELL = 0.012
        p, c = bc.aligned_buffered_case(work, f"couto_{model}_{mesh}", model, check,
                                        nProcs=NP)
    elif mesh == "buffered":
        # STL shell inside a larger mesh, rounded corners (Sia, 2026-11-20)
        import buffered_case as bc
        CELL = 0.007
        p, c = bc.buffered_case(work, f"couto_{model}_{mesh}", model, check, nProcs=NP)
    else:
        p, c = hc.couto_case(work, f"couto_{model}_{mesh}", model, CELL, check,
                             nProcs=NP)
    dp = f"{p}/system/decomposeParDict"
    txt = open(dp).read().replace(f"n ({NP} 1 1)", f"n (1 {NP} 1)")
    open(dp, "w").write(txt)
    open(f"{p}/system/controlDict", "a").write("gravityRamp 0.1;\ndamping 100;\n")
    wall = vc.run(p, NP)
    rows = table()
    fails = []

    def check(ok, text):
        print(f"  {text}: {'PASS' if ok else 'FAIL'}")
        if not ok:
            fails.append(text)

    h = np.loadtxt(f"{p}/history.dat")
    t, m = h[:, 1], h[:, 2] * math.pi / (math.sin(math.radians(hc.THETA))
                                         * math.cos(math.radians(hc.THETA)))
    i, j = np.searchsorted(t, AVG[0]), np.searchsorted(t, AVG[1]) - 1
    rate = -(m[j] - m[i]) / (t[j] - t[i])
    print(f"[1] Couto silo, {model.upper()}, {mesh} mesh, {CELL * 1e3:.1f} mm cells ({wall:.0f} s): "
          f"discharge {rate:.1f} kg/s over {AVG[0] - T_OPEN:.2f} to "
          f"{AVG[1] - T_OPEN:.2f} s after opening; largest speed after opening "
          f"{h[t > T_OPEN, 9].max():.1f} m/s")
    if mesh in ("buffered", "aligned_buffered"):
        import buffered_case as bc
        print("[leak] material outside the STL shell (full-silo kg; where)")
        for w in WRITES[::4]:
            a = read_cells(p, f"cells_{w:g}")
            out, r, y = bc.leak_report(c, a)
            vol = leak_volume(p, a, out)
            print(f"  t {w:g}: {out.sum()} cells, {vol * c['rho']:.4f} kg"
                  + (f", r {r[out].min():.3f}-{r[out].max():.3f} m, y {y[out].min():.3f}-"
                     f"{y[out].max():.3f} m" if out.any() else ""))
    rest = pressures(p, c, f"cells_{T_OPEN:g}")
    H = hc.hopper_geometry(c)[1] + c["Hc"]
    K, mu = 0.338, c["mu"]
    print("[2] At rest")
    for n in ["GA0", "GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7"]:
        lo, hi = exp_range(rows, "static", n)
        v = rest[n]
        ok = lo / 1.15 <= v <= hi * 1.15
        txt = f"{n} {v:.2f} kPa, experiment {lo:.2f} to {hi:.2f}"
        if n in ("GA4", "GA5", "GA6", "GA7"):
            pj = c["rho"] * 9.81 / (4 * mu) * (1 - math.exp(-4 * mu * K * (H - GA[n]))) / 1e3
            ok = ok or abs(v / pj - 1) <= 0.15
            txt += f", Janssen K_A {pj:.2f}"
        check(ok, txt)
    print(f"[3] Discharge, mean of {len([w for w in WRITES if AVG[0] <= w <= AVG[1]])} "
          "writes")
    dis = [pressures(p, c, f"cells_{w:g}") for w in WRITES if AVG[0] <= w <= AVG[1]]
    for n in ["GA0", "GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7"]:
        v = np.array([d[n] for d in dis])
        lo, hi = exp_range(rows, "discharge", n)
        txt = (f"{n} {v.mean():.2f} kPa (min {v.min():.2f}, max {v.max():.2f}), "
               f"experiment {lo:.2f} to {hi:.2f}")
        if n in ("GA0", "GA1", "GA2"):
            check(lo / 1.15 <= v.mean() <= hi * 1.15, txt)
        else:
            print(f"  report: {txt}")
    print(f"{len(fails)} failures ({time.perf_counter() - t0:.0f} s)")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
