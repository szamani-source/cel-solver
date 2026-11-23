"""Couto silo, MC, started with an empty gap along the wall (Sia,
2026-11-20): the material is placed inside the silo but the GAP_CELLS cells
next to the wall stay empty, and the material settles under gravity to fill
them; the at-rest pressures are read after settling, then the outlet opens
as in hopper_pressure_checks. Mesh: the wall-aligned mesh with the STL shell
inside it and rounded corners (buffered_case.aligned_buffered_case).

- Gap: the inner wall offset into the material by the width of the wall
  cells (measured on the mesh) times GAP_CELLS, normal to the wall, over the
  spout, both fillets, the cone and the cylinder.
- Mass: the fill top is raised so the Couto fill mass is placed; once the
  gap has filled, the surface sits near the silo top. Headroom for that.
- Timing: gravity ramp 0.1 s, no damping before T_DAMP so the material can
  slump into the gap freely, damping 100 1/s from T_DAMP, outlet opens at
  T_OPEN, run to T_END; discharge averaged over the last 0.75 s.

    python3 gap_case.py WORKDIR picture   # build, write the start, picture
    CEL_NPROCS=4 python3 gap_case.py WORKDIR run

Environment: CEL_H_TOP (m, 0.025): row height in the upper cylinder (the
columns there are about 12 mm already); the headroom stays 0.3 m.
CEL_NCOLS (44), CEL_WALL_RATIO (0.6, wall column / axis column) and
CEL_H_HOP (m, 0.012, row height at the top of the hopper block) set the
rest of the grading (decision 0115 addendum: the fine "ideal" mesh).
CEL_OUTLET_WALL_RATIO (default CEL_WALL_RATIO): the column grading at the
outlet only (decision 0117: 1 keeps the outlet cells at 3 mm or more).
"""
import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hopper_case as hc  # noqa: E402
import hopper_pressure_checks as hp  # noqa: E402
import buffered_case as bc  # noqa: E402

vc = hp.vc
NP = int(os.environ.get("CEL_NPROCS", "1"))
GAP_CELLS = 2
T_DAMP, T_OPEN, T_END = 0.4, 1.0, 2.0
H_TOP = float(os.environ.get("CEL_H_TOP", 0.025))
NCOLS = int(os.environ.get("CEL_NCOLS", 44))
WALL_RATIO = float(os.environ.get("CEL_WALL_RATIO", 0.6))
H_HOP = float(os.environ.get("CEL_H_HOP", 0.012))
OWR = os.environ.get("CEL_OUTLET_WALL_RATIO")      # decision 0117
HEAD_ROWS = int(round(12 * 0.025 / H_TOP))
CELLS = ["GA0", "GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7"]


def wall_cell_width(c):
    """Width across the wall of the open cell next to the wall in the
    cylinder: the column spacing at r = R from the aligned mesh's grading."""
    # columns axis -> M (R + TS/2), ncols 44, simpleGrading 0.6 (last/first)
    n, L, g = NCOLS, c["R"] + 0.5 * bc.TS, WALL_RATIO
    q = g ** (1.0 / (n - 1))
    a = L * (1 - q) / (1 - q ** n)
    edges = np.concatenate([[0.0], np.cumsum(a * q ** np.arange(n))])
    k = np.searchsorted(edges, c["R"]) - 1        # the column the wall cuts
    return float(edges[k + 1] - edges[k]), edges


def gap_profile(c, gap, y_top):
    """(r, y) of the fill's right edge: the inner wall offset by gap into
    the material."""
    W = bc.wall_path(c, y_top)
    G = bc.offset_path(W, -gap)
    G = G[G[:, 0] > 0]
    return G[np.argsort(G[:, 1])]


def build(work, check):
    c = dict(hc.COUTO)
    b, ht = hc.hopper_geometry(c)
    H = ht + c["Hc"]
    w, _ = wall_cell_width(c)
    gap = GAP_CELLS * w
    Hm = H + HEAD_ROWS * H_TOP
    G = gap_profile(c, gap, Hm)
    # Couto fill mass in the full silo (incl. the spout) and in the gap fill
    prof = np.vstack([[c["Rout"], bc.Y0], bc.inner_profile(c)])
    ys = np.linspace(bc.Y0, H, 4001)
    full = getattr(np, "trapezoid", getattr(np, "trapz", None))(math.pi * np.interp(ys, prof[:, 1], prof[:, 0]) ** 2, ys)
    yy = np.linspace(bc.Y0, Hm, 8001)
    area = math.pi * np.clip(np.interp(yy, G[:, 1], G[:, 0]), 0, None) ** 2
    cum = np.concatenate([[0], np.cumsum(0.5 * (area[1:] + area[:-1]) * np.diff(yy))])
    y_fill = float(np.interp(full, cum, yy))
    ct = math.cos(math.radians(hc.THETA))
    check = dict(check, fill=[-1.0, bc.Y0, 1.0, y_fill],
                 fillRight=[[float(y), float(r * ct)] for r, y in G[::2]]
                 + [[Hm + 1.0, float(G[-1, 0] * ct)]])
    p, c = bc.aligned_buffered_case(work, "couto_mc_gap", "mc", check, c, nProcs=NP,
                                    head_rows=HEAD_ROWS, h_top=H_TOP,
                                    ncols=NCOLS, wall_ratio=WALL_RATIO, h_hop=H_HOP,
                                    **({} if OWR is None else dict(outlet_wall_ratio=float(OWR))))
    dp = f"{p}/system/decomposeParDict"
    txt = open(dp).read().replace(f"n ({NP} 1 1)", f"n (1 {NP} 1)")
    open(dp, "w").write(txt)
    open(f"{p}/system/controlDict", "a").write(
        f"gravityRamp 0.1;\ndamping 100;\ndampingStart {T_DAMP!r};\n")
    info = dict(gap=gap, w=w, y_fill=y_fill, H=H, mass=full * c["rho"])
    return p, c, info


def picture(work):
    """Run one step on one core and draw the start: alpha per cell."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    global NP
    NP = 1
    p, c, info = build(work, dict(tEnd=1e-6, maxCo=0.4, writeTimes=[]))
    vc.run(p)
    a = np.loadtxt(f"{p}/cells.dat")
    sys.path.insert(0, os.path.join(HERE, "..", "angled"))
    import angled_checks as ac
    m = f"{p}/constant/polyMesh"
    pts = ac.foam_list(f"{m}/points", "point")
    faces = ac.foam_list(f"{m}/faces", "face")
    t = math.radians(hc.THETA)
    polys = [np.column_stack([pts[f][:, 0] / math.cos(t), pts[f][:, 1]])
             for f in faces if np.all(pts[f][:, 2] <= 1e-12) and len(f) == 4]
    cen = np.array([q.mean(0) for q in polys])
    # alpha x open of the nearest cell centre (cells.dat order differs)
    r = a[:, 0] / math.cos(t)
    P = np.column_stack([r, a[:, 1]])
    k = np.array([np.argmin(((P - q) ** 2).sum(1)) for q in cen])
    val = a[k, 2] * a[k, 10]
    poly = bc.shell_polygon(c, y_top=info["H"] + HEAD_ROWS * H_TOP)
    fig, axs = plt.subplots(1, 4, figsize=(19, 9), gridspec_kw=dict(width_ratios=[1.1, 1, 1, 1]))
    views = [((-0.02, 0.6), (-0.05, 2.85), "Whole silo"),
             ((0.10, 0.32), (-0.03, 0.16), "Outlet edge"),
             ((0.36, 0.60), (0.36, 0.62), "Transition"),
             ((0.36, 0.60), (1.40, 1.66), f"Upper cylinder ({H_TOP * 1000:.0f} mm rows)")]
    for ax, (xl, yl, title) in zip(axs, views):
        pc = PolyCollection(polys, array=val, cmap="YlOrBr", clim=(0, 1),
                            edgecolors="#9aa4b1", linewidths=0.25)
        ax.add_collection(pc)
        pp = np.vstack([poly, poly[:1]])
        ax.plot(pp[:, 0], pp[:, 1], color="#c0392b", lw=1.2, label="STL shell (20 mm)")
        ax.axhline(info["H"], color="#1f6fb4", lw=0.8, ls="--", label="silo top (Couto fill level)")
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("r (m)")
    axs[0].set_ylabel("y (m)")
    axs[0].legend(loc="upper left", fontsize=8)
    fig.colorbar(pc, ax=axs, shrink=0.6, label="material fraction at the start")
    fig.suptitle(f"Couto silo start, MC: material placed with a {info['gap'] * 1000:.0f} mm "
                 f"empty gap along the wall ({GAP_CELLS} wall cells of {info['w'] * 1000:.1f} mm); "
                 f"fill top {info['y_fill']:.3f} m so the Couto mass ({info['mass']:.0f} kg) is in")
    png = os.environ.get("CEL_PNG", "/mnt/project-files/hopper/couto-gap-start.png")
    fig.savefig(png, dpi=130)
    print(png, info)


def run(work):
    t0 = time.perf_counter()
    writes = [T_OPEN - 0.2, T_OPEN] + [round(T_OPEN + 0.05 * i, 2) for i in range(1, 21)]
    os.environ.setdefault("CEL_HOPPER_RULES", hp.MU_I + hp.RULES + "fullBlend yes;\n")
    p, c, info = build(work, dict(tEnd=T_END, maxCo=0.4, writeTimes=writes, tOpen=T_OPEN))
    print(info, flush=True)
    vc.run(p, NP)
    hp.CELL = 0.012
    h = vc.history(p)
    t = h[:, vc.T]
    f = math.pi / (math.sin(math.radians(hc.THETA)) * math.cos(math.radians(hc.THETA)))
    m = h[:, vc.MASS] * f
    i = (t >= T_OPEN + 0.25) & (t <= T_END)
    print(f"run {time.perf_counter() - t0:.0f} s; mass at opening {m[np.argmin(abs(t - T_OPEN))]:.0f} kg; "
          f"discharge {-np.polyfit(t[i], m[i], 1)[0]:.1f} kg/s over 0.25-1.0 s after opening; "
          f"largest speed before opening at the end of settling "
          f"{h[(t > T_OPEN - 0.05) & (t < T_OPEN), vc.UMAX].max():.3f} m/s")
    rows = hp.table()
    for w in (T_OPEN - 0.2, T_OPEN):
        r = hp.pressures(p, c, f"cells_{w:g}")
        print(f"  rest t {w:g}: " + " ".join(f"{n} {r[n]:.2f}" for n in CELLS))
    print("  experiment at rest: " + " ".join(
        "{} {:.2f}-{:.2f}".format(n, *hp.exp_range(rows, "static", n)) for n in CELLS))
    dis = [hp.pressures(p, c, f"cells_{w:g}") for w in writes if T_OPEN + 0.25 <= w <= T_END]
    for n in CELLS:
        v = np.array([d[n] for d in dis])
        lo, hi = hp.exp_range(rows, "discharge", n)
        print(f"  discharge {n} {v.mean():.2f} (min {v.min():.2f}, max {v.max():.2f}), "
              f"experiment {lo:.2f}-{hi:.2f}")
    for w in (0.0001, T_OPEN, T_END):
        try:
            a = hp.read_cells(p, f"cells_{w:g}")
        except Exception:
            continue
        out, rr, yy = bc.leak_report(c, a)
        print(f"  leak t {w:g}: {out.sum()} cells, {hp.leak_volume(p, a, out) * c['rho']:.4f} kg")
    a = hp.read_cells(p, f"cells_{T_OPEN:g}")
    top = a[(a[:, 2] > 0.5) & (a[:, 10] > 0) & (a[:, 0] < 0.3)]
    print(f"  surface near the axis at opening {top[:, 1].max():.3f} m (silo top {info['H']:.3f} m)")


if __name__ == "__main__":
    work = os.path.abspath(sys.argv[1])
    (picture if sys.argv[2] == "picture" else run)(work)
