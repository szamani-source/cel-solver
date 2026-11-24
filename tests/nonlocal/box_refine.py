"""WP N.4 (decision 0123, Sia 2026-11-23): mesh refinement of a small 2D
box discharge with and without the nonlocal law, in place of the Couto
hopper (too slow for a refinement study).

The T0 tank of tests/outlet/outlet_checks.py (W 0.4 m, fill H 0.6 m,
outlet D 0.06 m centred in the floor, walls and floor 20 deg, Mohr-Coulomb
rho 1800, phi 35, c 0, psi 0, E 10 MPa as tutorials/siloDischarge, outlet
opens at 0.3 s), with the mu(I) set of decision 0037 and the grain size
d = 5 mm (D = 12 d, so Beverloo applies; the nonlocal band of a few d is
resolved on 5 and 2.5 mm cells, decision 0111 N.2: cells of 1 d or finer).

Laws: "local" = mu(I) with flowImpedance (the law the solver had before
the nonlocal one); "NL" = the same plus nonlocalA 0.48 (Henann and Kamrin
2013). Meshes: 10, 5 and 2.5 mm square cells (outlet 6, 12, 24 cells).

Measured per run: the steady rate (kg/s per metre, tank mass 97 to 88 %
of the mass at opening, as WP 3.3) and Beverloo C; at 0.3 and 0.5 s after
opening, the vertical speed across the row at y = 0.15 m (2.5 D above the
outlet): centreline speed and flowing-zone width (where the speed is at
least 10 % of the row's largest).

    python3 box_refine.py WORKDIR [cell_mm ...]       # default 10 5 2.5
    python3 box_refine.py WORKDIR --picture PNG        # draw the meshes only

CEL_LAWS (default "local NL"), CEL_JOBS (default 4), CEL_CHECKPOINT (s,
default 0.05). Put WORKDIR on /mnt/project-files for long runs: run the
same command again after a lost container and each case resumes.
"""

import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "discharge"))
import discharge_checks as dc  # noqa: E402

oc = dc.oc
E_MPA = 10.0
D_GRAIN = 0.005
MUI = (f"mu2 0.961813;\nI0 0.279;\nd {D_GRAIN!r};\nrhoS 2650.0;\n"
       "flowImpedance yes;\nnormalCap yes;\n")
LAWS = {"local": MUI, "NL": MUI + "nonlocalA 0.48;\n"}
SNAP = [0.3, 0.5]                  # s after opening
Y_ROW = 0.15
CKPT = float(os.environ.get("CEL_CHECKPOINT", 0.05))


def row_profile(path, t, h):
    c = np.loadtxt(f"{path}/cells_{oc.T_OPEN + t:g}.dat")
    j = np.floor(c[:, 1] / h + 1e-9).astype(int)
    s = (j == int(Y_ROW / h)) & (c[:, 2] > 0.5)
    x, v = c[s, 0], -c[s, 5]
    o = np.argsort(x)
    x, v = x[o], v[o]
    vc_ = np.interp(oc.W / 2, x, v)
    flow = x[v >= 0.1 * v.max()]
    return vc_, (flow.max() - flow.min() + h) if flow.size else 0.0


def run_case(sub):
    """dc.case with checkpoints every CKPT s of simulated time (binary,
    exact restart, decision 0117): run again after a lost container and it
    resumes from the latest complete checkpoint."""
    oc.E = E_MPA * 1e6
    oc.T_DIS = dc.T_DIS
    path = f"{sub}/E{E_MPA:g}"
    if oc.vc.latest_time(path) > 0:
        sec = oc.vc.resume(path)
    else:
        oc.CHECK_EXTRA = {"checkpointEvery": CKPT}
        path = oc.setup(sub, f"E{E_MPA:g}")
        open(f"{path}/system/controlDict", "a").write("writeFormat binary;\n")
        sec = oc.vc.run(path)
    # keep the latest two checkpoints
    ts = sorted((t for t in os.listdir(path) if t.replace(".", "").isdigit() and float(t) > 0),
                key=float)
    import shutil
    for t in ts[:-2]:
        shutil.rmtree(f"{path}/{t}", ignore_errors=True)
    h = np.loadtxt(f"{path}/history.dat")
    t, m, out = h[:, oc.T], h[:, oc.MASS] / oc.CELL, h[:, oc.OUT] / oc.CELL
    i_open = int(np.where(t <= oc.T_OPEN + 1e-12)[0][-1])
    m_open = m[i_open]
    rate, win = oc.steady_rate(t[i_open:] - oc.T_OPEN, m[i_open:], m_open)
    C = rate / (oc.RHO * math.sqrt(oc.G_ACC) * oc.D ** 1.5)
    return dict(rate=rate, C=C, win=win, m_open=m_open,
                cons=np.abs((m + out) / m[0] - 1.0).max(), sec=sec)


def one(job):
    work, h, law = job
    oc.CELL = oc.vc.CELL = h
    oc.WRITE_TIMES = [oc.T_OPEN + t for t in SNAP]
    os.environ["CEL_EXTRA_MATERIAL"] = LAWS[law]
    sub = f"{work}/{law}_{h * 1000:g}mm"
    os.makedirs(sub, exist_ok=True)
    r = run_case(sub)
    path = f"{sub}/E{E_MPA:g}"
    for t in SNAP:
        r[f"vc{t}"], r[f"w{t}"] = row_profile(path, t, h)
    r.update(law=law, h=h)
    return r


def picture(png, sizes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, len(sizes), figsize=(4.2 * len(sizes), 9),
                            gridspec_kw=dict(height_ratios=[2.2, 1]))
    H = oc.H + oc.HEADROOM * max(sizes)
    for k, h in enumerate(sizes):
        nx = int(round(oc.W / h))
        n_out = int(round(oc.D / h))
        x0 = (nx - n_out) // 2 * h
        for ax, (xl, yl) in zip(axs[:, k], [((0, oc.W), (0, H)),
                                            ((0.12, 0.28), (0, 0.10))]):
            xs = np.arange(0, oc.W + h / 2, h)
            ys = np.arange(0, H + h / 2, h)
            ax.add_patch(plt.Rectangle((0, 0), oc.W, oc.H, color="#e8d9b5", zorder=0))
            ax.vlines(xs, 0, ys[-1], lw=0.3 if h < 0.004 else 0.5, color="#555555")
            ax.hlines(ys, 0, oc.W, lw=0.3 if h < 0.004 else 0.5, color="#555555")
            ax.plot([0, x0], [0, 0], color="#000000", lw=3)
            ax.plot([x0 + oc.D, oc.W], [0, 0], color="#000000", lw=3)
            ax.plot([x0, x0 + oc.D], [0, 0], color="#d62728", lw=3, label="outlet")
            ax.plot([0, 0], [0, H], color="#000000", lw=3)
            ax.plot([oc.W, oc.W], [0, H], color="#000000", lw=3)
            ax.set_xlim(xl[0] - 0.005, xl[1] + 0.005)
            ax.set_ylim(yl[0] - 0.005, yl[1] + 0.005)
            ax.set_aspect("equal")
            ax.set_xlabel("x (m)")
        axs[0, k].set_title(f"{h * 1000:g} mm cells: {nx} x {int(math.ceil(oc.H / h)) + oc.HEADROOM}"
                            f" = {nx * (int(math.ceil(oc.H / h)) + oc.HEADROOM)} cells")
        axs[1, k].set_title(f"outlet close up: {n_out} cells across D = {oc.D * 1000:g} mm")
    axs[0, 0].set_ylabel("y (m)")
    axs[1, 0].legend(loc="upper left")
    fig.suptitle(f"Box discharge for the nonlocal refinement check: W {oc.W} m, fill {oc.H} m "
                 f"(sand colour), outlet {oc.D * 1000:g} mm, grain d {D_GRAIN * 1000:g} mm")
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    print(png)


def main():
    args = sys.argv[1:]
    work = args.pop(0)
    if "--picture" in args:
        png = args[args.index("--picture") + 1]
        picture(png, [0.01, 0.005, 0.0025])
        return 0
    sizes = [float(a) / 1000 for a in args] or [0.01, 0.005, 0.0025]
    laws = os.environ.get("CEL_LAWS", "local NL").split()
    jobs = [(work, h, law) for h in sorted(sizes) for law in laws]
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=int(os.environ.get("CEL_JOBS", 4))) as ex:
        res = list(ex.map(one, jobs))
    res.sort(key=lambda r: (r["law"], -r["h"]))
    print(f"Box discharge, E {E_MPA:g} MPa, d {D_GRAIN * 1000:g} mm, outlet {oc.D * 1000:g} mm; "
          f"row at y = {Y_ROW} m")
    print("law    cell     rate kg/s/m   C      centre speed 0.3/0.5 s   zone width 0.3/0.5 s"
          "   mass      run")
    for r in res:
        print(f"{r['law']:6s} {r['h'] * 1000:4.1f}mm  {r['rate']:8.1f}   {r['C']:.3f}   "
              f"{r['vc0.3']:.3f} / {r['vc0.5']:.3f} m/s      "
              f"{r['w0.3'] * 1000:5.1f} / {r['w0.5'] * 1000:5.1f} mm    "
              f"{r['cons']:.0e}  {r['sec']:6.0f} s")
    print(f"total {time.perf_counter() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
