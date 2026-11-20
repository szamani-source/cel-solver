"""Discharge snapshots of a Couto silo run: bulk density (alpha x rho, the
material per cell volume) and speed, cell by cell on the mesh, mirrored
about the axis so the silo reads whole (Sia, 2026-11-20).

    python3 plot_snapshots.py CASE OUT.png "title" t1 t2 t3 ...
"""
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "angled"))
import hopper_case as hc  # noqa: E402
import hopper_pressure_checks as hp  # noqa: E402
import buffered_case as bc  # noqa: E402
import angled_checks as ac  # noqa: E402


def mesh_polys(p):
    m = f"{p}/constant/polyMesh"
    pts = ac.foam_list(f"{m}/points", "point")
    faces = ac.foam_list(f"{m}/faces", "face")
    t = math.radians(hc.THETA)
    polys = [np.column_stack([pts[f][:, 0] / math.cos(t), pts[f][:, 1]])
             for f in faces if np.all(pts[f][:, 2] <= 1e-12) and len(f) == 4]
    return polys, np.array([q.mean(0) for q in polys])


def match(cen, a):
    """Row of a (cells file) for each polygon centre: nearest cell centre."""
    t = math.radians(hc.THETA)
    P = np.column_stack([a[:, 0] / math.cos(t), a[:, 1]])
    idx = np.empty(len(cen), int)
    for s in range(0, len(cen), 500):
        d = ((cen[s:s + 500, None, :] - P[None, :, :]) ** 2).sum(-1)
        idx[s:s + 500] = d.argmin(1)
    return idx


def main():
    p, out, title = sys.argv[1], sys.argv[2], sys.argv[3]
    times = [float(x) for x in sys.argv[4:]]
    c = dict(hc.COUTO)
    polys, cen = mesh_polys(p)
    mirror = [q * np.array([-1.0, 1.0]) for q in polys]
    shell = bc.shell_polygon(c, y_top=2.55)
    fig, axs = plt.subplots(2, len(times), figsize=(3.2 * len(times) + 1.5, 13))
    axs = np.atleast_2d(axs).reshape(2, len(times))
    for j, t in enumerate(times):
        a = hp.read_cells(p, f"cells_{t:g}")
        k = match(cen, a)
        rho = a[k, 2] * a[k, 3] * (a[k, 10] > 0)
        spd = np.hypot(a[k, 4], a[k, 5]) * (a[k, 2] > 1e-3)
        for i, (val, cmap, lim, lab) in enumerate(
                ((rho, "cividis", (700, 880), "bulk density (kg/m³); void and below 700 dark"),
                 (spd, "inferno", (0, 2.0), "speed (m/s)"))):
            ax = axs[i, j]
            for pp in (polys, mirror):
                pc = PolyCollection(pp, array=val, cmap=cmap, clim=lim, edgecolors="none",
                                    antialiased=False)
                ax.add_collection(pc)
            for sg in (1, -1):
                sp = np.vstack([shell, shell[:1]]) * np.array([sg, 1.0])
                ax.fill(sp[:, 0], sp[:, 1], color="#888888", lw=0)
            ax.set_xlim(-0.56, 0.56)
            ax.set_ylim(-0.05, 2.55)
            ax.set_aspect("equal")
            ax.set_xticks([-0.5, 0, 0.5])
            if j == 0:
                ax.set_ylabel("y (m)")
            else:
                ax.set_yticklabels([])
            ax.set_title(f"t = {t:g} s ({t - hp.T_OPEN:.2f} s after opening)" if i == 0
                         else "", fontsize=9)
            if j == len(times) - 1:
                fig.colorbar(pc, ax=axs[i, :], shrink=0.8, label=lab)
    fig.suptitle(title, fontsize=12)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
