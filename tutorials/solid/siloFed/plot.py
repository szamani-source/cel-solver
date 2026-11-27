"""Plots and rates for the siloFed tutorial (WP 3.7b), from siloDischarge.

Run after ./Allrun (needs numpy and matplotlib; OpenFOAM sourced for the
cell centres):

    python3 plot.py

Writes snapshots.png (material and speed at six times), mass.png (mass in
the tank, fed and discharged against time) and rate.txt (steady rates out
and in). Rates are least-squares slopes from 0.35 to 2.0 s after opening,
as tests/fedsilo/fed_silo_checks.py.
"""

import math
import os
import re
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
T_OPEN, D, RHO, G, DEPTH = 0.3, 0.06, 1800.0, 9.81, None
SNAPS = [0.3, 0.6, 1.0, 1.4, 1.8, 2.3]
INK, MUTED, SURFACE, SERIES = "#1f1f1e", "#6b6b68", "#fcfcfb", "#2a78d6"


def read_internal(path, n):
    txt = open(path).read()
    m = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    if m:
        v = m.group(1).strip("() ").split()
        return np.tile(np.array(v, float), (n, 1)).squeeze()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*\n?(\d+)\s*\n?\(", txt)
    body = txt[m.end():]
    rows = re.findall(r"\(([^()]*)\)|([-+0-9.eE]+)", body[:body.index("\n)")])
    if rows and rows[0][0]:
        return np.array([r[0].split() for r in rows[:n]], float)
    return np.array([r[1] for r in rows[:n]], float)


def time_dirs():
    out = []
    for d in os.listdir(HERE):
        try:
            out.append((float(d), d))
        except ValueError:
            pass
    return sorted(out)


def main():
    if not os.path.exists(f"{HERE}/0/C"):
        subprocess.run(["postProcess", "-func", "writeCellCentres", "-time", "0"],
                       cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    C = read_internal(f"{HERE}/0/C", 10 ** 9)
    n = len(C)
    x, y = C[:, 0], C[:, 1]
    h = np.diff(np.unique(np.round(x, 9)))[0]
    depth = h                                   # one cell thick (2D)

    # tank, discharged and fed mass at every write time from the log
    log = open(f"{HERE}/log.celFoam").read()
    blocks = re.findall(r"Time = ([0-9.e+-]+)  deltaT[^\n]*\n(?:(?!Time = ).)*?"
                        r"Total mass\s+= ([0-9.e+-]+).*?Discharged mass = ([0-9.e+-]+)"
                        r"\s*\nFed mass\s+= ([0-9.e+-]+)", log, re.S)
    t = np.array([float(b[0]) for b in blocks])
    m, out, fed = (np.array([float(b[i]) for b in blocks]) / depth for i in (1, 2, 3))
    m0 = float(re.search(r"Total mass\s+= ([0-9.e+-]+)", log).group(1)) / depth
    t, m, out, fed = np.r_[0.0, t], np.r_[m0, m], np.r_[0.0, out], np.r_[0.0, fed]
    E = float(re.search(r"^E\s+([0-9.e+-]+);", open(f"{HERE}/constant/materialProperties")
                        .read(), re.M).group(1))
    i_open = int(np.argmin(abs(t - T_OPEN)))
    m_open = m[i_open]
    sel = (t >= T_OPEN + 0.35) & (t <= T_OPEN + 2.0)
    r_out = np.polyfit(t[sel], out[sel], 1)[0]
    r_in = np.polyfit(t[sel], fed[sel], 1)[0]
    bal = np.abs((m - fed + out) / m0 - 1).max()
    rate = r_out
    Cb = rate / (RHO * math.sqrt(G) * D ** 1.5)
    msg = (f"E = {E / 1e6:g} MPa, {h * 1000:g} mm cells: out {r_out:.1f} kg/s/m "
           f"(Beverloo C = {Cb:.2f}), fed {r_in:.1f} kg/s/m; tank mass at opening "
           f"{m_open:.1f} kg/m, at the end {m[-1] / m_open:.3f} of it; mass balance "
           f"{bal:.0e}\n")
    open(f"{HERE}/rate.txt", "w").write(msg)
    print(msg, end="")

    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.edgecolor": MUTED, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE})

    # snapshots: material (alpha) shaded by speed
    dirs = dict(time_dirs())
    fig, axs = plt.subplots(1, len(SNAPS), figsize=(12, 3.6), sharey=True)
    nx, ny = int(round(0.4 / h)), int(round(0.7 / h))
    ix = np.floor(x / h + 1e-6).astype(int)
    iy = np.floor(y / h + 1e-6).astype(int)
    im = None
    for ax, ts in zip(axs, SNAPS):
        key = min(dirs, key=lambda k: abs(k - ts))
        a = read_internal(f"{HERE}/{dirs[key]}/alpha", n)
        U = read_internal(f"{HERE}/{dirs[key]}/U", n)
        sp = np.hypot(U[:, 0], U[:, 1])
        grid = np.full((ny, nx), np.nan)
        grid[iy, ix] = np.where(a > 0.5, sp, np.nan)
        cmap = matplotlib.colormaps["Blues"].copy()
        cmap.set_bad("#d6d6d3")
        im = ax.imshow(grid, origin="lower", extent=(0, 0.4, 0, 0.7), cmap=cmap,
                       vmin=0, vmax=1.5, interpolation="nearest")
        ax.plot([0, 0.17], [0, 0], color=INK, lw=2)
        ax.plot([0.23, 0.4], [0, 0], color=INK, lw=2)
        ax.plot([0, 0.17], [0.7, 0.7], color=INK, lw=2)
        ax.plot([0.23, 0.4], [0.7, 0.7], color=INK, lw=2)
        ax.set_title(f"t = {key:g} s" + (" (opens)" if abs(key - T_OPEN) < 1e-9 else ""),
                     fontsize=9)
        ax.set_xlabel("x (m)")
        ax.set_xticks([0, 0.2, 0.4])
    axs[0].set_ylabel("y (m)")
    cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.01)
    cb.set_label("speed (m/s); grey = empty")
    fig.suptitle(f"Fed tank discharge, E = {E / 1e6:g} MPa, {h * 1000:g} mm cells", fontsize=10)
    fig.savefig(f"{HERE}/snapshots.png", dpi=150, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.plot(t, m, color=SERIES, lw=2, label="in the tank")
    ax.plot(t, fed, color="#e07b39", lw=1.5, label="fed in")
    ax.plot(t, out, color=MUTED, lw=1.5, ls="--", label="discharged")
    ax.axvline(T_OPEN, color=MUTED, lw=1, ls=":")
    ax.legend(frameon=False, loc="center right")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mass (kg per m depth)")
    ax.set_title(f"Out {r_out:.0f} kg/s/m, fed {r_in:.0f} kg/s/m", fontsize=10)
    ax.grid(color="#e4e4e1", lw=0.6)
    ax.set_ylim(0, None)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(f"{HERE}/mass.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
