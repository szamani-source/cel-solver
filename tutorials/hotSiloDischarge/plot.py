"""Plots for the hotSiloDischarge tutorial (WP T.8, decision 0135).

Run after ./Allrun (needs numpy and matplotlib; OpenFOAM sourced for the
cell centres):

    python3 plot.py

Writes snapshots.png (temperature of the material at six times) and
outlet.png (temperature of the sand leaving through the outlet, and the
share of the heat still in the tank), and outlet.txt (the same numbers).
The outlet temperature between two writes is the heat carried out over
c times the mass discharged (both printed in log.celFoam at each write).
"""

import os
import re
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
T_OPEN, C_SAND, T_HOT, T_COLD = 0.3, 800.0, 423.15, 293.15
SNAPS = [0.3, 0.7, 1.1, 1.5, 1.9, 2.3]
INK, MUTED, SURFACE, SERIES, SERIES2 = "#1f1f1e", "#6b6b68", "#fcfcfb", "#c2410c", "#2a78d6"


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

    # counters printed at every write: time, discharged mass, heat carried
    # in, thermal energy (one block per write time)
    log = open(f"{HERE}/log.celFoam").read()
    e0 = float(re.search(r"^Thermal energy at start\s*=\s*(\S+)", log, re.M).group(1))
    t, out, carried, left = [], [], [], []
    for chunk in log.split("\nTime = ")[1:]:
        if "Heat carried in" not in chunk:
            continue
        val = lambda k: float(re.search(rf"^{k}\s*=\s*(\S+)", chunk, re.M).group(1))
        t.append(float(chunk.split()[0]))
        out.append(val("Discharged mass"))
        carried.append(val("Heat carried in"))
        left.append(val("Thermal energy") / e0)
    t, out, carried, left = map(np.array, (t, out, carried, left))
    dm = np.diff(np.r_[0.0, out])
    dq = -np.diff(np.r_[0.0, carried])
    with np.errstate(invalid="ignore", divide="ignore"):
        T_out = np.where(dm > 1e-9, dq / (C_SAND * dm), np.nan)
    tm = t - 0.5 * np.diff(np.r_[0.0, t])
    rows = [f"{a:.3f} {b - 273.15:.1f} {c:.3f}" for a, b, c in zip(tm, T_out, left)]
    open(f"{HERE}/outlet.txt", "w").write(
        "# mid time (s), outlet temperature (C), share of the heat left in the tank\n"
        + "\n".join(rows) + "\n")

    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.edgecolor": MUTED, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE})

    dirs = dict(time_dirs())
    fig, axs = plt.subplots(1, len(SNAPS), figsize=(12, 3.6), sharey=True)
    nx, ny = int(round(0.4 / h)), int(round(0.62 / h))
    ix = np.floor(x / h + 1e-6).astype(int)
    iy = np.floor(y / h + 1e-6).astype(int)
    im = None
    for ax, ts in zip(axs, SNAPS):
        key = min(dirs, key=lambda k: abs(k - ts))
        a = read_internal(f"{HERE}/{dirs[key]}/alpha", n)
        T = read_internal(f"{HERE}/{dirs[key]}/T", n)
        grid = np.full((ny, nx), np.nan)
        grid[iy, ix] = np.where(a > 0.5, T - 273.15, np.nan)
        cmap = matplotlib.colormaps["inferno"].copy()
        cmap.set_bad("#d6d6d3")
        im = ax.imshow(grid, origin="lower", extent=(0, 0.4, 0, 0.62), cmap=cmap,
                       vmin=T_COLD - 273.15, vmax=T_HOT - 273.15, interpolation="nearest")
        ax.plot([0, 0.17], [0, 0], color=INK, lw=2)
        ax.plot([0.23, 0.4], [0, 0], color=INK, lw=2)
        ax.set_title(f"t = {key:g} s" + (" (opens)" if abs(key - T_OPEN) < 1e-9 else ""),
                     fontsize=9)
        ax.set_xlabel("x (m)")
        ax.set_xticks([0, 0.2, 0.4])
    axs[0].set_ylabel("y (m)")
    cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.01)
    cb.set_label("sand temperature (C); grey = empty")
    fig.suptitle("Hot sand under cold sand, discharging", fontsize=10)
    fig.savefig(f"{HERE}/snapshots.png", dpi=150, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    sel = np.isfinite(T_out)
    ax.plot(tm[sel], T_out[sel] - 273.15, color=SERIES, lw=2, label="leaving the outlet (C)")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("outlet temperature (C)")
    ax.set_ylim(T_COLD - 273.15 - 5, T_HOT - 273.15 + 5)
    ax2 = ax.twinx()
    ax2.plot(t, 100 * left, color=SERIES2, lw=2, ls="--", label="heat left in the tank (%)")
    ax2.set_ylabel("heat left in the tank (%)")
    ax2.set_ylim(0, 105)
    ax.axvline(T_OPEN, color=MUTED, lw=1, ls=":")
    lines = ax.get_lines()[:1] + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], frameon=False, loc="lower left")
    ax.grid(color="#e4e4e1", lw=0.6)
    for a in (ax, ax2):
        a.spines["top"].set_visible(False)
    ax.set_title("Funnel flow draws the cold top layer through the hot charge", fontsize=10)
    fig.savefig(f"{HERE}/outlet.png", dpi=150, bbox_inches="tight")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
