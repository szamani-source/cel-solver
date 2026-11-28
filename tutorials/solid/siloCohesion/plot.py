"""Pictures for the siloCohesion tutorial (WP 3.5).

Run after ./Allrun (needs numpy and matplotlib; OpenFOAM sourced for the
cell centres):

    python3 plot.py

Reads every D<width> folder and writes snapshots.png (material and speed at
three times for each outlet width), mass.png (mass in the tank against time)
and result.txt (discharged fraction and the rate over the last 0.2 s against
the Beverloo rate with C = 1, and Jenike's critical width).
"""

import math
import os
import re
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "siloDischarge"))
from plot import read_internal  # noqa: E402

T_OPEN, RHO, G, W = 0.3, 1800.0, 9.81, 0.4
SNAPS = [0.3, 0.5, 1.3]
INK, MUTED, SURFACE = "#1f1f1e", "#6b6b68", "#fcfcfb"
SERIES = ["#2a78d6", "#d6602a", "#3a9a5b", "#8a5ad6"]


def time_dirs(case):
    out = {}
    for d in os.listdir(case):
        try:
            out[float(d)] = d
        except ValueError:
            pass
    return out


def history(case, depth):
    log = open(f"{case}/log.celFoam").read()
    blocks = re.findall(r"Time = ([0-9.e+-]+)  deltaT[^\n]*\n(?:(?!Time = ).)*?"
                        r"Total mass\s+= ([0-9.e+-]+)", log, re.S)
    m0 = float(re.search(r"Total mass\s+= ([0-9.e+-]+)", log).group(1)) / depth
    t = np.r_[0.0, [float(a) for a, _ in blocks]]
    m = np.r_[m0, [float(b) for _, b in blocks]] / np.r_[1.0, np.full(len(blocks), depth)]
    return t, m


def main():
    cases = sorted((int(d[1:]), f"{HERE}/{d}") for d in os.listdir(HERE)
                   if re.fullmatch(r"D\d+", d) and os.path.exists(f"{HERE}/{d}/log.celFoam"))
    if not cases:
        sys.exit("no D<width> folders with log.celFoam: run ./Allrun first")
    mp = open(f"{cases[0][1]}/constant/materialProperties").read()
    c = float(re.search(r"^c\s+([0-9.e+-]+);", mp, re.M).group(1))
    phi = float(re.search(r"^phi\s+([0-9.e+-]+);", mp, re.M).group(1))
    E = float(re.search(r"^E\s+([0-9.e+-]+);", mp, re.M).group(1))
    sp = math.sin(math.radians(phi))
    sigma_c = 2 * c * math.cos(math.radians(phi)) / (1 - sp)
    B_c = sigma_c / (RHO * G)

    plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED,
                         "axes.edgecolor": MUTED, "figure.facecolor": SURFACE,
                         "axes.facecolor": SURFACE})
    fig, axs = plt.subplots(len(cases), len(SNAPS), figsize=(2.2 * len(SNAPS) + 1.2,
                            3.3 * len(cases)), sharey=True, squeeze=False)
    lines = [f"c = {c:g} Pa, phi {phi:g} deg, E = {E / 1e6:g} MPa: unconfined strength "
             f"{sigma_c:.0f} Pa, Jenike critical width sigma_c / (rho g) = {B_c * 1e3:.0f} mm"]
    mass = []
    im = None
    for D_mm, case in cases:
        D = D_mm / 1000
        if not os.path.exists(f"{case}/0/C"):
            subprocess.run(["postProcess", "-func", "writeCellCentres", "-time", "0"],
                           cwd=case, check=True, stdout=subprocess.DEVNULL)
        C = read_internal(f"{case}/0/C", 10 ** 9)
        n = len(C)
        x, y = C[:, 0], C[:, 1]
        h = np.diff(np.unique(np.round(x, 9)))[0]
        t, m = history(case, h)
        i_open = int(np.argmin(abs(t - T_OPEN)))
        sel = t >= t[-1] - 0.2 - 1e-9
        rate = -np.polyfit(t[sel], m[sel], 1)[0]
        bev = RHO * math.sqrt(G) * D ** 1.5
        state = "stopped" if rate < 0.05 * bev else "flows"
        lines.append(f"D = {D_mm} mm: {state}; discharged {1 - m[-1] / m[i_open]:.3f} of the "
                     f"fill in {t[-1] - T_OPEN:g} s, rate over the last 0.2 s {rate:.1f} kg/s/m "
                     f"= {rate / bev:.3f} of Beverloo C = 1")
        mass.append((D_mm, t, m, state))
        dirs = time_dirs(case)
        nx, ny = int(round(W / h)), int(round(0.62 / h))
        ix = np.floor(x / h + 1e-6).astype(int)
        iy = np.floor(y / h + 1e-6).astype(int)
        row = axs[[d for d, _ in cases].index(D_mm)]
        for ax, ts in zip(row, SNAPS):
            key = min(dirs, key=lambda k: abs(k - ts))
            a = read_internal(f"{case}/{dirs[key]}/alpha", n)
            U = read_internal(f"{case}/{dirs[key]}/U", n)
            grid = np.full((ny, nx), np.nan)
            grid[iy, ix] = np.where(a > 0.5, np.hypot(U[:, 0], U[:, 1]), np.nan)
            cmap = matplotlib.colormaps["Blues"].copy()
            cmap.set_bad("#d6d6d3")
            im = ax.imshow(grid, origin="lower", extent=(0, W, 0, 0.62), cmap=cmap,
                           vmin=0, vmax=1.5, interpolation="nearest")
            ax.plot([0, (W - D) / 2], [0, 0], color=INK, lw=2)
            ax.plot([(W + D) / 2, W], [0, 0], color=INK, lw=2)
            ax.set_title(f"D = {D_mm} mm, t = {key:g} s"
                         + (" (opens)" if abs(key - T_OPEN) < 1e-9 else ""), fontsize=9)
            ax.set_xticks([0, 0.2, 0.4])
        row[0].set_ylabel("y (m)")
    for ax in axs[-1]:
        ax.set_xlabel("x (m)")
    cb = fig.colorbar(im, ax=axs, shrink=0.6, pad=0.01)
    cb.set_label("speed (m/s); grey = empty")
    fig.suptitle(f"Tank with cohesion c = {c:g} Pa (Jenike critical width "
                 f"{B_c * 1e3:.0f} mm), E = {E / 1e6:g} MPa", fontsize=10)
    fig.savefig(f"{HERE}/snapshots.png", dpi=150, bbox_inches="tight")

    fig, ax = plt.subplots(figsize=(5, 3.2))
    for k, (D_mm, t, m, state) in enumerate(mass):
        ax.plot(t, m, color=SERIES[k % len(SERIES)], lw=2, label=f"D = {D_mm} mm ({state})")
    ax.axvline(T_OPEN, color=MUTED, lw=1, ls=":")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mass in tank (kg per m depth)")
    ax.set_title("Mass in the tank", fontsize=10)
    ax.grid(color="#e4e4e1", lw=0.6)
    ax.set_ylim(0, None)
    ax.legend(frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(f"{HERE}/mass.png", dpi=150, bbox_inches="tight")
    msg = "\n".join(lines) + "\n"
    open(f"{HERE}/result.txt", "w").write(msg)
    print(msg, end="")


if __name__ == "__main__":
    main()
