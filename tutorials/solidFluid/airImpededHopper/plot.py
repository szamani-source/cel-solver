"""Plots for the airImpededHopper tutorial (WP F.9; case of WP F.5,
decision 0130).

Run after ./Allrun (needs numpy and matplotlib; OpenFOAM sourced for the
cell centres):

    python3 plot.py

Writes discharge.png (mass discharged per metre of slot, with gas and dry),
snapshots.png (grain speed and gas pressure at four times) and
discharge.txt (the steady rates).
"""

import math
import os
import re
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
R, RO, BETA = 0.10, 0.03, 30.0
HT = (R - RO) / math.tan(math.radians(BETA))
DEPTH = 0.0125                          # mesh thickness [m] (one cell)
T_OPEN, WINDOW = 0.15, (0.3, 0.5)
SNAPS = [0.15, 0.25, 0.35, 0.45]
INK, MUTED, SURFACE, GAS, DRY = "#1f1f1e", "#6b6b68", "#fcfcfb", "#c2410c", "#2a78d6"


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


def discharge(case):
    """(t, discharged mass per metre of slot) at each write."""
    log = open(f"{case}/log.celFoam").read()
    t, m = [0.0], [0.0]
    for chunk in log.split("\nTime = ")[1:]:
        k = re.search(r"^Discharged mass = (\S+)", chunk, re.M)
        if k:
            t.append(float(chunk.split()[0]))
            m.append(float(k.group(1)) / DEPTH)
    return np.array(t), np.array(m)


def rate(t, m):
    s = (t >= WINDOW[0] - 1e-9) & (t <= WINDOW[1] + 1e-9)
    return float(np.polyfit(t[s], m[s], 1)[0])


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def main():
    if not os.path.exists(f"{HERE}/0/C"):
        subprocess.run(["postProcess", "-func", "writeCellCentres", "-time", "0"],
                       cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    C = read_internal(f"{HERE}/0/C", 10 ** 9)
    n = len(C)
    x, y = C[:, 0], C[:, 1]

    # discharge, with gas and dry
    tg, mg = discharge(HERE)
    td, md = discharge(f"{HERE}/dry")
    wg, wd = rate(tg, mg), rate(td, md)
    with open(f"{HERE}/discharge.txt", "w") as f:
        f.write(f"# steady rate {WINDOW[0]}-{WINDOW[1]} s [kg/(m s)]\n")
        f.write(f"with gas {wg:.2f}\ndry {wd:.2f}\nratio {wg / wd:.3f}\n")
    fig, ax = plt.subplots(figsize=(6.4, 3.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ax.plot(td, md, color=DRY, lw=2, label=f"dry, {wd:.1f} kg/(m s)")
    ax.plot(tg, mg, color=GAS, lw=2, label=f"with gas, {wg:.1f} kg/(m s)")
    ax.axvline(T_OPEN, color=MUTED, lw=0.8, ls=":")
    ax.text(T_OPEN + 0.005, ax.get_ylim()[1] * 0.92, "outlet opens", color=MUTED, fontsize=8)
    ax.set_xlabel("time [s]", color=INK)
    ax.set_ylabel("discharged mass [kg per m of slot]", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(f"{HERE}/discharge.png", facecolor=SURFACE)

    # snapshots: material fraction and gas pressure
    tri = mtri.Triangulation(x, y)
    xm, ym = x[tri.triangles].mean(axis=1), y[tri.triangles].mean(axis=1)
    half = np.where(ym < HT, RO + (R - RO) * ym / HT, R)
    tri.set_mask(np.abs(xm) > half)
    times = {float(d): d for d in os.listdir(HERE) if re.fullmatch(r"[0-9.e+-]+", d)}
    snaps = [times[min(times, key=lambda t: abs(t - s))] for s in SNAPS]
    pmax = max(np.abs(read_internal(f"{HERE}/{t}/p", n)).max() for t in snaps)
    umax = max(np.linalg.norm(read_internal(f"{HERE}/{t}/U", n), axis=1).max() for t in snaps)
    fig, axs = plt.subplots(2, len(snaps), figsize=(10, 6.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    for k, t in enumerate(snaps):
        a = read_internal(f"{HERE}/{t}/alpha", n)
        u = np.where(a > 0.5, np.linalg.norm(read_internal(f"{HERE}/{t}/U", n), axis=1), np.nan)
        p = read_internal(f"{HERE}/{t}/p", n)
        for row, (f, cmap, lo, hi) in enumerate(((u, "viridis", 0, umax), (p, "RdBu_r", -pmax, pmax))):
            ax = axs[row, k]
            im = ax.tripcolor(tri, f, cmap=cmap, vmin=lo, vmax=hi, shading="gouraud")
            ax.set_aspect("equal")
            ax.axis("off")
            if row == 0:
                ax.set_title(f"t = {float(t):.2f} s", color=INK, fontsize=10)
        if k == len(snaps) - 1:
            fig.colorbar(axs[0, k].collections[0], ax=axs[0, :], shrink=0.8, label="grain speed |U| [m/s]")
            fig.colorbar(im, ax=axs[1, :], shrink=0.8, label="gas pressure p - p0 [Pa]")
    fig.savefig(f"{HERE}/snapshots.png", facecolor=SURFACE)
    print(f"with gas {wg:.2f}, dry {wd:.2f} kg/(m s), ratio {wg / wd:.3f}")


if __name__ == "__main__":
    main()
