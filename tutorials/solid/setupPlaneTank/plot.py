"""Figures for the setupPlaneTank tutorial (run after ./Allrun; numpy,
matplotlib, OpenFOAM sourced for the cell centres).

snapshots.png  material (alpha > 0.5) coloured by speed at six times, the half
               tank mirrored about its symmetry plane
mass.png       material in the tank against time, steady rate fitted
rate.txt       the steady rate and the Beverloo constant (plane flow)
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
CASE = os.path.join(HERE, "case")
T_OPEN, D, RHO, G = 0.5, 0.1, 1500.0, 9.81
INK, MUTED, SERIES, REF = "#1f1f1e", "#6b6b68", "#2a78d6", "#b8322a"


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


def times():
    out = []
    for d in os.listdir(CASE):
        try:
            if float(d) > 0:
                out.append((float(d), d))
        except ValueError:
            pass
    return sorted(out)


def mass_curve():
    log = open(f"{CASE}/log.celFoam").read()
    blocks = re.findall(r"Time = ([0-9.e+-]+)  deltaT[^\n]*\n(?:(?!Time = ).)*?"
                        r"Total mass\s+= ([0-9.e+-]+)", log, re.S)
    t = np.array([float(a) for a, _ in blocks])
    m = np.array([float(b) for _, b in blocks])
    t, i = np.unique(t, return_index=True)          # resumed runs repeat times
    return t, m[i]


def main():
    if not os.path.exists(f"{CASE}/0/C"):
        subprocess.run(["postProcess", "-func", "writeCellCentres", "-time", "0"],
                       cwd=CASE, check=True, stdout=subprocess.DEVNULL)
    C = read_internal(f"{CASE}/0/C", 10 ** 9)
    n = len(C)
    x, y, zc = C[:, 0], C[:, 1], C[:, 2]
    # plane case, one cell deep: mass per metre of depth = mass / depth
    pts = open(f"{CASE}/constant/polyMesh/points").read()
    zs = np.array(re.findall(r"\(\S+ \S+ (\S+)\)", pts), float)
    depth = zs.max() - zs.min()

    t, m = mass_curve()
    m = 2 * m / depth                                  # whole tank, per m of depth
    m_open = np.interp(T_OPEN, t, m)
    sel = (t > T_OPEN) & (m <= 0.97 * m_open) & (m >= 0.85 * m_open)
    rate = -np.polyfit(t[sel], m[sel], 1)[0] if sel.sum() >= 3 else float("nan")
    cb = rate / (RHO * math.sqrt(G) * D ** 1.5)
    msg = (f"steady rate {rate:.1f} kg/s per m of depth (whole tank), "
           f"Beverloo C = {cb:.2f} (plane flow, W = {D} m)")
    open(f"{HERE}/rate.txt", "w").write(msg + "\n")
    print(msg)

    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot(t, m, color=SERIES, lw=1.4)
    ax.axvline(T_OPEN, color=MUTED, lw=0.8, ls="--")
    ax.text(T_OPEN, ax.get_ylim()[0], " outlet opens", color=MUTED, fontsize=8, va="bottom")
    if sel.sum() >= 3:
        tt = np.array([t[sel][0], t[sel][-1]])
        a = np.polyfit(t[sel], m[sel], 1)
        ax.plot(tt, np.polyval(a, tt), color=REF, lw=2.5, alpha=0.5,
                label=f"fit: {rate:.0f} kg/s/m")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mass in tank (kg per m depth)")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{HERE}/mass.png", dpi=150)

    ts = times()
    want = [T_OPEN, 0.7, 1.0, 1.3, 1.6, ts[-1][0]]
    snaps = []
    for w in want:
        k = int(np.argmin([abs(a - w) for a, _ in ts]))
        if ts[k] not in snaps:
            snaps.append(ts[k])
    xx = np.r_[x, -x]
    yy = np.r_[y, y]
    tri = mtri.Triangulation(xx, yy)
    fig, axs = plt.subplots(1, len(snaps), figsize=(1.6 * len(snaps) + 0.9, 3.0),
                            constrained_layout=True, sharey=True)
    for ax, (tv, d) in zip(np.atleast_1d(axs), snaps):
        al = read_internal(f"{CASE}/{d}/alpha", n)
        U = read_internal(f"{CASE}/{d}/U", n)
        sp_ = np.hypot(U[:, 0], U[:, 1])
        a2, s2 = np.r_[al, al], np.r_[sp_, sp_]
        tri.set_mask(np.min(a2[tri.triangles], axis=1) < 0.5)
        im = ax.tripcolor(tri, s2, shading="gouraud", cmap="inferno", vmin=0, vmax=1.5)
        ax.set_title(f"t = {tv:g} s", fontsize=8)
        ax.set_aspect("equal")
        ax.set_xlim(xx.min(), xx.max())
        ax.set_ylim(yy.min(), yy.max())
        ax.set_xticks([])
        ax.set_facecolor("#ececea")
    np.atleast_1d(axs)[0].set_ylabel("y (m)")
    cbar = fig.colorbar(im, ax=list(np.atleast_1d(axs)), shrink=0.8)
    cbar.set_label("speed (m/s)")
    fig.savefig(f"{HERE}/snapshots.png", dpi=150)
    print("snapshots.png, mass.png, rate.txt written")


if __name__ == "__main__":
    main()
