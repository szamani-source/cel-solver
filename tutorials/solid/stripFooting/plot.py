"""Figures for the stripFooting tutorial (run after ./Allrun; numpy, matplotlib).

q-s.png       footing pressure q / cu against settlement s / B, with Prandtl's
              limit 2 + pi and the 10 % band of tests/test_footing.sh
snapshot.png  mobilised shear tau / cu in the clay at three settlements, with
              the footing outline
"""
import glob
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "run", "footing")
CU, B, NC = 10e3, 0.1, 2 + math.pi
SPEED, T_RAMP, X0, Y0 = 0.025, 0.05, 0.2525, 0.2975
INK, MUTED, SERIES, REF = "#1f1f1e", "#6b6b68", "#2a78d6", "#b8322a"


def settlement(t):
    return np.where(t < T_RAMP, 0.5 * SPEED * t ** 2 / T_RAMP,
                    0.5 * SPEED * T_RAMP + SPEED * (t - T_RAMP))


def q_s():
    s, q = np.loadtxt(os.path.join(CASE, "q.dat"), unpack=True)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.axhspan(0.9 * NC, 1.1 * NC, color=REF, alpha=0.10, lw=0)
    ax.axhline(NC, color=REF, lw=1.2, label="Prandtl 2 + $\\pi$ = 5.14 ($\\pm$10 %)")
    ax.plot(s, q, color=SERIES, lw=1.4, label="celFoam, 10 mm cells")
    ax.axvspan(0.10, 0.18, color=MUTED, alpha=0.08, lw=0)
    ax.text(0.14, 0.6, "averaging\nwindow", ha="center", color=MUTED, fontsize=8)
    ax.set_xlabel("settlement s / B")
    ax.set_ylabel("footing pressure q / c$_u$")
    ax.set_xlim(0, s.max())
    ax.set_ylim(0, 6.5)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "q-s.png"), dpi=150)


def snapshots():
    files = sorted(glob.glob(os.path.join(CASE, "cells_*.dat")),
                   key=lambda f: float(os.path.basename(f)[6:-4]))
    files = [f for f in files if not f.endswith(".proc0")]
    if not files:
        return
    fig, axs = plt.subplots(1, len(files), figsize=(3.0 * len(files) + 0.8, 2.6),
                            constrained_layout=True)
    axs = np.atleast_1d(axs)
    for ax, f in zip(axs, files):
        t = float(os.path.basename(f)[6:-4])
        d = np.loadtxt(f)
        x, y, al = d[:, 0], d[:, 1], d[:, 2]
        sxx, syy, sxy = d[:, 6], d[:, 7], d[:, 8]
        tau = np.sqrt((0.5 * (sxx - syy)) ** 2 + sxy ** 2) / CU
        h = np.diff(np.unique(np.round(x, 9)))[0]
        nx, ny = len(np.unique(np.round(x, 9))), len(np.unique(np.round(y, 9)))
        i = np.round((x - x.min()) / h).astype(int)
        j = np.round((y - y.min()) / h).astype(int)
        img = np.full((ny, nx), np.nan)
        img[j, i] = np.where(al > 0.5, tau, np.nan)
        im = ax.imshow(img, origin="lower", cmap="viridis", vmin=0, vmax=1.0,
                       extent=(x.min() - h / 2, x.max() + h / 2, y.min() - h / 2, y.max() + h / 2))
        s = float(settlement(np.array(t)))
        ax.add_patch(plt.Rectangle((X0, Y0 - s), B, 0.06, fc="#d9d9d6", ec=INK, lw=0.8))
        ax.set_xlim(0.05, 0.55)
        ax.set_ylim(0.05, 0.36)
        ax.set_aspect("equal")
        ax.set_title(f"t = {t:g} s, s / B = {s / B:.2f}", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    cb = fig.colorbar(im, ax=list(axs), shrink=0.8)
    cb.set_label("$\\tau$ / c$_u$")
    fig.savefig(os.path.join(HERE, "snapshot.png"), dpi=150)


if __name__ == "__main__":
    q_s()
    snapshots()
    print("q-s.png, snapshot.png written")
