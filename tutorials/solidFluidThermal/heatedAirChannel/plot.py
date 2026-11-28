"""Plots for the heatedAirChannel tutorial (WP TU.3, decision 0136).

Run after ./Allrun (needs numpy and matplotlib):

    python3 plot.py

Writes temperature.png (air temperature at the end) and bulk.png (the bulk,
flow-weighted temperature difference to the plates along the channel against
the developed decay exp(-lambda x) with Nu = 7.5407, Shah and London 1978),
and results.txt (Nu from the slope between x = 20 and 55 mm).
"""

import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, SURFACE, SERIES, SERIES2 = "#1f1f1e", "#6b6b68", "#fcfcfb", "#c2410c", "#2a78d6"


def read_internal(path, n):
    """Internal field of an ascii OpenFOAM field file (n values)."""
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


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.grid(color="#e4e4e0", lw=0.6)

# The case (as system/blockMeshDict, initFields.py and 0.orig/Tg)
H, LC = 2.0e-3, 60.0e-3
NX, NY = 200, 16
K_G, C_G = 0.02569, 1006.0
T0, DT0, P0, END = 293.15, 10.0, 101325.0, 0.3
NU_DEV = 7.5407                 # developed Nu_Dh, plates at a fixed temperature


def main():
    Tg = read_internal(f"{HERE}/{END:g}/Tg", NX * NY).reshape(NY, NX)
    U = read_internal(f"{HERE}/{END:g}/Ug", NX * NY)[:, 0].reshape(NY, NX)
    p = read_internal(f"{HERE}/{END:g}/p", NX * NY).reshape(NY, NX)
    x = (np.arange(NX) + 0.5) * LC / NX
    th = ((Tg - T0) * U).sum(axis=0) / U.sum(axis=0) / DT0
    sel = (x > 0.020) & (x < 0.055)
    lam = -np.polyfit(x[sel], np.log(th[sel]), 1)[0]
    rho = (P0 + p.mean()) / (287.05 * T0)
    G = rho * U.mean(axis=0)[sel].mean()
    nu = lam * G * C_G * H / 2 * 2 * H / K_G
    lam_dev = NU_DEV * K_G / (2 * H) * 2 / (G * C_G * H)
    x0 = x[sel][0]
    th0 = th[sel][0]

    fig, ax = plt.subplots(figsize=(8, 2.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    im = ax.imshow(Tg - 273.15, origin="lower", aspect="auto", cmap="inferno",
                   extent=(0, LC * 1e3, 0, H * 1e3))
    ax.set_xlabel("distance from the inlet [mm] (flow left to right)", color=INK)
    ax.set_ylabel("height [mm]", color=INK)
    ax.tick_params(colors=MUTED, labelsize=9)
    cb = fig.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("Tg [C]", color=INK)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_title("Air temperature at 0.3 s: enters 10 K hot, plates held at 20 C",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(f"{HERE}/temperature.png")

    fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ax.semilogy(x * 1e3, th, color=SERIES, lw=1.8, label="celFoam, bulk temperature")
    ax.semilogy(x * 1e3, th0 * np.exp(-lam_dev * (x - x0)), color=INK, lw=1.2, ls="--",
                label=f"developed flow, Nu = {NU_DEV} (Shah and London 1978)")
    ax.axvspan(20, 55, color="#e4e4e0", alpha=0.5, lw=0)
    ax.set_xlabel("distance from the inlet x [mm]")
    ax.set_ylabel(r"$(T_b - T_{wall})/\Delta T$")
    ax.set_xlim(0, LC * 1e3)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(f"Nu from the slope (shaded): {nu:.3f} against {NU_DEV} ({(nu / NU_DEV - 1) * 100:+.2f} %)",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(f"{HERE}/bulk.png")
    msg = f"developed Nu_Dh {nu:.4f} against {NU_DEV} ({(nu / NU_DEV - 1) * 100:+.2f} %)\n"
    open(f"{HERE}/results.txt", "w").write(msg)
    print(msg, end="")


if __name__ == "__main__":
    main()
