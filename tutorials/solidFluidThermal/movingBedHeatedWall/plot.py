"""Plots for the movingBedHeatedWall tutorial (WP TU.3, decision 0136).

Run after ./Allrun (needs numpy and matplotlib):

    python3 plot.py

Writes temperature.png (the bed temperature at the end, the gap stretched
across) and nusselt.png (local Nusselt number along the heated plate against
the slug-flow solution with the near-wall gap, Adapa et al. 2024, eqs. 4 and
9), and results.txt (mean h against eq. 9).
"""

import math
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

# The case (as constant/ and system/blockMeshDict)
B, LEN = 2.5e-3, 0.5            # half gap, heated length [m]
NX, NY = 200, 10
K_EFF, D_AIR, T_BED = 0.28, 35e-6, 753.15   # Adapa et al. Table 5 (480 C, 10 mm/s)
RHO, CP, U_REAL = 2090.0, 1100.0, 0.010     # bulk density, specific heat, real speed
Q, END = 5000.0, 0.7
DH = 4 * B


def air_conductivity(T):        # Sutherland form for air
    return 0.0241 * (T / 273.15) ** 1.5 * (273.15 + 194.0) / (T + 194.0)


def nu_local(gz, gap, terms=2000):
    """Local Nu_Dh of slug flow between plates, uniform flux, with the gap
    resistance R_nw k_eff/(4 b) (Adapa et al. 2024, eq. 4 with eq. 2)."""
    s = sum(math.exp(-16.0 * (n * math.pi) ** 2 / gz) / (n * math.pi) ** 2
            for n in range(1, terms + 1))
    return 1.0 / (1.0 / 12.0 - 0.5 * s + gap)


def nu_mean(dh_over_l, pe, gap, terms=2000):
    """Mean Nu_Dh over the heated length (Adapa et al. 2024, eq. 9)."""
    s = sum(dh_over_l * pe / (32.0 * (n * math.pi) ** 4)
            * (math.exp(-16.0 * (n * math.pi) ** 2 / (dh_over_l * pe)) - 1.0)
            for n in range(1, terms + 1))
    return 1.0 / (1.0 / 12.0 + s + gap)


def main():
    dy = B / NY
    r_nw = D_AIR / air_conductivity(T_BED)
    gap = r_nw * K_EFF / (4 * B)
    pe = U_REAL * DH * RHO * CP / K_EFF
    T = read_internal(f"{HERE}/{END:g}/T", NX * NY).reshape(NY, NX)
    x = (np.arange(NX) + 0.5) * LEN / NX
    Tm = T.mean(axis=0)
    Tw = T[0] + Q * (dy / 2) / K_EFF + Q * r_nw    # plate: wall cell + half cell + gap
    nu = Q / (Tw - Tm) * DH / K_EFF
    want = np.array([nu_local(DH / xi * pe, gap) for xi in x])
    h = Q / (Tw - Tm).mean()
    hw = nu_mean(DH / LEN, pe, gap) * K_EFF / DH
    sel = x > 0.05
    worst = np.abs(nu[sel] / want[sel] - 1).max()

    fig, ax = plt.subplots(figsize=(8, 2.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    im = ax.imshow(T - 273.15, origin="lower", aspect="auto", cmap="inferno",
                   extent=(0, LEN * 1e3, 0, B * 1e3))
    ax.set_xlabel("distance along the plate [mm] (flow left to right)", color=INK)
    ax.set_ylabel("from the plate [mm]", color=INK)
    ax.tick_params(colors=MUTED, labelsize=9)
    cb = fig.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("T [C]", color=INK)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_title("Bed temperature at 0.7 s: heated plate at the bottom, mid-plane at the top",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(f"{HERE}/temperature.png")

    fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    ax.plot(x * 1e3, want, color=INK, lw=1.2, ls="--",
            label="slug flow with the wall gap (Adapa et al. 2024)")
    ax.plot(x * 1e3, nu, color=SERIES, lw=1.8, label="celFoam")
    ax.set_xlabel("distance along the plate x [mm]")
    ax.set_ylabel(r"local Nu$_{D_h}$ = q D$_h$ / (k (T$_w$ - T$_m$))")
    ax.set_xlim(0, LEN * 1e3)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(f"Mean h {h:.1f} W/(m$^2$ K) against {hw:.1f} (eq. 9)",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(f"{HERE}/nusselt.png")
    msg = (f"Peclet number {pe:.0f}, near-wall gap R_nw {r_nw:.3e} m^2 K/W\n"
           f"mean h {h:.1f} W/(m^2 K) against {hw:.1f} (eq. 9): {(h / hw - 1) * 100:+.2f} %\n"
           f"local Nu past x = 50 mm within {worst * 100:.2f} % of eq. 4\n")
    open(f"{HERE}/results.txt", "w").write(msg)
    print(msg, end="")


if __name__ == "__main__":
    main()
