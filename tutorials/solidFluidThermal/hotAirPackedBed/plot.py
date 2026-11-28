"""Plots for the hotAirPackedBed tutorial (WP TU.3, decision 0136).

Run after ./Allrun (needs numpy and matplotlib):

    python3 plot.py

Writes profiles.png (gas and grain temperature along the column at each
write against the closed form of Schumann (1929) with Gunn's (1978) heat
transfer coefficient) and results.txt (the largest differences).
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

# The case (as constant/ and 0.orig/)
L, N = 0.6, 300                 # column length [m], cells
E, D, U = 0.4, 0.02, 1.0        # gas fraction, grain diameter [m], superficial speed [m/s]
RHOP, C_S = 2650.0, 5.0         # grain density [kg/m^3], specific heat [J/(kg K)]
C_G, K_G = 1006.0, 0.02569      # air
T0, DT0, P0 = 293.15, 100.0, 101325.0


def air_viscosity(T):           # Sutherland
    return 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)


def gunn_nu(e, re, pr):         # Gunn (1978)
    return (7.0 - 10.0 * e + 5.0 * e * e) * (1.0 + 0.7 * re ** 0.2 * pr ** (1.0 / 3.0)) \
        + (1.33 - 2.4 * e + 1.2 * e * e) * re ** 0.7 * pr ** (1.0 / 3.0)


def bessel_i0(z):
    term, s, k = 1.0, 1.0, 0
    q = 0.25 * z * z
    while True:
        k += 1
        term *= q / (k * k)
        s += term
        if term < 1e-17 * s:
            return s


def schumann_J(x, y, n=1000):
    """J(x, y) = 1 - exp(-y) int_0^x exp(-s) I0(2 sqrt(y s)) ds (Simpson)."""
    if x <= 0.0:
        return 1.0
    h = x / n
    s = 0.0
    for i in range(n + 1):
        si = i * h
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        s += w * math.exp(-si - y) * bessel_i0(2.0 * math.sqrt(y * si))
    return 1.0 - s * h / 3.0


def main():
    mu = air_viscosity(T0)
    G = P0 / (287.05 * T0) * U                 # gas mass flux [kg/(m^2 s)]
    re, pr = G * D / mu, mu * C_G / K_G
    h = gunn_nu(E, re, pr) * K_G / D
    ha = h * 6 * (1 - E) / D                   # per bed volume [W/(m^3 K)]
    cs = (1 - E) * RHOP * C_S
    u_i = U / E
    x = (np.arange(N) + 0.5) * L / N
    times = sorted(float(d) for d in os.listdir(HERE)
                   if re_num(d) and float(d) > 0)

    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    worst_g = worst_s = 0.0
    xs = np.linspace(0, L, 61)[1:]
    for k, t in enumerate(times):
        Tg = read_internal(f"{HERE}/{t:g}/Tg", N)
        Ts = read_internal(f"{HERE}/{t:g}/T", N)
        tg, ts = (Tg - T0) / DT0, (Ts - T0) / DT0
        eg, es = [], []
        for xi_ in xs:
            xi = ha * xi_ / (G * C_G)
            tau = ha * (t - xi_ / u_i) / cs
            eg.append(schumann_J(xi, tau) if tau > 0 else 0.0)
            es.append(1.0 - schumann_J(tau, xi) if tau > 0 else 0.0)
        eg, es = np.array(eg), np.array(es)
        worst_g = max(worst_g, np.abs(np.interp(xs, x, tg) - eg).max())
        worst_s = max(worst_s, np.abs(np.interp(xs, x, ts) - es).max())
        lab = k == 0
        ax.plot(x, tg, color=SERIES, lw=1.6, label="gas, celFoam" if lab else None)
        ax.plot(x, ts, color=SERIES2, lw=1.6, label="grains, celFoam" if lab else None)
        ax.plot(xs, eg, "o", ms=3, mfc="none", mec=INK, mew=0.7,
                label="Schumann (1929)" if lab else None)
        ax.plot(xs, es, "o", ms=3, mfc="none", mec=INK, mew=0.7)
        i = np.argmin(np.abs(tg - 0.5))
        ax.annotate(f"{t:g} s", (x[i], 0.5), xytext=(4, 6), textcoords="offset points",
                    fontsize=8, color=MUTED)
    ax.set_xlabel("distance from the inlet x [m]")
    ax.set_ylabel(r"$(T - T_0)/\Delta T$")
    ax.set_xlim(0, L)
    ax.set_ylim(-0.02, 1.05)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title(f"Hot air into a cold packed bed: h = {h:.1f} W/(m$^2$ K) (Gunn), Re = {re:.0f}",
                 fontsize=10, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(f"{HERE}/profiles.png")
    msg = (f"Gunn h {h:.2f} W/(m^2 K) at Re {re:.1f}, Pr {pr:.3f}\n"
           f"largest difference from Schumann: gas {worst_g * 100:.2f} %, grains "
           f"{worst_s * 100:.2f} % of the inlet step\n")
    open(f"{HERE}/results.txt", "w").write(msg)
    print(msg, end="")


def re_num(d):
    return re.fullmatch(r"[0-9.e+-]+", d) is not None


if __name__ == "__main__":
    main()
