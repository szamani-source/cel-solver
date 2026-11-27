"""Plots for the fluidisedDamBreak tutorial (WP F.9; case of WP F.8,
decision 0133).

Run after ./Allrun (needs numpy and matplotlib):

    python3 plot.py

Writes front.png (front of the deposit against time, fluidised and dry,
with Roche et al. 2008), snapshots.png (material fraction at four times)
and front.txt (run-out and front speed).

The front is the end of the continuous deposit from the gate: the first
column holding less than half a cell of material (grains sprayed ahead of
the flow do not count). Front speed: least squares over Roche's
constant-velocity phase, t/sqrt(h0/g) from 1.5 to 4.
"""

import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
G, X0, H0, LX, LY = 9.81, 0.2, 0.4, 3.0, 0.45
TS = math.sqrt(H0 / G)                    # sqrt(h0/g) [s]
U_ROCHE = math.sqrt(2 * G * H0)           # Roche's front speed [m/s]
RUNOUT = (5.5, 6.0)                       # Roche's run-out from the gate [h0]
PHASE = (1.5, 4.0)
SNAPS = [0.0, 2.0, 4.0, 9.0]              # t/sqrt(h0/g)
INK, MUTED, SURFACE, GAS, DRY = "#1f1f1e", "#6b6b68", "#fcfcfb", "#c2410c", "#2a78d6"


def read_scalar(path):
    txt = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?(\d+)\s*\n?\(", txt)
    if not m:
        n = None
        u = re.search(r"internalField\s+uniform\s+(\S+);", txt)
        return float(u.group(1)), n
    n = int(m.group(1))
    return np.array(txt[m.end():].split(")")[0].split(), float), n


def times(case):
    return sorted((float(d), d) for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d))


def history(case, cell, nx, ny):
    rows = []
    for t, d in times(case):
        a, _ = read_scalar(f"{case}/{d}/alpha")
        hc = np.broadcast_to(a, (nx * ny,)).reshape(ny, nx).sum(axis=0) * cell
        k = int(round(X0 / cell))
        while k < nx and hc[k] >= 0.5 * cell:
            k += 1
        rows.append((t, k * cell))
    return np.array(rows)


def speed(h):
    tau = h[:, 0] / TS
    s = (tau >= PHASE[0]) & (tau <= PHASE[1])
    return float(np.polyfit(h[s, 0], h[s, 1], 1)[0])


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def main():
    cell = float(open(f"{HERE}/cellSize").read()) / 1000
    nx, ny = int(round(LX / cell)), int(round(LY / cell))
    cases = [("fluidised", HERE, GAS)]
    if os.path.isdir(f"{HERE}/dry") and times(f"{HERE}/dry"):
        cases.append(("dry", f"{HERE}/dry", DRY))
    hist = {name: history(path, cell, nx, ny) for name, path, _ in cases}

    lines = [f"# {cell * 1e3:g} mm cells; run-out from the gate [h0], front speed [m/s]"]
    fig, ax = plt.subplots(figsize=(6.4, 3.8), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    style(ax)
    tau = np.linspace(0, 9, 50)
    ax.fill_between([0, 9], RUNOUT[0], RUNOUT[1], color=MUTED, alpha=0.12, lw=0)
    ax.plot(tau, np.minimum(tau * U_ROCHE * TS / H0, RUNOUT[1]), color=MUTED, lw=1, ls="--",
            label="Roche et al. (2008): sqrt(2 g h0), run-out 5.5-6 h0")
    for name, _, colour in cases:
        h = hist[name]
        ro, u = (h[-1, 1] - X0) / H0, speed(h)
        lines.append(f"{name}: run-out {ro:.2f}, front speed {u:.3f} ({100 * (u / U_ROCHE - 1):+.1f} % of sqrt(2 g h0))")
        ax.plot(h[:, 0] / TS, (h[:, 1] - X0) / H0, color=colour, lw=2,
                label=f"{name}: run-out {ro:.2f} h0, front {u:.2f} m/s")
    ax.set_xlabel("t / sqrt(h0/g)", color=INK)
    ax.set_ylabel("front from the gate [h0]", color=INK)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right")
    fig.tight_layout()
    fig.savefig(f"{HERE}/front.png", facecolor=SURFACE)
    open(f"{HERE}/front.txt", "w").write("\n".join(lines) + "\n")

    fig, axs = plt.subplots(len(SNAPS), len(cases), figsize=(5.2 * len(cases), 1.3 * len(SNAPS) + 0.6),
                            dpi=150, squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for c, (name, path, _) in enumerate(cases):
        ts = times(path)
        for r, s in enumerate(SNAPS):
            t, d = min(ts, key=lambda td: abs(td[0] - s * TS))
            a, _ = read_scalar(f"{path}/{d}/alpha")
            a = np.broadcast_to(a, (nx * ny,)).reshape(ny, nx)
            ax = axs[r, c]
            ax.imshow(a, origin="lower", cmap="Greys", vmin=0, vmax=1, extent=(0, LX, 0, LY),
                      aspect="equal", interpolation="nearest")
            ax.set_xlim(0, LX)
            ax.set_yticks([])
            ax.tick_params(colors=MUTED, labelsize=7)
            ax.text(0.99, 0.75, f"{name}, t = {t / TS:.0f} sqrt(h0/g)", transform=ax.transAxes,
                    ha="right", color=INK, fontsize=8)
            if r < len(SNAPS) - 1:
                ax.set_xticklabels([])
        axs[-1, c].set_xlabel("x [m]", color=INK)
    fig.tight_layout()
    fig.savefig(f"{HERE}/snapshots.png", facecolor=SURFACE)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
