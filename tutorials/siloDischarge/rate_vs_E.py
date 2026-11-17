"""Discharge rate against stiffness E, WP 3.6d (decision 0045).

The numbers are from tests/manual/test_discharge_stiffness.sh (celVoidCheck,
rate over 97 to 88 percent of the mass at opening) on 2026-11-17. Writes
rate_vs_E.png. Rerun the sweep to refresh them.
"""
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

E = [0.1, 1.0, 10.0]
RUNS = [  # label, E values, rates (kg/s per m depth), colour, marker, style
    ("mu(I) + flow impedance, 10 mm (this package)", E, [65.5, 53.8, 67.7], "#2a78d6", "o", "-"),
    ("mu(I) + flow impedance, 5 mm (this package)", E[:2], [77.6, 61.2], "#1baf7a", "s", "-"),
    ("flow impedance before the fix, 10 mm", E, [62.3, 49.8, 12.4], "#eb6834", "^", "--"),
]
RHO, G, D = 1800.0, 9.81, 0.06
scale = RHO * math.sqrt(G) * D ** 1.5            # rate for C = 1

INK, MUTED, SURFACE = "#1f1f1e", "#6b6b68", "#fcfcfb"
plt.rcParams.update({"font.size": 9, "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": MUTED,
                     "figure.facecolor": SURFACE, "axes.facecolor": SURFACE})
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.axhspan(1.0 * scale, 1.9 * scale, color="#e4e4e1", lw=0)
ax.axhline(1.45 * scale, color=MUTED, lw=1, ls=":")
ax.text(0.085, 1.45 * scale + 2, "Beverloo, continuum mu(I) refs: C = 1.45 (band 1.0 to 1.9)",
        color=MUTED, fontsize=8)
for label, e, r, c, mk, ls in RUNS:
    ax.plot(e, r, color=c, marker=mk, ms=7, lw=2, ls=ls, label=label,
            markeredgecolor=SURFACE, markeredgewidth=1.5)
ax.set_xscale("log")
ax.set_xticks(E)
ax.set_xticklabels(["0.1", "1", "10"])
ax.set_xlim(0.07, 14)
ax.set_ylim(0, 2.0 * scale)
ax.set_xlabel("Young's modulus E (MPa)")
ax.set_ylabel("steady discharge rate (kg/s per m depth)")
sec = ax.secondary_yaxis("right", functions=(lambda w: w / scale, lambda c: c * scale))
sec.set_ylabel("Beverloo C")
ax.grid(color="#e4e4e1", lw=0.6, axis="y")
ax.legend(frameon=False, loc="lower left", fontsize=8)
ax.set_title("T0 tank, outlet 0.06 m: rate against stiffness", fontsize=10)
for s in ("top",):
    ax.spines[s].set_visible(False)
fig.savefig("rate_vs_E.png", dpi=150, bbox_inches="tight")
