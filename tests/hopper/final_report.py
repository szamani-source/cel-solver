"""Final Couto silo comparison on the 4 mm mesh (decision 0117).

Reads the MC and DP cases of tests/hopper/run_final_4mm.sh and compares, as
the papers report them (decision 0117, item 4):

- at rest (t = 0.5 s, before opening): GA0-GA7 against Couto's static range;
- in discharge: our largest value over the writes 0.25-1.0 s after opening
  (the Couto table gives the maximum during discharge) against Couto's
  discharge range; the mean is printed too;
- the wall-pressure profile 0.5 s after opening (t = 1.0 s) against Zheng
  and Yu's Abaqus CEL (Fig. 11);
- the discharge rate (mass slope 0.25-1.0 s after opening) against Beverloo
  (110-131 kg/s, docs/spec/conical-hopper.md).

Cylinder cells: normal pressure = -sigma_xx; hopper (GA0, profile below the
transition): n.sigma.n with the cone normal (hopper_case.wall_pressures).

    python3 tests/hopper/final_report.py OUTDIR MC_CASE [DP_CASE]

Writes OUTDIR/couto-final-wall-pressure.png, OUTDIR/couto-final-mass.png and
OUTDIR/couto-final-table.md.
"""
import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hopper_pressure_checks as hp  # noqa: E402
import hopper_case as hc  # noqa: E402

CELL = 0.012                       # wall band of the aligned_buffered cases
T_OPEN = hp.T_OPEN
DIS = [w for w in hp.WRITES if hp.AVG[0] <= w <= hp.AVG[1]]
NAMES = ["GA0", "GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7"]
BEV = (110.0, 131.0)
COLS = {"MC": "#c0392b", "DP": "#1f6fb4"}
THETA = math.radians(hc.THETA)
F_WEDGE = math.pi / (math.sin(THETA) * math.cos(THETA))


def has(p, name):
    return os.path.exists(f"{p}/{name}.dat") or os.path.exists(f"{p}/{name}.dat.proc0")


def profile(p, c, name):
    np.savetxt(f"{p}/cells.dat", hp.read_cells(p, name))
    cyl, hop = hc.wall_pressures(p, c, CELL)
    zs = np.unique(np.round(hop[:, 0], 3))
    h = np.array([[z, hop[np.abs(hop[:, 0] - z) < 1e-3, 1].mean()] for z in zs])
    return np.vstack([h, cyl])


def analyse(label, p):
    c = dict(hc.COUTO)
    hp.CELL = CELL
    # columns step t mass ... Umax (9); a cut-off last line is skipped
    h = np.genfromtxt(f"{p}/history.dat", usecols=(0, 1, 2, 9), invalid_raise=False)
    h = h[~np.isnan(h).any(1)]
    t, m = h[:, 1], h[:, 2] * F_WEDGE
    i = (t >= hp.AVG[0]) & (t <= hp.AVG[1])
    rate = -np.polyfit(t[i], m[i], 1)[0] if i.sum() > 10 else float("nan")
    rest = hp.pressures(p, c, f"cells_{T_OPEN:g}")
    ws = [w for w in DIS if has(p, f"cells_{w:g}")]
    dis = [hp.pressures(p, c, f"cells_{w:g}") for w in ws]
    snap = f"cells_{T_OPEN + 0.5:g}"
    return dict(label=label, p=p, c=c, t=t, m=m, rate=rate, rest=rest, ws=ws,
                dmax={n: max(d[n] for d in dis) for n in NAMES} if dis else {},
                dmean={n: np.mean([d[n] for d in dis]) for n in NAMES} if dis else {},
                prof_rest=profile(p, c, f"cells_{T_OPEN:g}"),
                prof_snap=profile(p, c, snap) if has(p, snap) else None,
                umax=h[t > T_OPEN, 3].max() if (t > T_OPEN).any() else float("nan"))


def verdict(v, lo, hi, tol=0.15):
    if lo <= v <= hi:
        return "inside"
    return f"{(v / hi - 1) * 100:+.0f} %" if v > hi else f"{(v / lo - 1) * 100:+.0f} %"


def main():
    out = sys.argv[1]
    runs = [analyse(l, p) for l, p in zip(["MC", "DP"], sys.argv[2:]) if os.path.isdir(p)]
    rows = hp.table()
    os.makedirs(out, exist_ok=True)
    fem = np.loadtxt(os.path.join(hp.ROOT, "docs", "spec", "hopper-data",
                                  "zheng2015_fig11_fem.csv"), delimiter=",", comments=("#", "z_m"))
    zfem = dict(zip(np.round(fem[:, 0], 2), fem[:, 1]))

    # table
    L = ["| Quantity | Reference | " + " | ".join(r["label"] for r in runs) + " |",
         "|---|---|" + "---|" * len(runs)]
    L.append(f"| Discharge rate, kg/s | Beverloo {BEV[0]:.0f}-{BEV[1]:.0f} | "
             + " | ".join(f"{r['rate']:.1f} ({verdict(r['rate'], *BEV)})" for r in runs) + " |")
    for n in NAMES:
        lo, hi = hp.exp_range(rows, "static", n)
        L.append(f"| {n} at rest, kPa | Couto {lo:.2f}-{hi:.2f} | "
                 + " | ".join(f"{r['rest'][n]:.2f} ({verdict(r['rest'][n], lo, hi)})" for r in runs) + " |")
    for n in NAMES:
        lo, hi = hp.exp_range(rows, "discharge", n)
        L.append(f"| {n} discharge max (mean), kPa | Couto max {lo:.2f}-{hi:.2f} | "
                 + " | ".join(f"{r['dmax'].get(n, float('nan')):.2f} ({r['dmean'].get(n, float('nan')):.2f}); "
                              f"{verdict(r['dmax'].get(n, float('nan')), lo, hi)}" for r in runs) + " |")
    for r in runs:
        if r["prof_snap"] is not None:
            pr = r["prof_snap"]
            pk = pr[pr[:, 0] < 0.6]
            k = np.argmax(pk[:, 1])
            L.append(f"| Peak at 0.5 s after opening, {r['label']} | Abaqus 27.5 kPa at z 0.46 m | "
                     f"{pk[k, 1]:.1f} kPa at z {pk[k, 0]:.2f} m |" + " |" * (len(runs) - 1))
    L.append("| Largest speed after opening, m/s | | " + " | ".join(f"{r['umax']:.2f}" for r in runs) + " |")
    L.append("")
    L.append("Discharge writes used: " + "; ".join(f"{r['label']} {len(r['ws'])}" for r in runs))
    txt = "\n".join(L)
    open(f"{out}/couto-final-table.md", "w").write(txt + "\n")
    print(txt)

    # wall pressure figure
    cells = {n: z for n, z in hp.GA.items()}
    cells["GA0"] = 0.41
    fig, axs = plt.subplots(1, 3, figsize=(16, 7.5), dpi=150, sharey=True)
    for ax, phase, title in ((axs[0], "static", "At rest (t = 0.5 s, before opening)"),
                             (axs[1], "discharge", "Discharge: max over 0.25-1.0 s after opening"),
                             (axs[2], None, "0.5 s after opening (snapshot)")):
        if phase:
            for k, (n, z) in enumerate(cells.items()):
                lo, hi = hp.exp_range(rows, phase, n)
                ax.plot([lo, hi], [z, z], color="k", lw=5, alpha=0.35, solid_capstyle="butt",
                        label="Couto et al. 2013, 7 tests" + (" (max)" if phase == "discharge" else "")
                        if k == 0 else None)
                ax.text(hi + 0.4, z, n, va="center", fontsize=8, color="#444")
        for r in runs:
            col = COLS[r["label"]]
            if phase == "static":
                pr = r["prof_rest"]
                ax.plot(pr[:, 1], pr[:, 0], color=col, lw=1.6, label=f"celFoam {r['label']}")
            elif phase == "discharge":
                if not r["dmax"]:
                    continue
                zz = [cells[n] for n in NAMES]
                ax.plot([r["dmax"][n] for n in NAMES], zz, "o", color=col, ms=6,
                        label=f"celFoam {r['label']}: max at the load cells")
                ax.plot([r["dmean"][n] for n in NAMES], zz, "x", color=col, ms=6, alpha=0.7,
                        label=f"celFoam {r['label']}: mean")
            elif r["prof_snap"] is not None:
                pr = r["prof_snap"]
                ax.plot(pr[:, 1], pr[:, 0], color=col, lw=1.6, label=f"celFoam {r['label']}")
        if phase is None:
            ax.plot(fem[:, 1], fem[:, 0], "--", color="#2a9d55", lw=1.6,
                    label="Zheng & Yu 2015, Abaqus CEL")
        ax.axhline(0.4764, color="#999", lw=0.6, ls=":")
        ax.text(0.3, 0.495, "transition", fontsize=7, color="#777")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("normal wall pressure (kPa)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, loc="upper right")
    axs[0].set_ylabel("height above the outlet z (m)")
    axs[0].set_ylim(0, 2.5)
    axs[0].set_xlim(0, 16)
    axs[1].set_xlim(0, 35)
    axs[2].set_xlim(0, 35)
    fig.suptitle("Couto silo on the 4 mm mesh (20,398 cells): normal wall pressure "
                 "(cylinder horizontal, hopper normal to the cone)", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{out}/couto-final-wall-pressure.png")

    # mass figure
    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)
    for r in runs:
        ax.plot(r["t"] - T_OPEN, r["m"], color=COLS[r["label"]], lw=1.6,
                label=f"celFoam {r['label']}: {r['rate']:.1f} kg/s")
    m0 = runs[0]["m"][np.argmin(abs(runs[0]["t"] - T_OPEN))]
    tt = np.linspace(0, 1.0, 11)
    ax.fill_between(tt, m0 - BEV[0] * tt, m0 - BEV[1] * tt, color="#999", alpha=0.3,
                    label=f"Beverloo {BEV[0]:.0f}-{BEV[1]:.0f} kg/s")
    ax.axvspan(0.25, 1.0, color="#f3e6c4", alpha=0.4, lw=0, label="rate fitted here")
    ax.set_xlim(-0.1, 1.0)
    ax.set_xlabel("time after opening (s)")
    ax.set_ylabel("mass in the silo (kg)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    ax.set_title("Couto silo discharge, 4 mm mesh")
    fig.tight_layout()
    fig.savefig(f"{out}/couto-final-mass.png")


if __name__ == "__main__":
    main()
