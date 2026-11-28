"""Figure for the smokeBox tutorial (run after ./Allrun; numpy, matplotlib,
OpenFOAM sourced for the cell centres): profiles.png, the x velocity and the
density along x at t = 0 and at the last write, averaged over y and z."""
import os
import re
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, MUTED, SERIES = "#1f1f1e", "#6b6b68", "#2a78d6"


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
    for d in os.listdir(HERE):
        try:
            out.append((float(d), d))
        except ValueError:
            pass
    return sorted(out)


def main():
    if not os.path.exists(f"{HERE}/0/C"):
        subprocess.run(["postProcess", "-func", "writeCellCentres", "-time", "0"],
                       cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    C = read_internal(f"{HERE}/0/C", 10 ** 9)
    n = len(C)
    xs = np.unique(np.round(C[:, 0], 9))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for (t, d), col, ls in zip([times()[0], times()[-1]], (MUTED, SERIES), ("--", "-")):
        U = read_internal(f"{HERE}/{d}/U", n)
        rho = read_internal(f"{HERE}/{d}/rho", n)
        ux = [U[np.isclose(C[:, 0], x), 0].mean() for x in xs]
        r = [rho[np.isclose(C[:, 0], x)].mean() for x in xs]
        axs[0].plot(xs, ux, ls, color=col, marker="o", ms=3, label=f"t = {t:g} s")
        axs[1].plot(xs, r, ls, color=col, marker="o", ms=3, label=f"t = {t:g} s")
    axs[0].set_ylabel("U$_x$ (m/s)")
    axs[1].set_ylabel("density (kg/m$^3$)")
    for ax in axs:
        ax.set_xlabel("x (m)")
        ax.legend(frameon=False, fontsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{HERE}/profiles.png", dpi=150)
    print("profiles.png written")


if __name__ == "__main__":
    main()
