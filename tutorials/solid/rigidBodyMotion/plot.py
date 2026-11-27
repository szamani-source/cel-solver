"""Figure for the rigidBodyMotion tutorial (run after ./Allrun; numpy,
matplotlib): bodies.png, the three bodies seen from the front (x-z) and from
above (x-y) at t = 0 and at the end, from postProcessing/rigidBodies."""
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "postProcessing", "rigidBodies")
COL = {"hopperWall": "#6b6b68", "pusher": "#2a78d6", "spinner": "#b8322a"}


def read_stl(path):
    txt = open(path).read()
    v = np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", txt), float)
    return v.reshape(-1, 3, 3)


def main():
    dirs = sorted((float(os.path.basename(d)), d) for d in glob.glob(f"{PP}/*")
                  if os.path.isdir(d) and re.match(r"^[0-9.e+-]+$", os.path.basename(d)))
    first, last = dirs[0], dirs[-1]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for ax, (i, j), lab in zip(axs, ((0, 2), (0, 1)), (("x", "z"), ("x", "y"))):
        for (t, d), alpha, ls in ((first, 0.15, "--"), (last, 0.45, "-")):
            for f in sorted(glob.glob(f"{d}/*.stl")):
                name = os.path.basename(f)[:-4]
                tri = read_stl(f)
                pc = PolyCollection(tri[:, :, [i, j]], facecolor=COL.get(name, "#999"),
                                    edgecolor=COL.get(name, "#999"), alpha=alpha,
                                    linewidths=0.3, linestyles=ls)
                ax.add_collection(pc)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        ax.set_xlabel(f"{lab[0]} (m)")
        ax.set_ylabel(f"{lab[1]} (m)")
        ax.set_title(f"{'front' if j == 2 else 'top'} view: t = {first[0]:g} s pale, "
                     f"t = {last[0]:g} s dark", fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.5) for c in COL.values()]
    fig.legend(handles, COL.keys(), loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(f"{HERE}/bodies.png", dpi=150)
    print("bodies.png written")


if __name__ == "__main__":
    main()
