#!/usr/bin/env python3
"""WP 3.3p: figures and animation of the prototype silo discharge.

Reads .npz files written by ``prototype/discharge2d.py --run --save`` and
writes PNG figures and a GIF animation into an output folder. Needs
matplotlib (not used by the tests).

    python3 tools/plot_discharge.py --frames run.npz --out docs/results/silo-prototype
    python3 tools/plot_discharge.py --mass a.npz b.npz ... --labels ... --out ...
"""

import argparse
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import animation  # noqa: E402

RHO = 1800.0
G = 9.81
C_REF = 1.45


def load(p):
    return dict(np.load(p))


def field_axes(ax, r, title):
    ax.set_xlim(0, float(r["lx"]))
    ax.set_ylim(-0.02, float(r["ly"]))
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=9)
    ax.set_xticks([0, 0.2, 0.4])
    ax.tick_params(labelsize=7)
    # tank outline and outlet
    lx = float(r["lx"])
    cell = float(r["cell"])
    out = r["outlet"]
    x0 = np.argmax(out) * cell
    x1 = x0 + out.sum() * cell
    ax.plot([0, 0, x0], [r["ly"], 0, 0], "k-", lw=1.5)
    ax.plot([x1, lx, lx], [0, 0, r["ly"]], "k-", lw=1.5)


def frame_image(ax, r, k, kind):
    al = r["f_al"][k]
    ext = (0, float(r["lx"]), 0, float(r["ly"]))
    if kind == "speed":
        v = np.hypot(r["f_vx"][k], r["f_vy"][k])
        q = np.ma.masked_where(al < 0.5, v)
        im = ax.imshow(q, origin="lower", extent=ext, cmap="viridis", vmin=0.0,
                       vmax=float(np.sqrt(2 * G * r["D"] / 2)) * 1.5)
    elif kind == "epq":
        q = np.ma.masked_where(al < 0.5, np.log10(np.maximum(r["f_epq"][k], 1e-4)))
        im = ax.imshow(q, origin="lower", extent=ext, cmap="magma", vmin=-4, vmax=0.5)
    else:
        im = ax.imshow(al, origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=1)
    return im


def fig_fields(r, out, times=(0.0, 0.25, 0.5, 1.0, 1.5)):
    ft = r["f_t"]
    ks = [int(np.argmin(np.abs(ft - t))) for t in times if t <= ft[-1] + 1e-9]
    for kind, label, fname in (("speed", "speed (m/s)", "speed"),
                               ("epq", "log10 equivalent plastic strain", "plastic_strain")):
        fig, axs = plt.subplots(1, len(ks), figsize=(2.2 * len(ks), 3.6), constrained_layout=True)
        for ax, k in zip(np.atleast_1d(axs), ks):
            im = frame_image(ax, r, k, kind)
            field_axes(ax, r, f"t = {ft[k]:.2f} s after opening")
        fig.colorbar(im, ax=axs, shrink=0.8, label=label)
        p = os.path.join(out, f"fields_{fname}.png")
        fig.savefig(p, dpi=130)
        plt.close(fig)
        print("wrote", p)


def fig_surface(r, out, times=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)):
    ft = r["f_t"]
    fig, ax = plt.subplots(figsize=(5, 3.5), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    ts = [t for t in times if t <= ft[-1] + 1e-9]
    for i, t in enumerate(ts):
        k = int(np.argmin(np.abs(ft - t)))
        h = r["f_al"][k].sum(axis=0) * float(r["cell"])
        ax.step(r["xc"], h, where="mid", color=cmap(i / max(1, len(ts) - 1)),
                label=f"{ft[k]:.2f} s")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("material height (m)")
    ax.set_title("Free surface after opening")
    ax.legend(fontsize=7, ncol=2)
    p = os.path.join(out, "free_surface.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print("wrote", p)


def gif(r, out, fps=12):
    ft = r["f_t"]
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 5.2), constrained_layout=True)
    ims = []
    for ax, kind, label in ((axs[0], "speed", "speed (m/s)"),
                            (axs[1], "epq", "log10 plastic strain")):
        im = frame_image(ax, r, 0, kind)
        field_axes(ax, r, label)
        fig.colorbar(im, ax=ax, shrink=0.7, orientation="horizontal")
        ims.append(im)
    txt = fig.suptitle("")

    def update(k):
        for im, kind in zip(ims, ("speed", "epq")):
            al = r["f_al"][k]
            if kind == "speed":
                q = np.ma.masked_where(al < 0.5, np.hypot(r["f_vx"][k], r["f_vy"][k]))
            else:
                q = np.ma.masked_where(al < 0.5, np.log10(np.maximum(r["f_epq"][k], 1e-4)))
            im.set_data(q)
        txt.set_text(f"T0 silo discharge, D = {1000 * float(r['D']):.0f} mm: "
                     f"t = {ft[k]:.2f} s after opening")
        return ims + [txt]

    anim = animation.FuncAnimation(fig, update, frames=len(ft), blit=False)
    p = os.path.join(out, "discharge.gif")
    anim.save(p, writer=animation.PillowWriter(fps=fps), dpi=80)
    plt.close(fig)
    print("wrote", p)


def fig_mass(runs, labels, out, name="mass"):
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    for r, lab in zip(runs, labels):
        h = r["hist"]
        axs[0].plot(h[:, 0], h[:, 1], label=lab)
        t = h[:, 0]
        m = h[:, 1]
        # rate over 0.05 s windows
        dt = 0.05
        tb = np.arange(0.0, t[-1] - dt, dt)
        rate = [-(np.interp(a + dt, t, m) - np.interp(a, t, m)) / dt for a in tb]
        line, = axs[1].plot(tb + 0.5 * dt, rate, label=lab)
        D = float(r["D"])
        axs[1].axhline(C_REF * RHO * math.sqrt(G) * D ** 1.5, color=line.get_color(), ls=":",
                       lw=1)
    axs[0].set_xlabel("time after opening (s)")
    axs[0].set_ylabel("mass in tank (kg per m depth)")
    axs[1].set_xlabel("time after opening (s)")
    axs[1].set_ylabel("discharge rate (kg/s per m depth)")
    axs[1].set_title("dotted: Beverloo, C = 1.45", fontsize=9)
    axs[0].legend(fontsize=7)
    p = os.path.join(out, f"{name}.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print("wrote", p)


def fig_law(runs, out):
    Ds = np.array([float(r["D"]) for r in runs])
    rates = np.array([float(r["rate"]) for r in runs])
    fig, ax = plt.subplots(figsize=(4.5, 3.6), constrained_layout=True)
    ax.loglog(Ds, rates, "o", label="prototype (steady rate)")
    dd = np.linspace(0.03, 0.1, 50)
    for C, ls in ((1.0, ":"), (C_REF, "-"), (1.9, ":")):
        ax.loglog(dd, C * RHO * math.sqrt(G) * dd ** 1.5, "k" + ls, lw=1,
                  label=f"C = {C}" if C == C_REF else None)
    ax.set_xlabel("outlet width D (m)")
    ax.set_ylabel("rate (kg/s per m depth)")
    n = np.polyfit(np.log(Ds), np.log(rates), 1)[0]
    ax.set_title(f"Beverloo 2D: fitted exponent {n:.2f} (1.5)", fontsize=9)
    ax.legend(fontsize=7)
    p = os.path.join(out, "rate_law.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print("wrote", p)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--frames", help="run with frames: field figures, surface, GIF")
    p.add_argument("--mass", nargs="*", default=[])
    p.add_argument("--labels", nargs="*", default=None)
    p.add_argument("--mass-name", default="mass")
    p.add_argument("--law", nargs="*", default=[])
    p.add_argument("--out", required=True)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.frames:
        r = load(a.frames)
        fig_fields(r, a.out)
        fig_surface(r, a.out)
        gif(r, a.out)
    if a.mass:
        runs = [load(f) for f in a.mass]
        fig_mass(runs, a.labels or a.mass, a.out, a.mass_name)
    if a.law:
        fig_law([load(f) for f in a.law], a.out)


if __name__ == "__main__":
    main()
