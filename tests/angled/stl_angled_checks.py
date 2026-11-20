"""WP G.5 done-when checks: STL (cut-cell) walls with a free surface on
angled cells (decision 0107). About 2 minutes.

A 0.5 x 1.0 m box of 25 mm rows, meshed as two blocks split by an inclined
line so every column leans (tests/angled/free_surface_checks.py); the same
box as one Cartesian block is the reference. T0 Mohr-Coulomb
(void_checks.MC), gravity ramp 0.1 s, damping 100 1/s, 0.5 s.

[1] Lake in an STL tank: vertical STL walls at x = 0.06 and 0.44 m (they cut
    the leaning cells at all fractions), filled to 0.6 m: mass to 1e-10, no
    material beyond the walls (<= 1e-9 of the total), surface flat (top
    row of material at one height), s_yy hydrostatic within 5 % of rho g H
    at most and 0.5 % on average below 0.5 m, as on the Cartesian mesh.
    "Beyond" counts cells whose centre lies more than 1.5 cells past a wall
    face (cut cells at the wall hold material legitimately).
[2] Lake against an inclined STL wall: a slab turned 30 deg from vertical
    forms the right side (a hopper wall), filled to 0.6 m: mass to 1e-10,
    no material beyond the wall (<= 1e-9), and the mean wall-normal stress
    of the cells along the wall within 5 % of the Cartesian run's.

    python3 stl_angled_checks.py WORKDIR
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("void", "lagrange", "stlwalls", "silo", "angled"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_wall_checks as sw  # noqa: E402
import free_surface_checks as fs  # noqa: E402

W, H, NX, NY, FILL = 0.5, 1.0, 20, 40, 0.6
CUT = 0.0375    # a cell centre this far past a wall face is wholly outside
fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def add_bodies(p, bodies):
    body = ""
    for name, (lo, hi, rot) in bodies.items():
        cs.write(f"{p}/constant/triSurface/{name}.stl", sw.box_stl(lo, hi, rot))
        body += f"    {name}\n    {{\n        surface {name}.stl;\n        friction 0;\n    }}\n"
    cs.write(f"{p}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def lake(work, name, leaning, bodies):
    p = vc.make_case(f"{work}/{name}", NX, NY, vc.MC,
                     dict(fill=[-1, -1, 1, FILL], U0=[0, 0, 0], tEnd=0.5, maxCo=0.4))
    X0, X1 = (0.15, 0.35) if leaning else (W / 2, W / 2)
    fs.split_mesh(p, W, H, NX, NY, X0, X1, W / NX)
    open(f"{p}/system/controlDict", "a").write("gravityRamp 0.1;\ndamping 100;\n")
    add_bodies(p, bodies)
    vc.run(p)
    h = vc.history(p)
    c = np.loadtxt(f"{p}/cells.dat")
    return h, c


def beyond(c, inside):
    """Material outside the tank: alpha x open fraction summed (cells)."""
    m = ~inside(c[:, 0], c[:, 1])
    return float((c[m, 2] * c[m, 10]).sum())


def check_tank(work):
    print("[1] Lake in an STL tank on leaning cells", flush=True)
    bodies = {"left": ((-1.0, -1.0, -1.0), (0.06, 2.0, 1.0), 0.0),
              "right": ((0.44, -1.0, -1.0), (2.0, 2.0, 1.0), 0.0)}
    ref = None
    for name, leaning in (("tank_cartesian", False), ("tank_leaning", True)):
        h, c = lake(work, name, leaning, bodies)
        x, y, a, syy = c[:, 0], c[:, 1], c[:, 2], c[:, 7]
        dm = abs(h[-1, vc.MASS] / h[0, vc.MASS] - 1)
        out = beyond(c, lambda xx, yy: (xx > 0.06 - CUT) & (xx < 0.44 + CUT)) \
            / (a * c[:, 10]).sum()
        # flat: nothing above the top row that is at least half full
        tops = [y[(a > 0.5) & (c[:, 10] > 0)].max(), y[(a > 1e-6) & (c[:, 10] > 0)].max()]
        flat = tops[1] - tops[0] < 1e-3
        band = (a > 0.99) & (y > 0.1) & (y < 0.5) & (c[:, 10] > 0.99)
        e = np.abs(syy[band] + vc.RHO * vc.G_ACC * (FILL - y[band])) / (vc.RHO * vc.G_ACC * FILL)
        report(f"{name}", dm <= 1e-10 and out <= 1e-9 and flat
               and e.max() <= 0.05 and e.mean() <= 0.005,
               f"mass {dm:.1e}, beyond the walls {out:.1e}, top half-full row {tops[0]:.4f} m, "
               f"highest material {tops[1]:.4f} m, s_yy error largest "
               f"{e.max() * 100:.2f} %, mean {e.mean() * 100:.2f} %")


def check_inclined(work):
    print("[2] Lake against an inclined STL wall (30 deg from vertical)", flush=True)
    th = math.radians(30.0)
    # slab whose left face is the line through (0.30, 0) leaning right by
    # 30 deg: a box turned by -30 deg about the origin (sw.box_stl), its left
    # face through the pre-image of (0.30, 0)
    x0 = 0.30 * math.cos(th)
    y0 = 0.30 * math.sin(th)
    lo, hi, rot = (x0, y0 - 2.0, -1.0), (x0 + 1.0, y0 + 2.0, 1.0), -30.0
    bodies = {"left": ((-1.0, -1.0, -1.0), (0.06, 2.0, 1.0), 0.0),
              "slope": (lo, hi, rot)}
    # inside: right of the left wall and left of the slope face
    nrm = np.array([math.cos(th), -math.sin(th)])        # out of the material

    def inside(xx, yy):
        return (xx > 0.06 - CUT) & ((xx - 0.30) * nrm[0] + yy * nrm[1] < CUT)
    res = {}
    for name, leaning in (("slope_cartesian", False), ("slope_leaning", True)):
        h, c = lake(work, name, leaning, bodies)
        x, y, a = c[:, 0], c[:, 1], c[:, 2]
        dm = abs(h[-1, vc.MASS] / h[0, vc.MASS] - 1)
        out = beyond(c, inside) / (a * c[:, 10]).sum()
        d = -((x - 0.30) * nrm[0] + y * nrm[1])               # depth from the wall face
        sel = (d > 0) & (d < 0.04) & (a > 0.99) & (y > 0.1) & (y < 0.5)
        sn = c[:, 6] * nrm[0] ** 2 + 2 * c[:, 8] * nrm[0] * nrm[1] + c[:, 7] * nrm[1] ** 2
        res[name] = sn[sel].mean()
        report(f"{name}: mass and no leak", dm <= 1e-10 and out <= 1e-9,
               f"mass {dm:.1e}, beyond the wall {out:.1e}")
    r = res["slope_leaning"] / res["slope_cartesian"]
    report("wall-normal stress against Cartesian", abs(r - 1) <= 0.05,
           f"leaning {res['slope_leaning']:.0f} Pa, Cartesian {res['slope_cartesian']:.0f} Pa, "
           f"ratio {r:.3f}")


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    check_tank(work)
    check_inclined(work)
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
