"""WP 4.5 done-when checks: silo and hopper discharge with STL walls
(decisions 0047, 0051). Run by tests/test_stl_silo.sh after celVoidCheck is
built:

    python3 stl_silo_checks.py WORKDIR [tank] [flatstl] [hopper]

flatstl runs one STL flat tank at maxCo 0.40 for the hopper to compare with
(decision 0091); CEL_CELL=0.005 runs everything at 5 mm cells.

The T0 tank of WP 3.3/3.6d (tests/outlet/outlet_checks.py): W 0.4 m, fill
0.6 m, outlet D 0.06 m, Mohr-Coulomb c 0, phi 35 with the T0 mu(I) law and
the flow impedance (decisions 0037, 0045), walls 20 deg, 10 mm cells, gravity
ramp 0.1 s, damping 250 1/s, outlet opens at 0.3 s. E = 0.1 MPa (the rate
does not depend on E within the scatter, decision 0045; 0.1 MPa runs are
about 25 times shorter than 10 MPa).

A single 10 mm rate scatters by about 15 % with round-off (decision 0046),
so each tank is run three times with maxCo 0.40, 0.38, 0.36 (a different
round-off history each) and the mean rates are compared.

1. Flat tank with STL walls: side walls and floor are STL bodies (the
   tutorials/solid/stlSilo tank: walls cut cells at a quarter and three quarters,
   floor at y = 0.0125 with a 0.06 m gap over the mesh outlet patch). No
   material leaves while the outlet is closed (mass to 1e-12), tank plus
   discharged constant to 1e-10, and the mean rate within 5 % of the mean
   rate of the mesh-wall tank.
2. Hopper (no pass band until Sia agrees one): the same tank with STL hopper
   walls at 60 deg to the horizontal (half angle 30 deg) from the gap edges
   up to the side walls, filled to the same top, runs until it is empty.
   Records the rate and the Beverloo constant against the flat tank, and
   the Rose-Tanaka hopper factor (tan 30 deg)^-0.35 = 1.21 for comparison.
   Checks only mass balance and that it empties (below 2 % of the fill).

Prints PASS/FAIL per check, exits non-zero on any failure.
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("outlet", "void", "lagrange", "silo"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import outlet_checks as oc  # noqa: E402
import cases as cs  # noqa: E402

MUI = "mu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650.0;\nflowImpedance yes;\n"
CO = [0.40, 0.38, 0.36]
TOL = 0.05
FLOOR = 0.0125                     # STL floor top
X0, NX = -0.0275, 46               # STL meshes: x = X0 .. X0 + NX cells
GAP = (0.17, 0.23)
HOPPER = 60.0                      # wall angle to the horizontal [deg]
T_HOP = 6.0                        # hopper: run this long after opening (empties in about 5 s)


def prism_stl(name, poly, z=(-0.01, 0.02)):
    """ASCII STL of the prism of a convex polygon (x, y) listed
    anticlockwise, between z[0] and z[1], outward normals."""
    out = [f"solid {name}"]

    def tri(a, b, c):
        a, b, c = (np.array(p, float) for p in (a, b, c))
        n = np.cross(b - a, c - a)
        n /= np.linalg.norm(n)
        out.append(f"facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}\nouter loop")
        for v in (a, b, c):
            out.append(f"vertex {v[0]:.12e} {v[1]:.12e} {v[2]:.12e}")
        out.append("endloop\nendfacet")
    lo = [(x, y, z[0]) for x, y in poly]
    hi = [(x, y, z[1]) for x, y in poly]
    for i in range(1, len(poly) - 1):
        tri(lo[0], lo[i + 1], lo[i])            # bottom, normal -z
        tri(hi[0], hi[i], hi[i + 1])            # top, normal +z
    for i in range(len(poly)):
        j = (i + 1) % len(poly)
        tri(lo[i], lo[j], hi[j])
        tri(lo[i], hi[j], hi[i])
    out.append(f"endsolid {name}")
    return "\n".join(out) + "\n"


def box(lo, hi):
    return [(lo[0], lo[1]), (hi[0], lo[1]), (hi[0], hi[1]), (lo[0], hi[1])]


def stl_mesh(path, ny):
    """Box mesh x X0 .. X0 + NX cells, y 0 .. ny cells, three blocks so that
    the cells under the floor gap (x 0.1625 to 0.2325) are the outlet."""
    c = oc.CELL
    k = int(round(0.01 / c))           # cells per 10 mm (CEL_CELL 0.005: 2)
    n1, n2, nx = 19 * k, 26 * k, NX * k
    xs = [X0, X0 + n1 * c, X0 + n2 * c, X0 + nx * c]
    Hd = ny * c
    verts = ([f"({x!r} 0 0)" for x in xs] + [f"({x!r} {Hd!r} 0)" for x in xs]
             + [f"({x!r} 0 {c!r})" for x in xs] + [f"({x!r} {Hd!r} {c!r})" for x in xs])
    ncell = [n1, n2 - n1, nx - n2]
    blocks = [f"    hex ({b} {b + 1} {5 + b} {4 + b} {8 + b} {9 + b} {13 + b} {12 + b}) "
              f"({ncell[b]} {ny} 1) simpleGrading (1 1 1)" for b in range(3)]

    def ymin(b):
        return f"({b} {b + 1} {9 + b} {8 + b})"
    patches = {
        "left": "(0 8 12 4)", "right": "(3 7 15 11)",
        "bottom": ymin(0) + " " + ymin(2), "outlet": ymin(1),
        "top": " ".join(f"({4 + b} {12 + b} {13 + b} {5 + b})" for b in range(3)),
    }
    pt = [f"    {n}\n    {{\n        type patch;\n        faces ({f});\n    }}"
          for n, f in patches.items()]
    empty = " ".join(f"({b} {4 + b} {5 + b} {b + 1}) ({8 + b} {9 + b} {13 + b} {12 + b})"
                     for b in range(3))
    pt.append(f"    frontAndBack\n    {{\n        type empty;\n        faces ({empty});\n    }}")
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + "scale 1;\nvertices\n(\n    " + "\n    ".join(verts) + "\n);\nblocks\n(\n"
             + "\n".join(blocks) + "\n);\nboundary\n(\n" + "\n".join(pt) + "\n);\n")


def add_bodies(path, bodies):
    """bodies: name -> polygon; friction DELTA_W on each."""
    txt = ""
    for name, poly in bodies.items():
        cs.write(f"{path}/constant/triSurface/{name}.stl", prism_stl(name, poly))
        txt += (f"    {name}\n    {{\n        surface {name}.stl;\n"
                f"        friction {oc.DELTA_W!r};\n    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + txt + "}\n")


def set_check(path, **kw):
    fn = f"{path}/system/voidCheckDict"
    s = open(fn).read()
    for k, v in kw.items():
        s = re.sub(rf"^{k} .*;$", f"{k} {cs.foam_value(v)};", s, flags=re.M)
    open(fn, "w").write(s)


def make(work, name, kind, co, t_dis=oc.T_DIS):
    oc.E = 1.0e5
    oc.T_DIS = t_dis
    path = oc.setup(work, name)
    open(f"{path}/constant/materialProperties", "a").write(MUI)
    set_check(path, maxCo=co)
    if kind == "mesh":
        return path
    Hfill = oc.H + FLOOR
    ny = int(math.ceil(Hfill / oc.CELL - 1e-9)) + oc.HEADROOM
    stl_mesh(path, ny)
    # fill: everything open below the top (the bodies take the rest)
    set_check(path, fill=[-1.0, -1.0, 2.0, Hfill])
    top = ny * oc.CELL + 0.1
    bodies = {"leftWall": box((-0.1, -0.1), (0.0, top)),
              "rightWall": box((0.4, -0.1), (0.5, top)),
              "leftFloor": box((-0.1, -0.1), (GAP[0], FLOOR)),
              "rightFloor": box((GAP[1], -0.1), (0.5, FLOOR))}
    if kind == "hopper":
        r = GAP[0] * math.tan(math.radians(HOPPER))
        bodies["leftFloor"] = [(-0.1, -0.1), (GAP[0], -0.1), (GAP[0], FLOOR),
                               (0.0, FLOOR + r), (-0.1, FLOOR + r)]
        bodies["rightFloor"] = [(GAP[1], -0.1), (0.5, -0.1), (0.5, FLOOR + r),
                                (0.4, FLOOR + r), (GAP[1], FLOOR)]
    add_bodies(path, bodies)
    return path


def result(path):
    h = np.loadtxt(f"{path}/history.dat")
    t, m, out = h[:, oc.T], h[:, oc.MASS] / oc.CELL, h[:, oc.OUT] / oc.CELL
    closed = t <= oc.T_OPEN + 1e-12
    i_open = int(np.where(closed)[0][-1])
    m_open = m[i_open]
    rate, win = oc.steady_rate(t[i_open:] - oc.T_OPEN, m[i_open:], m_open)
    return dict(rate=rate, C=rate / (oc.RHO * math.sqrt(oc.G_ACC) * oc.D ** 1.5),
                win=win, m_open=m_open, m_end=m[-1] / m_open, t=t, m=m,
                leak=np.abs(m[closed] / m[0] - 1.0).max(),
                cons=np.abs((m + out) / m[0] - 1.0).max(), steps=int(h[-1, oc.STEP]))


def report(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    for d in detail:
        print(f"    {d}")
    return ok


def check_tank(work):
    print("[1] Flat tank: STL side walls and floor against mesh walls, E = 0.1 MPa, "
          f"maxCo {', '.join(f'{c:g}' for c in CO)}")
    res = {}
    for kind in ("mesh", "stl"):
        res[kind] = []
        for co in CO:
            p = make(work, f"{kind}{co:g}", kind, co)
            sec = oc.vc.run(p)
            r = result(p)
            r["sec"] = sec
            res[kind].append(r)
    mean = {k: float(np.mean([r["rate"] for r in v])) for k, v in res.items()}
    det = []
    for k, label in (("mesh", "mesh walls"), ("stl", "STL walls ")):
        rs = res[k]
        det.append(f"{label}: rates " + " / ".join(f"{r['rate']:.1f}" for r in rs)
                   + f" kg/s/m, mean {mean[k]:.1f} (C {mean[k] / (oc.RHO * math.sqrt(oc.G_ACC) * oc.D ** 1.5):.2f}), "
                   f"{sum(r['sec'] for r in rs):.0f} s")
    d = mean["stl"] / mean["mesh"] - 1.0
    leak = max(r["leak"] for r in res["stl"])
    cons = max(r["cons"] for r in res["stl"])
    det.append(f"mean STL rate {100 * d:+.1f} % against mesh walls (within {100 * TOL:g} %)")
    det.append(f"STL tank: mass change while closed {leak:.0e} (<= 1e-12), tank plus "
               f"discharged {cons:.0e} (<= 1e-10)")
    ok = abs(d) <= TOL and leak <= 1e-12 and cons <= 1e-10
    return report("Flat tank with STL walls", ok, det), mean


def check_hopper(work, flat):
    print(f"\n[2] Hopper, walls {HOPPER:g} deg to the horizontal, runs to empty")
    p = make(work, "hopper", "hopper", CO[0], T_HOP)
    sec = oc.vc.run(p)
    r = result(p)
    t, m = r["t"], r["m"]
    empty = np.where(m <= 0.02 * r["m_open"])[0]
    t_empty = t[empty[0]] - oc.T_OPEN if empty.size else float("nan")
    F = math.tan(math.radians(90.0 - HOPPER)) ** -0.35
    det = [f"rate {r['rate']:.1f} kg/s/m (C {r['C']:.2f}) over the first 3 to 12 % "
           f"discharged; flat tank mean {flat:.1f}: ratio {r['rate'] / flat:.2f} "
           f"(Rose-Tanaka hopper factor {F:.2f}; no pass band yet)",
           f"emptied to 2 % of {r['m_open']:.1f} kg/m in {t_empty:.2f} s after opening; at the "
           f"end {100 * r['m_end']:.1f} % left",
           f"mass while closed {r['leak']:.0e} (<= 1e-12), tank plus discharged {r['cons']:.0e} "
           f"(<= 1e-10); {r['steps']} steps, {sec:.0f} s"]
    np.savetxt(f"{work}/hopper_mass.dat", np.c_[t, m], header="t mass_kg_per_m")
    ok = np.isfinite(t_empty) and r["leak"] <= 1e-12 and r["cons"] <= 1e-10
    return report("Hopper runs to empty (rate recorded)", ok, det)


def main():
    work = sys.argv[1]
    if os.environ.get("CEL_CELL"):         # finer cells, e.g. CEL_CELL=0.005
        oc.CELL = oc.vc.CELL = float(os.environ["CEL_CELL"])
    which = sys.argv[2:] or ["tank", "hopper"]
    t0 = time.perf_counter()
    ok, flat = True, float("nan")
    if "tank" in which:
        good, mean = check_tank(work)
        ok &= good
        flat = mean["mesh"]
    if "flatstl" in which:            # one STL flat tank at maxCo 0.40, the
        p = make(work, "stl0.4", "stl", CO[0])   # hopper's twin (decision 0091)
        sec = oc.vc.run(p)
        r = result(p)
        flat = r["rate"]
        print(f"[1b] Flat STL tank, maxCo {CO[0]:g}: rate {flat:.1f} kg/s/m "
              f"(C {r['C']:.2f}), {sec:.0f} s")
    if "hopper" in which:
        ok &= check_hopper(work, flat)
    print(f"\nWP 4.5 checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
