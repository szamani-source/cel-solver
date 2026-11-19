"""WP 3.7b done-when: the T0 tank discharges while an inlet in its lid feeds
it (inlets, decision 0060).

Set-up as the outlet test (tests/outlet/outlet_checks.py, decision 0033) with
the T0 rule set (decision 0050): tank W = 0.4 m, filled to H = 0.6 m, 10 mm
cells, T0 Mohr-Coulomb with the mu(I) placeholders at E = 0.1 MPa (fast; the
rate at T0's E is the tutorial's), walls and floor 20 deg, outlet D = 0.06 m
in the middle of the floor. The lid is 0.1 m above the fill; an inlet of
the outlet's width sits in the middle of the lid. Stress-free fill, gravity
ramp 0.1 s, damping 250 1/s; at t = 0.3 s the outlet opens and the inlet
starts feeding (stress free, rho 1800) at a fixed speed.

Three runs: no feed (reference rate R), feed F = R/2, feed F = R (speed
F / (rho D)), each to 1.3 s.

Checks:
1. Mass balance, every run: tank mass = initial + fed - discharged to 1e-10
   of the initial mass, and nothing fed or discharged before the opening.
2. The discharge rate does not depend on the feed (Beverloo: the rate is set
   at the outlet, not by the head): the rate out with F = R/2 and F = R
   within 10 % of R (steady window 0.35 to 1.0 s after opening).
3. Matched feed: the tank mass stays within 3 % of its mass at the opening
   over the window (the level holds).
4. No stress in empty cells, every run.

Run by tests/test_fed_silo.sh after celVoidCheck is built:

    python3 fed_silo_checks.py WORKDIR
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
sys.path.insert(0, os.path.join(HERE, "..", "lagrange"))
import void_checks as vc  # noqa: E402  (make_case, run)
import cases as cs  # noqa: E402  (HEADER, write)

RHO, NU, PHI, E = 1800.0, 0.3, 35.0, 1.0e5
CELL = vc.CELL
W, H, D = 0.40, 0.60, 0.06
HEADROOM = 10                      # cells of void between fill and lid
DELTA_W = 20.0
T_RAMP, DAMPING, T_OPEN = 0.1, 250.0, 0.3
T_DIS = 1.0
WINDOW = (0.35, 1.0)               # after opening
MAXCO = 0.4
RULES = ("mu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650;\nshearCap yes;\n"
         "spreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n")
# history columns (voidCheck)
STEP, T, MASS, SEMPTY, OUT, FED = 0, 1, 2, 8, 14, 22


def block_mesh(path, nx, ny):
    """Three blocks along x: the outlet is the middle block's floor, the
    inlet the middle block's lid."""
    n_mid = int(round(D / CELL))
    n_side = (nx - n_mid) // 2
    xs = [0.0, n_side * CELL, (n_side + n_mid) * CELL, W]
    Hd = ny * CELL
    verts = ([f"({x!r} 0 0)" for x in xs] + [f"({x!r} {Hd!r} 0)" for x in xs]
             + [f"({x!r} 0 {CELL!r})" for x in xs] + [f"({x!r} {Hd!r} {CELL!r})" for x in xs])
    ncell = [n_side, n_mid, n_side]
    blocks = [f"    hex ({b} {b + 1} {5 + b} {4 + b} {8 + b} {9 + b} {13 + b} {12 + b}) "
              f"({ncell[b]} {ny} 1) simpleGrading (1 1 1)" for b in range(3)]

    def ymin(b):
        return f"({b} {b + 1} {9 + b} {8 + b})"

    def ymax(b):
        return f"({4 + b} {12 + b} {13 + b} {5 + b})"

    patches = {
        "left": "(0 8 12 4)",
        "right": "(3 7 15 11)",
        "bottom": ymin(0) + " " + ymin(2),
        "outlet": ymin(1),
        "top": ymax(0) + " " + ymax(2),
        "inlet": ymax(1),
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


def setup(work, name, feed):
    """feed: mass rate per metre of depth [kg/s/m] (0: inlet never opens)."""
    nx = int(round(W / CELL))
    ny = int(math.ceil(H / CELL - 1e-9)) + HEADROOM
    path = f"{work}/{name}"
    material = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
                "sigmaT 0;\n" + RULES)
    check = dict(fill=[0.0, 0.0, W, H], U0=[0.0, 0.0, 0.0], tEnd=T_OPEN + T_DIS,
                 maxCo=MAXCO, writeTimes=[T_OPEN])
    vc.make_case(path, nx, ny, material, check)
    block_mesh(path, nx, ny)
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    txt = txt.replace(f"rho {vc.RHO};", f"rho {RHO};")
    for w in ("left", "right", "bottom"):
        txt = txt.replace(f"    {w} {{ friction 0; }}", f"    {w} {{ friction {DELTA_W!r}; }}")
    txt += f"outlets\n{{\n    outlet {{ friction {DELTA_W!r}; openTime {T_OPEN!r}; }}\n}}\n"
    speed = float(feed) / (RHO * D)
    start = T_OPEN if feed > 0 else 1e30
    txt += (f"inlets\n{{\n    inlet\n    {{\n        speed {speed!r};\n"
            f"        startTime {start!r};\n        density {RHO!r};\n    }}\n}}\n")
    open(mp, "w").write(txt)
    open(f"{path}/system/controlDict", "a").write(
        f"gravityRamp {T_RAMP!r};\ndamping {DAMPING!r};\n")
    return path


def rate(t, y):
    """Least-squares slope of y over the steady window after opening."""
    sel = (t >= T_OPEN + WINDOW[0]) & (t <= T_OPEN + WINDOW[1])
    return np.polyfit(t[sel], y[sel], 1)[0]


def run_case(work, name, feed):
    path = setup(work, name, feed)
    sec = vc.run(path)
    h = np.loadtxt(f"{path}/history.dat")
    # per metre of depth (the mesh is one cell, CELL, thick)
    t, m, out, fed = h[:, T], h[:, MASS] / CELL, h[:, OUT] / CELL, h[:, FED] / CELL
    m0 = m[0]
    closed = t <= T_OPEN + 1e-12
    before = max(np.abs(out[closed]).max(), np.abs(fed[closed]).max())
    bal = np.abs((m - fed + out) / m0 - 1.0).max()
    i_open = int(np.where(closed)[0][-1])
    sel = (t >= T_OPEN + WINDOW[0]) & (t <= T_OPEN + WINDOW[1])
    return dict(name=name, feed=feed, sec=sec, steps=int(h[-1, STEP]), t=t, m=m,
                before=before, bal=bal, r_out=rate(t, out), r_tank=rate(t, m),
                m_open=m[i_open], dm_win=np.abs(m[sel] / m[i_open] - 1.0).max(),
                s_empty=h[:, SEMPTY].max(), fed_end=fed[-1], out_end=out[-1])


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    ref = run_case(work, "unfed", 0.0)
    R = float(ref["r_out"])
    half = run_case(work, "half", 0.5 * R)
    full = run_case(work, "matched", R)
    runs = (ref, half, full)

    print(f"T0 tank W = {W} m, H = {H} m, outlet D = {D} m, 10 mm cells, Mohr-Coulomb "
          f"phi {PHI:g}, E = {E / 1e6:g} MPa, T0 rule set;")
    print(f"lid {HEADROOM * CELL:g} m above the fill, inlet D wide in its middle; outlet "
          f"opens and feed starts at {T_OPEN} s; rates over {WINDOW[0]} to {WINDOW[1]} s "
          "after opening (kg/s/m)")
    for r in runs:
        print(f"  {r['name']:8s} feed {r['feed']:6.2f}: out {r['r_out']:6.2f}, tank "
              f"{r['r_tank']:+7.2f}; fed {r['fed_end']:.2f}, discharged {r['out_end']:.2f} "
              f"kg/m; {r['steps']} steps ({r['sec']:.0f} s)")
    print()
    res = []

    def report(name, good, detail):
        res.append(good)
        print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}")

    bal = max(r["bal"] for r in runs)
    before = max(r["before"] for r in runs)
    report("1 mass balance", bal <= 1e-10 and before == 0.0,
           f"tank = initial + fed - discharged to {bal:.0e} (<= 1e-10); fed or discharged "
           f"before opening {before:.0e} (0)")
    dev = [r["r_out"] / R - 1.0 for r in (half, full)]
    report("2 rate independent of the feed", max(abs(d) for d in dev) <= 0.10,
           f"rate out {100 * dev[0]:+.1f} % (half feed), {100 * dev[1]:+.1f} % (matched) "
           f"of the unfed {R:.1f} kg/s/m (within 10 %)")
    report("3 matched feed", full["dm_win"] <= 0.03,
           f"tank mass within {100 * full['dm_win']:.2f} % of its mass at opening (<= 3 %)")
    se = max(r["s_empty"] for r in runs)
    report("4 empty cells", se == 0.0, f"largest stress in empty cells {se:.0e} Pa (0)")
    ok = all(res)
    print(f"\nFed silo checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
