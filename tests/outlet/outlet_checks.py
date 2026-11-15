"""WP 3.3 done-when: the T0 tank in celFoam's cycle with an outlet that is a
closed frictional wall while the material settles and opens at a set time
(materialProperties outlets, decision 0033; prototype rule 10 and 10b,
decisions 0026 and 0027). Shear cap and rule 9 on (decision 0032), neither
rule 12 nor 13 (decision 0029).

Set-up as the prototype's outlet check (prototype/discharge2d.py, 10 mm):
tank W = 0.4 m, material to H = 0.6 m with 2 cells of void above, T0
Mohr-Coulomb (rho 1800, nu 0.3, phi 35, c 0, psi 0, no tension) at
E = 0.1 MPa, walls and floor friction 20 deg, outlet D = 0.06 m centred in
the floor. Stress-free fill, gravity ramp 0.1 s, damping 250 1/s; the outlet
opens at t = 0.3 s and damping is switched off.

Checks:
1. Nothing leaves before the opening time (discharged mass exactly 0) and the
   mass in the tank is constant to 1e-12 until then.
2. Material leaves after it opens: discharged mass > 0 within 0.05 s.
3. Mass in tank plus discharged mass constant to 1e-10 through the run.
4. No stress in empty cells.
5. Steady rate (least squares, tank mass from 97 to 88 percent of the mass at
   opening, as the prototype's check): Beverloo C at least 0.8 (the outlet is
   not held shut, rule 10b) and within 20 percent of the prototype at the
   same settings without rules 12 and 13 (PROTO_RATE).

The rate at E = 0.1, 1 and 10 MPa is WP 3.4.

Run by tests/test_outlet.sh after celVoidCheck is built:

    python3 outlet_checks.py WORKDIR
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
G_ACC = 9.81
CELL = vc.CELL
W, H, D = 0.40, 0.60, 0.06
HEADROOM = 2                       # cells of void above the fill
DELTA_W = 20.0
T_RAMP, DAMPING, T_OPEN = 0.1, 250.0, 0.3
T_DIS = 1.0                        # run this long after opening
WINDOW = (0.97, 0.88)
# prototype/discharge2d.py run_case(cell=0.01, H=0.6, D=0.06, E=1e5,
# window=(0.97, 0.88), new_rules=False), 2026-11-15: neither rule 12 nor 13, as here
PROTO_RATE = 81.35
C_MIN = 0.8
# history columns (voidCheck)
STEP, T, MASS, SEMPTY, KE, OUT = 0, 1, 2, 8, 12, 14


def block_mesh(path, nx, ny):
    """Three blocks along x so that the outlet is its own patch."""
    n_out = int(round(D / CELL))
    n_side = (nx - n_out) // 2
    xs = [0.0, n_side * CELL, (n_side + n_out) * CELL, W]
    Hd = ny * CELL
    verts = ([f"({x!r} 0 0)" for x in xs] + [f"({x!r} {Hd!r} 0)" for x in xs]
             + [f"({x!r} 0 {CELL!r})" for x in xs] + [f"({x!r} {Hd!r} {CELL!r})" for x in xs])
    ncell = [n_side, n_out, n_side]
    blocks = [f"    hex ({b} {b + 1} {5 + b} {4 + b} {8 + b} {9 + b} {13 + b} {12 + b}) "
              f"({ncell[b]} {ny} 1) simpleGrading (1 1 1)" for b in range(3)]

    def ymin(b):
        return f"({b} {b + 1} {9 + b} {8 + b})"

    patches = {
        "left": "(0 8 12 4)",
        "right": "(3 7 15 11)",
        "bottom": ymin(0) + " " + ymin(2),
        "outlet": ymin(1),
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


def setup(work, name="tank"):
    nx = int(round(W / CELL))
    ny = int(math.ceil(H / CELL - 1e-9)) + HEADROOM
    path = f"{work}/{name}"
    material = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
                "sigmaT 0;\nshearCap yes;\nspreadLayers yes;\n")
    check = dict(fill=[0.0, 0.0, W, H], U0=[0.0, 0.0, 0.0], tEnd=T_OPEN + T_DIS,
                 maxCo=0.4)
    vc.make_case(path, nx, ny, material, check)
    block_mesh(path, nx, ny)
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    txt = txt.replace(f"rho {vc.RHO};", f"rho {RHO};")
    for w in ("left", "right", "bottom"):
        txt = txt.replace(f"    {w} {{ friction 0; }}", f"    {w} {{ friction {DELTA_W!r}; }}")
    txt += f"outlets\n{{\n    outlet {{ friction {DELTA_W!r}; openTime {T_OPEN!r}; }}\n}}\n"
    open(mp, "w").write(txt)
    cd = f"{path}/system/controlDict"
    open(cd, "a").write(f"gravityRamp {T_RAMP!r};\ndamping {DAMPING!r};\n")
    return path


def steady_rate(t, m, m_open):
    sel = (m <= WINDOW[0] * m_open) & (m >= WINDOW[1] * m_open)
    if sel.sum() < 5:
        return float("nan"), (float("nan"), float("nan"))
    return -np.polyfit(t[sel], m[sel], 1)[0], (t[sel][0], t[sel][-1])


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    path = setup(work)
    sec = vc.run(path)
    h = np.loadtxt(f"{path}/history.dat")
    # per metre of depth (the mesh is one cell, CELL, thick)
    t, m, out = h[:, T], h[:, MASS] / CELL, h[:, OUT] / CELL
    m0 = m[0]
    closed = t <= T_OPEN + 1e-12
    leak = np.abs(out[closed]).max()
    dm_closed = np.abs(m[closed] / m0 - 1.0).max()
    i_open = int(np.where(closed)[0][-1])
    m_open = m[i_open]
    first = np.where(out > 0.0)[0]
    t_first = t[first[0]] - T_OPEN if first.size else float("inf")
    cons = np.abs((m + out) / m0 - 1.0).max()
    s_empty = h[:, SEMPTY].max()
    tau = t[i_open:] - T_OPEN
    rate, win = steady_rate(tau, m[i_open:], m_open)
    C = rate / (RHO * math.sqrt(G_ACC) * D ** 1.5) if np.isfinite(rate) else float("nan")
    m_end = m[-1] / m_open

    print(f"T0 tank W = {W} m, H = {H} m, outlet D = {D} m, 10 mm cells, Mohr-Coulomb "
          f"c = 0, phi {PHI:g}, E = {E / 1e6:g} MPa,")
    print(f"walls and floor {DELTA_W:g} deg, ramp {T_RAMP} s, damping {DAMPING:g} 1/s, outlet "
          f"opens at {T_OPEN} s; shear cap and rule 9 on, neither rule 12 nor 13")
    print(f"{int(h[-1, STEP])} steps to t = {t[-1]:.3f} s ({sec:.0f} s); tank mass at the end "
          f"{m_end:.3f} of the mass at opening\n")
    res = []

    def report(name, good, detail):
        res.append(good)
        print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}")

    report("1 closed outlet", leak == 0.0 and dm_closed <= 1e-12,
           f"discharged {leak:.1e} kg/m (0), tank mass change {dm_closed:.0e} (<= 1e-12)")
    report("2 opens", t_first <= 0.05,
           f"first material out {t_first * 1e3:.1f} ms after opening (<= 50 ms)")
    report("3 mass balance", cons <= 1e-10,
           f"tank plus discharged constant to {cons:.0e} (<= 1e-10)")
    report("4 empty cells", s_empty == 0.0,
           f"largest stress in empty cells {s_empty:.0e} Pa (0)")
    good = np.isfinite(rate) and C >= C_MIN and abs(rate / PROTO_RATE - 1.0) <= 0.2
    report("5 rate", good,
           f"{rate:.1f} kg/s/m over t = {win[0]:.3f} to {win[1]:.3f} s after opening, "
           f"C = {C:.3f} (>= {C_MIN}); prototype {PROTO_RATE} kg/s/m, "
           f"{100 * (rate / PROTO_RATE - 1.0):+.1f} % (within 20 %)")
    ok = all(res)
    print(f"\nOutlet (C++) checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
