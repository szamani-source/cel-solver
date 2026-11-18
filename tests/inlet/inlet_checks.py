"""WP 3.7a done-when: a solid inlet (materialProperties inlets, decision
0060) feeds fresh material into an empty box.

Set-up (10 mm cells, 2D): box W = 0.2 m by H = 0.3 m, walls and floor
friction 20 deg, a smooth lid with the inlet (w = 0.04 m) in its middle.
T0 Mohr-Coulomb (rho 1800, nu 0.3, phi 35, c 0, psi 0, no tension) at
E = 1 MPa with the T0 rule set (mu(I) placeholders, shear cap, rule 9,
flowImpedance, normalCap; decision 0050). Plain Mohr-Coulomb packs the
landing material to +4 % (overfilled cells, decision 0053), which is the
landing, not the inlet. Gravity from the
start, no damping. The inlet is closed until t = 0.05 s, then feeds material
at rho 1800 kg/m^3, stress free, at 0.5 m/s into the box until t = 0.35 s.
The run ends at t = 1.2 s.

Checks:
1. Closed inlet: before its start time nothing comes in (mass and fed mass
   exactly 0).
2. Exact rate: fed mass = rho s w (t - 0.05) to 1e-12 at every step while
   it is open, and rho s w 0.3 after it closes.
3. Mass balance: mass in the box equals the fed mass to 1e-12 throughout
   (nothing leaks, nothing is lost in the remap).
4. Empty cells carry no stress.
5. Density: in the falling stream (t = 0.35 s, above y = 0.15 m) every
   cell's material density rho/alpha is the inlet density to 1e-3 (exact
   until a cell passes through nearly full: the flowImpedance snap, decision
   0045, moves it by up to 6e-4 here); the
   deposit's mass over its material volume (sum alpha V) is the inlet
   density within 1 % (self weight compresses it by rho g h / K = 0.2 %,
   the landing packs a little more).
6. At rest: the largest speed at the end is below 0.02 m/s, and the
   deposit's centre of mass is in the lower third of the box.
7. Parallel: on 2 processors, split through the middle of the inlet, the
   fed mass and the mass in the box match serial to 1e-12 at the end.

Run by tests/test_inlet.sh after celVoidCheck is built:

    python3 inlet_checks.py WORKDIR
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

RHO, NU, PHI, E = 1800.0, 0.3, 35.0, 1.0e6
CELL = vc.CELL
W, H, WIN = 0.20, 0.30, 0.04
DELTA_W = 20.0
S_IN, T_START, T_STOP, T_END = 0.5, 0.05, 0.35, 1.2
MAXCO = 0.4
U_REST = 0.02
# history columns (voidCheck)
STEP, T, MASS, SEMPTY, UMAX, CY, CHARGED = 0, 1, 2, 8, 9, 6, 22


def block_mesh(path, nx, ny):
    """Three blocks along x so that the inlet is its own patch in the lid."""
    n_in = int(round(WIN / CELL))
    n_side = (nx - n_in) // 2
    xs = [0.0, n_side * CELL, (n_side + n_in) * CELL, W]
    Hd = ny * CELL
    verts = ([f"({x!r} 0 0)" for x in xs] + [f"({x!r} {Hd!r} 0)" for x in xs]
             + [f"({x!r} 0 {CELL!r})" for x in xs] + [f"({x!r} {Hd!r} {CELL!r})" for x in xs])
    ncell = [n_side, n_in, n_side]
    blocks = [f"    hex ({b} {b + 1} {5 + b} {4 + b} {8 + b} {9 + b} {13 + b} {12 + b}) "
              f"({ncell[b]} {ny} 1) simpleGrading (1 1 1)" for b in range(3)]

    def ymax(b):
        return f"({4 + b} {12 + b} {13 + b} {5 + b})"

    patches = {
        "left": "(0 8 12 4)",
        "right": "(3 7 15 11)",
        "bottom": " ".join(f"({b} {b + 1} {9 + b} {8 + b})" for b in range(3)),
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


def setup(work, name="box", nProcs=1):
    nx = int(round(W / CELL))
    ny = int(round(H / CELL))
    path = f"{work}/{name}"
    material = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
                "sigmaT 0;\nmu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650;\n"
                "shearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n")
    check = dict(fill=[0.0, 0.0, 0.0, 0.0], U0=[0.0, 0.0, 0.0], tEnd=T_END, maxCo=MAXCO,
                 writeTimes=[T_START, T_STOP])
    vc.make_case(path, nx, ny, material, check, nProcs=nProcs)
    block_mesh(path, nx, ny)
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    txt = txt.replace(f"rho {vc.RHO};", f"rho {RHO};")
    for w in ("left", "right", "bottom"):
        txt = txt.replace(f"    {w} {{ friction 0; }}", f"    {w} {{ friction {DELTA_W!r}; }}")
    txt += ("inlets\n{\n    inlet\n    {\n"
            f"        speed {S_IN!r};\n        startTime {T_START!r};\n"
            f"        endTime {T_STOP!r};\n        density {RHO!r};\n"
            "        friction 0;\n    }\n}\n")
    open(mp, "w").write(txt)
    return path


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    path = setup(work)
    sec = vc.run(path)
    h = np.loadtxt(f"{path}/history.dat")
    t, m, fed = h[:, T], h[:, MASS], h[:, CHARGED]
    rate = RHO * S_IN * WIN * CELL                  # kg/s through the inlet
    m_total = rate * (T_STOP - T_START)

    closed = t <= T_START + 1e-12
    before = max(np.abs(m[closed]).max(), np.abs(fed[closed]).max())
    expect = rate * np.clip(t - T_START, 0.0, T_STOP - T_START)
    err_rate = np.abs(fed - expect).max() / m_total
    err_bal = np.abs(m - fed).max() / m_total
    s_empty = h[:, SEMPTY].max()

    c35 = np.loadtxt(f"{path}/cells_{T_STOP!r}.dat")
    stream = (c35[:, 1] > 0.15) & (c35[:, 2] > 1e-6)
    err_stream = np.abs(c35[stream, 3] / c35[stream, 2] / RHO - 1.0).max()

    par = setup(work, "box2", nProcs=2)
    vc.run(par, nProcs=2)
    h2 = np.loadtxt(f"{par}/history.dat")
    err_par = max(abs(h2[-1, CHARGED] - fed[-1]), abs(h2[-1, MASS] - m[-1])) / m_total

    c = vc.cells(path)
    a, r = c[:, 2], c[:, 3]
    V = CELL ** 3
    rho_dep = (r * V).sum() / (a * V).sum()
    u_end = h[-1, UMAX]
    cy_end = h[-1, CY]
    full = a > 1 - 1e-9
    top = c[a > 1e-3, 1].max() + 0.5 * CELL if (a > 1e-3).any() else float("nan")

    print(f"Empty box W = {W} m, H = {H} m, 10 mm cells, inlet w = {WIN} m in the lid, "
          f"{S_IN} m/s from t = {T_START} to {T_STOP} s;")
    print(f"T0 Mohr-Coulomb E = {E / 1e6:g} MPa, walls and floor {DELTA_W:g} deg, "
          "T0 rule set (mu(I), flowImpedance, normalCap)")
    print(f"{int(h[-1, STEP])} steps to t = {t[-1]:.3f} s ({sec:.0f} s); fed "
          f"{fed[-1] / CELL:.4f} kg/m (expected {m_total / CELL:.4f}); {full.sum()} full cells, "
          f"deposit top y = {top:.3f} m\n")
    res = []

    def report(name, good, detail):
        res.append(good)
        print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}")

    report("1 closed inlet", before == 0.0,
           f"mass in before t = {T_START} s: {before:.1e} kg (0)")
    report("2 exact rate", err_rate <= 1e-12,
           f"fed mass against rho s w t: {err_rate:.0e} of the total (<= 1e-12)")
    report("3 mass balance", err_bal <= 1e-12,
           f"mass in the box against fed mass: {err_bal:.0e} (<= 1e-12)")
    report("4 empty cells", s_empty == 0.0,
           f"largest stress in empty cells {s_empty:.0e} Pa (0)")
    report("5 density", err_stream <= 1e-3 and abs(rho_dep / RHO - 1.0) <= 0.01,
           f"stream ({stream.sum()} cells) {err_stream:.0e} from the inlet density (<= 1e-3); "
           f"deposit mass over material volume {rho_dep:.1f} kg/m^3, "
           f"{100 * (rho_dep / RHO - 1.0):+.2f} % of the inlet density (within 1 %)")
    report("6 at rest", u_end <= U_REST and cy_end <= H / 3,
           f"largest speed at t = {t[-1]:.2f} s {u_end:.2e} m/s (<= {U_REST}), "
           f"centre of mass y = {cy_end:.3f} m (<= {H / 3:.3f})")
    report("7 parallel", err_par <= 1e-12,
           f"2 processors against serial, fed mass and mass in the box: {err_par:.0e} "
           "(<= 1e-12)")
    ok = all(res)
    print(f"\nInlet checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
