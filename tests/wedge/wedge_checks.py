"""WP 3.8a done-when: axisymmetric runs on an OpenFOAM wedge mesh (wedge
patches, decision 0061), without the free surface.

A wedge mesh is one cell thick in the hoop direction, 2 x 2.5 deg about the
y axis (gravity along -y); cells beside the axis are prisms (blockMesh
collapses the axis edge). Radial direction x, hoop direction z: s_zz of a
cell is its hoop stress.

Checks:
1. Lame thick cylinder (linear elastic, E 10 MPa, nu 0.3, plane strain:
   one cell high between symmetry planes): r = 0.1 to 0.3 m, 5 mm cells;
   the inner face moves out by 1e-5 m (1e-3 m/s for 0.01 s), the outer is
   held; damped to rest by 0.3 s. Radial, hoop and axial stress against
   Lame beyond two cells of the inner wall: within 3 % (5 mm) and 1 %
   (2.5 mm) of the largest |s_rr|; in those two cells the error is first order (decision
   0061): halving h must cut it to 0.6 of itself or less.
2. Lake at rest in a cylinder: a K0 column (s_rr = s_tt = K0 s_yy) between
   a smooth side wall and floor, R = 0.1 m, H = 0.6 m, axis cells included:
   largest speed over 1000 steps below 1e-10 m/s (discrete equilibrium of
   the hoop terms, as the Cartesian lake of WP 3.1).
3. Cylindrical Janssen tank: T0 Mohr-Coulomb (E 10 MPa), R = 0.2 m,
   H = 1.2 m, 10 mm cells, side wall and floor 20 deg, settled from zero
   stress (ramp 0.1 s, damping 250 1/s). Area-weighted mean s_v against
   Janssen for a circular bin, s_v = rho g R/(2 K mu) (1 - exp(-2 K mu z/R)),
   with K measured at the wall (wall s_rr / mean s_v, depth >= 2R): within
   5 % below depth 2R; K in 0.27 to 0.45 (decision 0019); weight = floor +
   wall forces within 1 %.
4. Two processors (split across the height) match serial to 1e-10 of
   rho g H over 300 steps of case 3.

Run by tests/test_wedge.sh after celSiloCheck is built:

    python3 wedge_checks.py WORKDIR
"""

import math
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "silo"))
sys.path.insert(0, os.path.join(HERE, "..", "lagrange"))
import silo_checks as sc  # noqa: E402  (make_case, run, summary, cells, walls)
import cases as cs  # noqa: E402

RHO, G_ACC = sc.RHO, sc.G_ACC
THETA = 2.5                       # half angle of the wedge [deg]
DELTA_W = sc.DELTA_W
K0 = sc.K0


def wedge_mesh(path, r0, r1, H, h, ny=None, inner="inner", top="top", bottom="bottom",
               top_type="patch", bottom_type="patch"):
    """Wedge about the y axis from r0 (0: the axis) to r1, height H, cells h
    (ny cells high if given)."""
    t = math.radians(THETA)
    c, s = math.cos(t), math.sin(t)
    nr = int(round((r1 - r0) / h))
    ny = ny or int(round(H / h))

    def P(r, y, sgn):
        return (r * c, y, sgn * r * s)
    if r0 > 0:
        v = [P(r0, 0, -1), P(r1, 0, -1), P(r1, H, -1), P(r0, H, -1),
             P(r0, 0, 1), P(r1, 0, 1), P(r1, H, 1), P(r0, H, 1)]
        hexv = "(0 1 2 3 4 5 6 7)"
        faces = {inner: "(0 4 7 3)", "right": "(1 2 6 5)", bottom: "(0 1 5 4)",
                 top: "(3 7 6 2)", "front": "(0 3 2 1)", "back": "(4 5 6 7)"}
    else:
        v = [P(0, 0, 0), P(r1, 0, -1), P(r1, H, -1), P(0, H, 0), P(r1, 0, 1), P(r1, H, 1)]
        hexv = "(0 1 2 3 0 4 5 3)"
        faces = {"right": "(1 2 5 4)", bottom: "(0 1 4 0)", top: "(3 5 2 3)",
                 "front": "(0 3 2 1)", "back": "(0 4 5 3)"}
    types = {"front": "wedge", "back": "wedge", top: top_type, bottom: bottom_type}
    pt = [f"    {n}\n    {{\n        type {types.get(n, 'patch')};\n        faces ({f});\n    }}"
          for n, f in faces.items()]
    if r0 == 0:
        pt.append("    axis\n    {\n        type empty;\n        faces ((0 3 3 0));\n    }")
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks ( hex {hexv} ({nr} {ny} 1) "
             "simpleGrading (1 1 1) );\nboundary\n(\n" + "\n".join(pt) + "\n);\n")


def set_walls(path, walls, material=sc.MC):
    cs.write(f"{path}/constant/materialProperties",
             cs.HEADER.format(cls="dictionary", obj="materialProperties")
             + f"rho {RHO};\n" + material + "walls\n{\n"
             + "".join(f"    {w} {{ friction {f!r}; }}\n" for w, f in walls.items()) + "}\n")


# ---------------------------------------------------------------- 1. Lame

A_IN, B_OUT, DELTA, U0, T_PUSH = 0.1, 0.3, 1e-5, 1e-3, 0.01
E_L, NU_L = 1e7, 0.3


def lame_exact(r):
    lam = E_L * NU_L / ((1 + NU_L) * (1 - 2 * NU_L))
    G = E_L / (2 * (1 + NU_L))
    A = DELTA * A_IN / (A_IN ** 2 - B_OUT ** 2)
    B = -A * B_OUT ** 2
    err, ett = A - B / r ** 2, A + B / r ** 2
    return ((lam + 2 * G) * err + lam * ett, lam * err + (lam + 2 * G) * ett,
            lam * (err + ett) + 0 * r)


def check_lame(work, h=0.005, tEnd=0.3):
    M = E_L * (1 - NU_L) / ((1 + NU_L) * (1 - 2 * NU_L))
    dt = 0.4 * h / math.sqrt(M / RHO)          # stable step (in-plane h)
    p = sc.make_case(f"{work}/lame", 0.2, h, False,
                     dict(mode="fixed", nSteps=int(math.ceil(tEnd / dt))), 0.0, 200.0)
    wedge_mesh(p, A_IN, B_OUT, h, h, ny=1, top_type="symmetryPlane",
               bottom_type="symmetryPlane")
    os.remove(f"{p}/constant/g")
    set_walls(p, {"right": 0.0}, f"type linearElastic;\nE {E_L!r};\nnu {NU_L!r};\n")
    cs.write(f"{p}/0/U", cs.HEADER.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
    inner {{ type uniformFixedValue; uniformValue table
        ((0 ({U0!r} 0 0)) ({T_PUSH!r} ({U0!r} 0 0)) ({T_PUSH * 1.0001!r} (0 0 0)) (100 (0 0 0))); }}
    right {{ type zeroGradient; }}
    top {{ type symmetryPlane; }}
    bottom {{ type symmetryPlane; }}
    front {{ type wedge; }}
    back {{ type wedge; }}
}}
""")
    sec = sc.run(p)
    s = sc.summary(p)
    c = sc.cells(p)
    r = c[:, 0]
    srr, stt, syy = lame_exact(r)
    ref = np.abs(srr).max()
    e = np.vstack([np.abs(c[:, 4] - srr), np.abs(c[:, 6] - stt), np.abs(c[:, 5] - syy)]) / ref
    near = r < A_IN + 2 * h
    return dict(h=h, steps=int(s["steps"]), sec=sec, ref=ref,
                e_wall=e[:, near].max(axis=1), e_rest=e[:, ~near].max(axis=1))


def report_lame(work):
    a = check_lame(work, 0.005)
    b = check_lame(work, 0.0025)
    rate = b["e_wall"].max() / a["e_wall"].max()
    ok = a["e_rest"].max() <= 0.03 and b["e_rest"].max() <= 0.01 and rate <= 0.6
    det = []
    for x in (a, b):
        det.append(f"{x['h'] * 1000:g} mm, {x['steps']} steps to 0.3 s ({x['sec']:.0f} s): "
                   f"beyond 2 cells of the inner wall s_rr / hoop / axial within "
                   f"{100 * x['e_rest'][0]:.2f} / {100 * x['e_rest'][1]:.2f} / "
                   f"{100 * x['e_rest'][2]:.2f} % of max |s_rr| {x['ref']:.0f} Pa "
                   f"(<= {3 if x['h'] > 0.003 else 1} %); "
                   f"in them {100 * x['e_wall'][0]:.1f} / {100 * x['e_wall'][1]:.1f} / "
                   f"{100 * x['e_wall'][2]:.1f} %")
    det.append(f"wall cells' largest error falls to {rate:.2f} of itself when h halves "
               "(<= 0.6: first order at the wall, decision 0061)")
    return sc.report("Lame thick cylinder (inner face pushed out 1e-5 m, outer held)", ok, det)


# ------------------------------------------------- 2. lake, 3. Janssen, 4.

def cylinder(work, name, R, H, friction, check, ramp=0.0, damping=0.0, nProcs=1):
    p = sc.make_case(f"{work}/{name}", R, H, friction, check, ramp, damping, nProcs)
    wedge_mesh(p, 0.0, R, H, sc.CELL)
    fr = DELTA_W if friction else 0.0
    set_walls(p, {"right": fr, "bottom": fr})
    return p


def check_lake(work):
    R, H = 0.1, 0.6
    p = cylinder(work, "lake", R, H, False, dict(mode="fixed", nSteps=1000, initK0=K0))
    sec = sc.run(p)
    s = sc.summary(p)
    ok = s["vmax"] <= 1e-10
    return sc.report("Lake at rest in a cylinder: K0 column, smooth wall and floor, 1000 steps",
                     ok, [f"largest speed {s['vmax']:.1e} m/s (<= 1e-10; Cartesian 3e-17), "
                          f"dt = {1e6 * s['dt']:.1f} us, {sec:.1f} s"])


def janssen_circ(z, R, K, mu):
    a = 2.0 * mu * K / R
    return RHO * G_ACC / a * (1.0 - np.exp(-a * z))


def check_janssen(work):
    R, H = 0.2, 1.2
    p = cylinder(work, "janssen", R, H, True,
                 dict(mode="settle", keTol=sc.KE_TOL, checkEvery=20, tMax=1.0,
                      continueFraction=0.1), sc.T_RAMP, sc.ALPHA)
    sec = sc.run(p)
    s = sc.summary(p)
    good, det = sc.settle_report(s, RHO * G_ACC * H)
    mu = math.tan(math.radians(DELTA_W))
    c = sc.cells(p)
    nr, ny = int(round(R / sc.CELL)), int(round(H / sc.CELL))
    g = {k: c[:, i].reshape(ny, nr) for i, k in
         enumerate(("x", "y", "ux", "uy", "sxx", "syy", "szz", "sxy", "rho", "ep"))}
    z = H - g["y"][:, 0]
    r = g["x"][0]
    w = r / r.sum()                                    # area weights ~ r
    sv = -(g["syy"] * w).sum(axis=1)
    deep = z >= 2 * R
    sh = -(1.5 * g["sxx"][:, -1] - 0.5 * g["sxx"][:, -2])   # extrapolated to the wall
    K_wall = sh[deep].sum() / sv[deep].sum()
    err = np.abs(janssen_circ(z, R, K_wall, mu)[deep] / sv[deep] - 1.0).max()
    ok1 = err <= 0.05 and 0.27 <= K_wall <= 0.45
    det += [f"K at the wall (wall s_rr / mean s_v, depth >= 2R) {K_wall:.4f} (0.27 to 0.45)",
            f"area-weighted mean s_v against circular Janssen with that K, depth >= 2R: "
            f"largest error {100 * err:.2f} % (<= 5 %)"]
    for frac in (0.25, 0.5, 1.0):
        i = int(np.argmin(np.abs(z - frac * H)))
        det.append(f"  depth {z[i]:.3f} m: s_v {sv[i]:8.1f} Pa, Janssen "
                   f"{janssen_circ(z[i], R, K_wall, mu):8.1f}, hydrostatic "
                   f"{RHO * G_ACC * z[i]:8.1f}")
    hoop = g["szz"][deep] / g["sxx"][deep]
    det.append(f"hoop / radial stress below depth 2R: {hoop.min():.3f} to {hoop.max():.3f}")
    fb, fr, wt = s["force_bottom"], s["force_right"], s["weight"]
    ef = abs(fb + fr - wt) / wt
    ok2 = ef <= 0.01
    det.append(f"weight = floor {fb / wt:.4f} + wall {fr / wt:.4f} of it: imbalance "
               f"{100 * ef:.3f} % (<= 1 %)")
    det.append(f"wall clock {sec:.0f} s")
    return sc.report("Cylindrical Janssen tank (R 0.2 m, H 1.2 m, wall and floor 20 deg)",
                     good and ok1 and ok2, det)


def check_parallel(work):
    R, H = 0.2, 1.2
    chk = dict(mode="fixed", nSteps=300)
    a = cylinder(work, "par1", R, H, True, chk, sc.T_RAMP, sc.ALPHA)
    b = cylinder(work, "par2", R, H, True, chk, sc.T_RAMP, sc.ALPHA, nProcs=2)
    sc.run(a)
    sc.run(b, 2)
    ca, cb = sc.cells(a), sc.cells(b, 2)
    d = np.abs(ca[:, 4:8] - cb[:, 4:8]).max() / (RHO * G_ACC * H)
    return sc.report("Two processors against serial (cylinder, 300 steps)", d <= 1e-10,
                     [f"largest stress difference {d:.1e} of rho g H (<= 1e-10)"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    res = [report_lame(work), check_lake(work), check_janssen(work),
           check_parallel(work)]
    ok = all(res)
    print(f"\nWedge checks: {'ALL PASS' if ok else 'FAILED'} ({time.perf_counter() - t0:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
