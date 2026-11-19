"""WP 3.8c done-when: a cylindrical silo discharging through a central
circular outlet, run as an axisymmetric wedge (decision 0058), against the
3D Beverloo correlation and against the same tank as a coarse 3D quarter
model.

Tank R = 0.1 m (D_t = 0.2 m), filled to H = 0.2 m, central outlet D = 0.06
m, 10 mm cells. T0 Mohr-Coulomb with the T0 rule set (decision 0050) at
E = 0.1 MPa, wall and floor 20 deg. Stress-free fill, gravity ramp 0.1 s,
damping 250 1/s; the outlet opens at 0.3 s (damping off).

Checks:
1. Closed outlet: nothing leaves before 0.3 s (0 exactly); tank plus
   discharged mass constant to 1e-10 throughout.
2. Beverloo (3D): steady mass flow W (least squares, tank mass 90 to 70 %
   of the mass at opening), full circle (wedge mass times 2 pi / wedge
   angle), against W = C rho sqrt(g) (D - k d)^2.5 with k = 1.5, d = 1 mm:
   C within 0.45 to 0.75 (Beverloo et al. 1961, C 0.58 for bulk density;
   Nedderman 1992: 0.55 to 0.65; widened by the 20 % rate uncertainty
   accepted for the 2D silo, decision 0045).
3. Quarter model: the same tank as a 3D quarter (symmetry planes x = 0 and
   z = 0, STL cylinder wall, 10 mm cells, square quarter outlet 0.03 x 0.03 m,
   so 3 cells across like the wedge's outlet radius; with 2 cells, 0.02 m,
   the quarter reads C 0.50 against the wedge's 0.69, decision 0062): C of the wedge within 10 % of C of
   the quarter model, both from their own outlet area.
4. Quarter model by the hydraulic diameter (decision 0090): Beverloo for a
   non-circular outlet is W = C rho (4/pi) A_k sqrt(g D_h,k), A and D_h = 4A/P
   taken after the empty annulus k d (Zhao et al., arXiv 2610.02543, eq. 15).
   For the square, side a: C = pi W / (4 rho sqrt(g)
   (a - k d)^2.5). The equal-area circle of check 3 reads the square 6 to 7 %
   low by construction; with D_h the wedge and the quarter agree within 3 %.

Run by tests/test_wedge_silo.sh after celVoidCheck is built:

    python3 wedge_silo_checks.py WORKDIR
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("void", "silo", "lagrange", "wedge", "wedgevoid"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import void_checks as vc  # noqa: E402

RHO, NU, PHI, E = 1800.0, 0.3, 35.0, 1.0e5
G_ACC = 9.81
CELL = float(os.environ.get("CEL_CELL", "0.01"))
R, H, R_OUT = 0.10, 0.20, 0.03
HEADROOM = 2
DELTA_W = 20.0
T_RAMP, DAMPING, T_OPEN, T_DIS = 0.1, 250.0, 0.3, 3.5
WINDOW = (0.90, 0.70)
THETA = 2.5                                   # wedge half angle [deg]
K_BEV, D_GRAIN = 1.5, 0.001
C_BAND = (0.45, 0.75)
MATERIAL = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
            "sigmaT 0;\nmu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650;\n"
            "shearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n")
# history columns (voidCheck)
STEP, T, MASS, OUT = 0, 1, 2, 14


def wedge_silo_mesh(path, ny):
    """Wedge about y: block 0 from the axis to R_OUT (axis prisms, its floor
    is the outlet), block 1 from R_OUT to R."""
    t = math.radians(THETA)
    c, s = math.cos(t), math.sin(t)
    Hm = ny * CELL

    def P(r, y, sg):
        return (r * c, y, sg * r * s)
    v = [P(0, 0, 0), P(R_OUT, 0, -1), P(R_OUT, Hm, -1), P(0, Hm, 0),   # 0-3
         P(R_OUT, 0, 1), P(R_OUT, Hm, 1),                              # 4-5
         P(R, 0, -1), P(R, Hm, -1), P(R, 0, 1), P(R, Hm, 1)]           # 6-9
    n0, n1 = int(round(R_OUT / CELL)), int(round((R - R_OUT) / CELL))
    blocks = (f"    hex (0 1 2 3 0 4 5 3) ({n0} {ny} 1) simpleGrading (1 1 1)\n"
              f"    hex (1 6 7 2 4 8 9 5) ({n1} {ny} 1) simpleGrading (1 1 1)\n")
    patches = {
        "outlet": ("patch", "(0 1 4 0)"),
        "bottom": ("patch", "(1 6 8 4)"),
        "right": ("patch", "(6 7 9 8)"),
        "top": ("patch", "(3 5 2 3) (2 5 9 7)"),
        "front": ("wedge", "(0 3 2 1) (1 2 7 6)"),
        "back": ("wedge", "(0 4 5 3) (4 8 9 5)"),
    }
    pt = [f"    {n}\n    {{\n        type {ty};\n        faces ({f});\n    }}"
          for n, (ty, f) in patches.items()]
    pt.append("    axis\n    {\n        type empty;\n        faces ((0 3 3 0));\n    }")
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n{blocks});\n"
             "boundary\n(\n" + "\n".join(pt) + "\n);\n")


def finish_case(path, walls):
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    i = txt.index("walls\n{")
    txt = txt[:i] + "walls\n{\n" + "".join(
        f"    {w} {{ friction {f!r}; }}\n" for w, f in walls.items()) + "}\n"
    txt = txt.replace(f"rho {vc.RHO};", f"rho {RHO};")
    txt += f"outlets\n{{\n    outlet {{ friction {DELTA_W!r}; openTime {T_OPEN!r}; }}\n}}\n"
    open(mp, "w").write(txt)
    open(f"{path}/system/controlDict", "a").write(
        f"gravityRamp {T_RAMP!r};\ndamping {DAMPING!r};\n")


def wedge_case(work, name="wedge"):
    ny = int(round(H / CELL)) + HEADROOM
    p = vc.make_case(f"{work}/{name}", 1, 1, MATERIAL,
                     dict(fill=[-1.0, 0.0, 1.0, H], U0=[0.0, 0.0, 0.0],
                          tEnd=T_OPEN + T_DIS, maxCo=0.4))
    wedge_silo_mesh(p, ny)
    finish_case(p, {"right": DELTA_W, "bottom": DELTA_W, "top": 0.0})
    return p


Q_OUT = 0.03                 # quarter model: square outlet side (quarter)
Q_BOX = 0.12                 # quarter model: mesh extent in x and z


def tube_stl(r0, r1, y0, y1, n=128):
    """ASCII STL of the closed tube r0 < r < r1, y0 < y < y1 about the y axis
    (outward normals)."""
    a = np.linspace(0.0, 2 * math.pi, n + 1)
    out = ["solid tube"]

    def tri(p, q, w):
        p, q, w = map(np.asarray, (p, q, w))
        nn = np.cross(q - p, w - p)
        nn /= np.linalg.norm(nn)
        out.append(f"facet normal {nn[0]:.9e} {nn[1]:.9e} {nn[2]:.9e}\nouter loop")
        for v in (p, q, w):
            out.append(f"vertex {v[0]:.12e} {v[1]:.12e} {v[2]:.12e}")
        out.append("endloop\nendfacet")

    def P(r, ang, y):
        return (r * math.cos(ang), y, -r * math.sin(ang))
    for i in range(n):
        a0, a1 = a[i], a[i + 1]
        # outer wall (normal out), inner wall (normal towards the axis)
        tri(P(r1, a0, y0), P(r1, a1, y0), P(r1, a1, y1))
        tri(P(r1, a0, y0), P(r1, a1, y1), P(r1, a0, y1))
        tri(P(r0, a0, y0), P(r0, a1, y1), P(r0, a1, y0))
        tri(P(r0, a0, y0), P(r0, a0, y1), P(r0, a1, y1))
        # bottom (normal -y) and top (normal +y) rings
        tri(P(r0, a0, y0), P(r1, a1, y0), P(r1, a0, y0))
        tri(P(r0, a0, y0), P(r0, a1, y0), P(r1, a1, y0))
        tri(P(r0, a0, y1), P(r1, a0, y1), P(r1, a1, y1))
        tri(P(r0, a0, y1), P(r1, a1, y1), P(r0, a1, y1))
    out.append("endsolid tube")
    return "\n".join(out) + "\n"


def quarter_mesh(path, ny):
    """Quarter box x, z in [0, Q_BOX], 2 x 2 blocks so that the floor under
    x, z < Q_OUT is the outlet; symmetry planes at x = 0 and z = 0."""
    Hm = ny * CELL
    xs = [0.0, Q_OUT, Q_BOX]
    n = [int(round(Q_OUT / CELL)), int(round((Q_BOX - Q_OUT) / CELL))]

    def vid(i, j, k):            # i: x, j: z, k: 0 bottom / 1 top
        return k * 9 + j * 3 + i
    verts = [f"    ({xs[i]!r} {k * Hm!r} {xs[j]!r})"
             for k in range(2) for j in range(3) for i in range(3)]
    blocks, pat = [], {"outlet": [], "bottom": [], "top": [], "symX": [], "symZ": [],
                       "outer": []}
    for j in range(2):
        for i in range(2):
            b = [vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j, 1), vid(i, j, 1),
                 vid(i, j + 1, 0), vid(i + 1, j + 1, 0), vid(i + 1, j + 1, 1),
                 vid(i, j + 1, 1)]
            blocks.append(f"    hex ({' '.join(map(str, b))}) ({n[i]} {ny} {n[j]}) "
                          "simpleGrading (1 1 1)")
            bot = f"({b[0]} {b[4]} {b[5]} {b[1]})"
            pat["outlet" if i == 0 and j == 0 else "bottom"].append(bot)
            pat["top"].append(f"({b[3]} {b[2]} {b[6]} {b[7]})")
            if i == 0:
                pat["symX"].append(f"({b[0]} {b[3]} {b[7]} {b[4]})")
            if j == 0:
                pat["symZ"].append(f"({b[0]} {b[1]} {b[2]} {b[3]})")
            if i == 1:
                pat["outer"].append(f"({b[1]} {b[5]} {b[6]} {b[2]})")
            if j == 1:
                pat["outer"].append(f"({b[4]} {b[7]} {b[6]} {b[5]})")
    ty = {"symX": "symmetryPlane", "symZ": "symmetryPlane"}
    pt = [f"    {k}\n    {{\n        type {ty.get(k, 'patch')};\n        faces ({' '.join(f)});"
          "\n    }" for k, f in pat.items()]
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + "scale 1;\nvertices\n(\n" + "\n".join(verts) + "\n);\nblocks\n(\n"
             + "\n".join(blocks) + "\n);\nboundary\n(\n" + "\n".join(pt) + "\n);\n")
    cs.write(f"{path}/constant/triSurface/tube.stl",
             tube_stl(R, R + 0.1, -0.05, Hm + 0.05))      # covers the box corner
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n    tube\n    {\n        surface tube.stl;\n"
             f"        friction {DELTA_W!r};\n    }}\n}}\n")


def quarter_case(work, name="quarter"):
    ny = int(round(H / CELL)) + HEADROOM
    p = vc.make_case(f"{work}/{name}", 1, 1, MATERIAL,
                     dict(fill=[-1.0, 0.0, 1.0, H], U0=[0.0, 0.0, 0.0],
                          tEnd=T_OPEN + T_DIS, maxCo=0.4))
    quarter_mesh(p, ny)
    finish_case(p, {"bottom": DELTA_W, "top": 0.0, "outer": DELTA_W})
    return p


def steady_rate(t, m, m_open):
    sel = (m <= WINDOW[0] * m_open) & (m >= WINDOW[1] * m_open)
    if sel.sum() < 5:
        return float("nan"), (float("nan"), float("nan"))
    return -np.polyfit(t[sel], m[sel], 1)[0], (t[sel][0], t[sel][-1])


def beverloo_C(W, D):
    return W / (RHO * math.sqrt(G_ACC) * (D - K_BEV * D_GRAIN) ** 2.5)


def beverloo_C_square(W, a):
    """Beverloo C for a square outlet of side a by the hydraulic diameter."""
    return math.pi * W / (4 * RHO * math.sqrt(G_ACC) * (a - K_BEV * D_GRAIN) ** 2.5)


def analyse(path, scale):
    h = np.loadtxt(f"{path}/history.dat")
    t, m, out = h[:, T], h[:, MASS] * scale, h[:, OUT] * scale
    i0 = int(np.searchsorted(t, T_OPEN))
    leak = abs(out[:i0]).max()
    cons = np.abs(m + out - (m[0] + out[0])).max() / m[0]
    W, win = steady_rate(t[i0:] - T_OPEN, m[i0:], m[i0])
    return dict(leak=leak, cons=cons, W=W, win=win, m0=m[0], end=m[-1] / m[i0])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    fails = []

    def check(ok, text):
        print(("PASS " if ok else "FAIL ") + text)
        if not ok:
            fails.append(text)

    p = wedge_case(work)
    sec = vc.run(p)
    a = analyse(p, 180.0 / THETA)          # 2 pi / (2 theta)
    C = beverloo_C(a["W"], 2 * R_OUT)
    print(f"Wedge: {sec:.0f} s; full tank {a['m0']:.3f} kg; W = {a['W']:.4f} kg/s "
          f"over t {a['win'][0]:.2f} to {a['win'][1]:.2f} s; tank at the end "
          f"{a['end']:.3f} of the mass at opening")
    q = quarter_case(work)
    sec = vc.run(q)
    b = analyse(q, 4.0)
    Deq = 2 * (2 * Q_OUT) / math.sqrt(math.pi)    # circle of the square's area
    Cq = beverloo_C(b["W"], Deq)
    Cqh = beverloo_C_square(b["W"], 2 * Q_OUT)   # full square side 2 Q_OUT
    print(f"Quarter: {sec:.0f} s; full tank {b['m0']:.3f} kg; W = {b['W']:.4f} kg/s "
          f"over t {b['win'][0]:.2f} to {b['win'][1]:.2f} s (D_eq {Deq:.4f}); end "
          f"{b['end']:.3f}")

    for name, r in (("wedge", a), ("quarter", b)):
        check(r["leak"] == 0.0, f"{name}: nothing leaves before opening ({r['leak']:.1e} kg)")
        check(r["cons"] < 1e-10, f"{name}: tank plus discharged constant ({r['cons']:.1e})")
    check(C_BAND[0] <= C <= C_BAND[1], f"wedge Beverloo C = {C:.3f} in {C_BAND}")
    check(C_BAND[0] <= Cq <= C_BAND[1], f"quarter Beverloo C = {Cq:.3f} in {C_BAND}")
    gap = abs(C - Cq) / Cq
    check(gap <= 0.10, f"wedge C within 10 % of the quarter's ({100 * gap:+.1f} %)")
    gaph = abs(C - Cqh) / Cqh
    check(gaph <= 0.03, f"wedge C within 3 % of the quarter's by D_h "
          f"(C_q = {Cqh:.3f}, {100 * gaph:+.1f} %)")
    print(f"WP 3.8c checks: {'ALL PASS' if not fails else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
