"""Conical hopper cases on the axisymmetric wedge (WP H.2, decision 0100).

Builds the Couto et al. (2013) test silo (docs/spec/conical-hopper.md, case A) as
a wedge about the y axis: the cylinder wall is the mesh patch `right` at r = 0.5 m,
the conical hopper wall is an STL body (a cut-cell wall, decisions 0047, 0051)
extruded across the wedge in the planes of the wedge's own flat faces, and the
outlet is the mesh floor under r < D0/2.

Used by hopper_checks.py; also runnable by hand:

    python3 hopper_case.py WORKDIR [mc|dp] [cell]
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("void", "silo", "lagrange", "wedge", "wedgevoid", "wedgesilo"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import void_checks as vc  # noqa: E402

THETA = 2.5                      # wedge half angle [deg]

# Case A: Couto silo (P2 Fig. 3) and wheat (P2, P1)
COUTO = dict(R=0.5, Hc=2.0, beta=34.3, Rout=0.175, rho=856.0, E=10.674e6, nu=0.3,
             phi=30.17, psi=1.0, mu=0.20)


def hopper_geometry(c):
    b = math.radians(c["beta"])
    ht = (c["R"] - c["Rout"]) / math.tan(b)
    return b, ht


def material(c, model):
    m = model.lower()
    txt = (f"type {'mohrCoulomb' if m == 'mc' else 'druckerPrager'};\n"
           f"E {c['E']!r};\nnu {c['nu']!r};\nphi {c['phi']!r};\nc 0;\n"
           f"psi {c['psi']!r};\n")
    if m == "mc":
        txt += "sigmaT 0;\n"
    else:
        txt += "fit " + os.environ.get("CEL_DP_FIT", "compression") + ";\n"
    return txt + os.environ.get("CEL_HOPPER_RULES",
                                "shearCap yes;\nspreadLayers yes;\n"
                                "normalCap yes;\n")


def wedge_mesh(path, c, cell, ny):
    """Wedge about y: block 0 from the axis to Rout (its floor is the outlet),
    block 1 from Rout to R."""
    t = math.radians(THETA)
    co, si = math.cos(t), math.sin(t)
    Hm = ny * cell
    R, Ro = c["R"], c["Rout"]

    def P(r, y, sg):
        return (r * co, y, sg * r * si)
    v = [P(0, 0, 0), P(Ro, 0, -1), P(Ro, Hm, -1), P(0, Hm, 0),
         P(Ro, 0, 1), P(Ro, Hm, 1),
         P(R, 0, -1), P(R, Hm, -1), P(R, 0, 1), P(R, Hm, 1)]
    n0, n1 = int(round(Ro / cell)), int(round((R - Ro) / cell))
    assert abs(n0 * cell - Ro) < 1e-9 and abs(n1 * cell - (R - Ro)) < 1e-9, \
        "cell must divide the outlet radius and the silo radius"
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


def hopper_stl(c, ext=6.0, seg=0.01):
    """Closed STL of the solid below and outside the cone: the (r, y) polygon
    (cone extended below the floor, out past R) extruded across the wedge.
    Each vertex (r, y) becomes (r cos t, y, +-ext r sin t): a straight line
    between the two lies in the plane of the wedge's flat cell faces, so the
    cut is the mesh's own conical face, and the body reaches past both wedge
    sides (no surface on a mesh face, decision 0093). ext puts the caps
    outside the mesh's hoop extent and seg splits the edges, so that no
    triangle's box overlaps cells deep inside the body (bodyCells treats a
    cell whose box meets a triangle as a cut candidate)."""
    b, ht = hopper_geometry(c)
    R, Ro = c["R"], c["Rout"]
    yb = -0.05
    poly = [(Ro + yb * math.tan(b), yb), (R + 0.1, ht + 0.1 / math.tan(b)),
            (R + 0.1, yb)]
    corners = list(poly)
    dense = []
    for i in range(len(poly)):
        a, bb = np.array(poly[i]), np.array(poly[(i + 1) % len(poly)])
        k = max(1, int(math.ceil(np.linalg.norm(bb - a) / seg)))
        dense += [tuple(a + (bb - a) * j / k) for j in range(k)]
    poly = dense
    t = math.radians(THETA)
    out = ["solid hopper"]

    def P(rv, sg):
        r, y = rv
        return np.array([r * math.cos(t), y, sg * ext * r * math.sin(t)])

    def tri(p, q, w):
        nn = np.cross(q - p, w - p)
        nn /= np.linalg.norm(nn)
        out.append(f"facet normal {nn[0]:.9e} {nn[1]:.9e} {nn[2]:.9e}\nouter loop")
        for x in (p, q, w):
            out.append(f"vertex {x[0]:.12e} {x[1]:.12e} {x[2]:.12e}")
        out.append("endloop\nendfacet")
    # polygon is counter-clockwise in (r, y); side quads, outward normals
    n = len(poly)
    area = sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1]
               for i in range(n))
    if area < 0:
        poly = poly[::-1]
    for i in range(n):
        a, bb = poly[i], poly[(i + 1) % n]
        tri(P(a, -1), P(bb, -1), P(bb, 1))
        tri(P(a, -1), P(bb, 1), P(a, 1))
    # caps: z > 0 cap normal +z-ish, z < 0 cap normal -z-ish
    # fan from the centroid of the corners (convex polygon, no collinear fans)
    cen = tuple(np.mean(np.array(corners), axis=0))
    for i in range(n):
        tri(P(cen, 1), P(poly[i], 1), P(poly[(i + 1) % n], 1))
        tri(P(cen, -1), P(poly[(i + 1) % n], -1), P(poly[i], -1))
    out.append("endsolid hopper")
    return "\n".join(out) + "\n"


def couto_case(work, name, model="mc", cell=0.025, check=None, c=None, nProcs=1):
    c = dict(COUTO, **(c or {}))
    b, ht = hopper_geometry(c)
    H = ht + c["Hc"]
    ny = int(math.ceil(H / cell - 1e-9)) + 2
    check = dict(check or {})
    check.setdefault("fill", [-1.0, 0.0, 1.0, H])
    check.setdefault("U0", [0.0, 0.0, 0.0])
    p = vc.make_case(f"{work}/{name}", 1, 1, material(c, model), check,
                     nProcs=nProcs)
    wedge_mesh(p, c, cell, ny)
    cs.write(f"{p}/constant/triSurface/hopper.stl", hopper_stl(c))
    dw = math.degrees(math.atan(c["mu"]))
    cs.write(f"{p}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n    hopper\n    {\n        surface hopper.stl;\n"
             f"        friction {dw!r};\n    }}\n}}\n")
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read()
    i = txt.index("walls\n{")
    txt = txt[:i] + ("walls\n{\n" f"    right {{ friction {dw!r}; }}\n"
                     f"    bottom {{ friction {dw!r}; }}\n    top {{ friction 0; }}\n}}\n")
    txt = txt.replace(f"rho {vc.RHO};", f"rho {c['rho']!r};")
    t_open = check.get("tOpen", 1e9)
    txt += f"outlets\n{{\n    outlet {{ friction {dw!r}; openTime {t_open!r}; }}\n}}\n"
    open(mp, "w").write(txt)
    return p, c


def aligned_mesh(path, c, cell, top_rows=0, wall_ratio=1.0, outlet_ratio=1.0):
    """Wall-aligned (body-fitted) wedge mesh (Sia, 2026-11-19): the hopper is
    one trapezoid block whose column lines fan from vertical at the axis to
    the wall slope (rows horizontal); the cylinder block above has square
    cells joined cell to cell. Patches: outlet, hopperWall, right, top,
    front/back (wedge). top_rows adds rows of headroom above the fill.
    Grading (Sia, 2026-11-19): wall_ratio is the wall column's width over
    the axis column's (both blocks, so fine along both walls); outlet_ratio
    is the outlet row's height over the hopper's top row, whose height stays
    about `cell` (the rows above stay uniform)."""
    b, ht = hopper_geometry(c)
    t = math.radians(THETA)
    co, si = math.cos(t), math.sin(t)
    R, Ro = c["R"], c["Rout"]
    n = int(round(R / cell))
    nh = max(1, int(round(ht / cell)))
    if outlet_ratio != 1.0:
        # rows growing from the outlet by a fixed ratio: the smallest count
        # whose top row is at most `cell`
        def top(m):
            q = (1 / outlet_ratio) ** (1 / max(m - 1, 1))
            return ht * (q - 1) / (q ** m - 1) * q ** (m - 1)
        nh = next(m for m in range(2, 400) if top(m) <= cell * 1.0001)
    gy = 1 / outlet_ratio
    nc = int(round(c["Hc"] / cell)) + top_rows
    Hm = ht + nc * c["Hc"] / int(round(c["Hc"] / cell))

    def P(r, y, sg):
        return (r * co, y, sg * r * si)
    v = [P(0, 0, 0), P(Ro, 0, -1), P(R, ht, -1), P(0, ht, 0), P(Ro, 0, 1),
         P(R, ht, 1), P(R, Hm, -1), P(0, Hm, 0), P(R, Hm, 1)]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n"
             f"    hex (0 1 2 3 0 4 5 3) ({n} {nh} 1) simpleGrading ({wall_ratio!r} {gy!r} 1)\n"
             f"    hex (3 2 6 7 3 5 8 7) ({n} {nc} 1) simpleGrading ({wall_ratio!r} 1 1)\n);\n"
             "boundary\n(\n"
             "    outlet { type patch; faces ((0 1 4 0)); }\n"
             "    hopperWall { type patch; faces ((1 2 5 4)); }\n"
             "    right { type patch; faces ((2 6 8 5)); }\n"
             "    top { type patch; faces ((7 8 6 7)); }\n"
             "    front { type wedge; faces ((0 3 2 1) (3 7 6 2)); }\n"
             "    back { type wedge; faces ((0 4 5 3) (3 5 8 7)); }\n"
             "    axis { type empty; faces ((0 3 3 0) (3 7 7 3)); }\n);\n")
    return Hm


def wall_pressures(path, c, cell):
    """Normal wall pressure (kPa) against height: cylinder from the wall column's
    -s_xx; hopper from n.s.n of the material cells whose centre lies within
    1.5 cells of the cone. Returns two arrays of (y, p)."""
    a = np.loadtxt(f"{path}/cells.dat")
    x, y, al, sxx, syy, sxy = a[:, 0], a[:, 1], a[:, 2], a[:, 6], a[:, 7], a[:, 8]
    b, ht = hopper_geometry(c)
    t = math.radians(THETA)
    r = x / math.cos(t)
    col = (r > c["R"] - cell) & (al > 0.5) & (y > ht)
    cyl = np.column_stack([y[col], -sxx[col] / 1e3])
    # distance from the cone line r = Rout + y tan b, positive inside
    nr, ny_ = -math.cos(b), math.sin(b)            # unit normal into the material
    d = (r - (c["Rout"] + y * math.tan(b))) * -math.cos(b)
    sel = (d > 0) & (d < 1.5 * cell) & (al > 0.5) & (y < ht)
    sn = sxx * nr * nr + 2 * sxy * nr * ny_ + syy * ny_ * ny_
    hop = np.column_stack([y[sel], -sn[sel] / 1e3])
    return cyl[np.argsort(cyl[:, 0])], hop[np.argsort(hop[:, 0])]


if __name__ == "__main__":
    work = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else "mc"
    cell = float(sys.argv[3]) if len(sys.argv) > 3 else 0.025
    p, c = couto_case(work, f"couto_{model}", model, cell,
                      dict(tEnd=0.5, maxCo=0.4))
    open(f"{p}/system/controlDict", "a").write("gravityRamp 0.1;\ndamping 100;\n")
    print("case", p)
