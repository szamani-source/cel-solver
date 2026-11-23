"""Couto silo as an STL shell inside a larger Eulerian mesh (Sia, 2026-11-20).

The silo wall is a closed STL body, a shell of thickness TS around the
material, with rounded corners (fillets at the cylinder-hopper transition
and at the outlet edge). The graded Cartesian wedge mesh reaches 2 to 3
cells past the shell's outer face, so any material that leaks through the
wall shows up in that buffer, and the leak monitor reports its mass and
where it is. The outlet is the mesh patch under the opening, as before.

    python3 buffered_case.py WORKDIR [picture.png]   # mesh and picture only
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hopper_case as hc  # noqa: E402

cs, vc = hc.cs, hc.vc

TS = 0.02           # shell thickness [m]
R_TRANS = 0.05      # fillet radius at the cylinder-hopper transition [m]
R_OUT = 0.015       # fillet radius at the outlet edge [m]
BUFFER = 3          # cells of buffer outside the shell


def inner_profile(c):
    """(r, y) of the silo's inner wall from just below the outlet to the top:
    a vertical lip r = Rout below y = 0 rounded into the cone (radius R_OUT;
    the arc meets the lip 4.6 mm below the mesh, so the opening at y = 0
    stays within 1 mm of Rout), the cone, and the cone rounded into the
    cylinder (radius R_TRANS)."""
    b, ht = hc.hopper_geometry(c)
    R, Ro, Htop = c["R"], c["Rout"], ht + c["Hc"]
    tb = math.tan(b)
    s = math.hypot(1.0, tb)
    nrm = np.array([1.0, -tb]) / s          # cone normal out of the material
    # outlet: centre outside the material, tangent to r = Ro and the cone
    yv = -R_OUT * (s - 1) / tb
    c1 = np.array([Ro + R_OUT, yv])
    a1 = math.atan2(tb, -1.0)               # direction of -nrm
    # transition: centre inside, tangent to the cone and to r = R
    yt = (R - R_TRANS - Ro + R_TRANS * s) / tb
    c2 = np.array([R - R_TRANS, yt])
    a2 = math.atan2(-tb, 1.0)               # direction of nrm
    pts = [c1 + R_OUT * np.array([math.cos(a), math.sin(a)])
           for a in np.linspace(math.pi, a1, 31)]
    pts += [c2 + R_TRANS * np.array([math.cos(a), math.sin(a)])
            for a in np.linspace(a2, 0.0, 41)]
    pts.append(np.array([R, Htop]))
    return np.array(pts)


def shell_polygon(c, yb=-0.05, y_top=None):
    """Closed (r, y) polygon of the shell: up the inner wall from below the
    mesh, across the top, down the outer face (the cone offset by TS, the
    cylinder at R + TS) to just below the mesh, and back under the lip."""
    b, ht = hc.hopper_geometry(c)
    R, Ro, Htop = c["R"], c["Rout"], ht + c["Hc"]
    tb = math.tan(b)
    inner = inner_profile(c)
    if y_top is not None:
        Htop = y_top
        inner[-1, 1] = Htop
    off = TS * math.hypot(1.0, tb)          # horizontal offset of the outer cone
    y_join = (R + TS - Ro - off) / tb       # outer cone meets r = R + TS
    y_low = -0.01
    r_low = Ro + off + y_low * tb
    poly = [(Ro, yb)] + [tuple(p) for p in inner] + \
        [(R + TS, Htop), (R + TS, y_join), (r_low, y_low), (r_low, yb)]
    return np.array(poly)


def graded(L, segs):
    """blockMesh multi-grading for one direction: segs = [(length, n, ratio)]."""
    tot = sum(s[0] for s in segs)
    ntot = sum(s[1] for s in segs)
    inner = " ".join(f"({s[0] / tot!r} {s[1] / ntot!r} {s[2]!r})" for s in segs)
    return ntot, f"({inner})"


def seg(length, h0, h1):
    """(length, n, ratio) for cells growing from h0 to h1 over length."""
    if abs(h1 - h0) < 1e-12:
        n = max(1, int(round(length / h0)))
        return (length, n, 1.0)
    n = max(1, int(round(length / (0.5 * (h0 + h1)))))
    return (length, n, h1 / h0)


Y0 = -0.02         # mesh floor: a 20 mm spout below the hopper, so the
                   # rounded outlet edge lies wholly inside the mesh


def buffered_mesh(path, c, h_out=0.004, h_wall=0.007, h_mid=0.012, h_top=0.025,
                  y_out=0.25, y_hop=0.6):
    """Graded wedge mesh about y. Radial: Rout block at h_out; then growing to
    h_mid and back to h_wall at the cylinder wall, h_wall through the shell
    and BUFFER cells past it. Vertical from the floor Y0 (outlet patch under
    the spout): h_out growing to h_wall by y_out, h_wall up to y_hop (past
    the transition fillet), then growing to h_top over 0.3 m; three top
    rows of headroom."""
    b, ht = hc.hopper_geometry(c)
    R, Ro, Htop = c["R"], c["Rout"], ht + c["Hc"]
    Rd = R + TS + BUFFER * h_wall
    t = math.radians(hc.THETA)
    co, si = math.cos(t), math.sin(t)
    n0 = int(round(Ro / h_out))
    rm = 0.5 * (Ro + R)
    nr, gr = graded(Rd - Ro, [seg(rm - Ro, h_out, h_mid), seg(R - 0.03 - rm, h_mid, h_wall),
                               seg(Rd - (R - 0.03), h_wall, h_wall)])
    Hm = Htop + 3 * h_top
    ny, gy = graded(Hm - Y0, [seg(y_out - Y0, h_out, h_wall),
                              seg(y_hop - y_out, h_wall, h_wall),
                              seg(0.3, h_wall, h_top), seg(Hm - y_hop - 0.3, h_top, h_top)])

    def P(r, y, sg):
        return (r * co, y, sg * r * si)
    y0 = Y0
    v = [P(0, y0, 0), P(Ro, y0, -1), P(Ro, Hm, -1), P(0, Hm, 0),
         P(Ro, y0, 1), P(Ro, Hm, 1),
         P(Rd, y0, -1), P(Rd, Hm, -1), P(Rd, y0, 1), P(Rd, Hm, 1)]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    blocks = (f"    hex (0 1 2 3 0 4 5 3) ({n0} {ny} 1) simpleGrading (1 {gy} 1)\n"
              f"    hex (1 6 7 2 4 8 9 5) ({nr} {ny} 1) simpleGrading ({gr} {gy} 1)\n")
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
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n{blocks});\n"
             "boundary\n(\n" + "\n".join(pt) + "\n);\n")
    return Rd, Hm


def shell_stl(c, ext=6.0, seg_len=0.005, y_top=None):
    """The shell polygon extruded across the wedge (as hc.hopper_stl)."""
    poly = shell_polygon(c, y_top=y_top)
    dense = []
    n = len(poly)
    for i in range(n):
        a, bb = poly[i], poly[(i + 1) % n]
        k = max(1, int(math.ceil(np.linalg.norm(bb - a) / seg_len)))
        dense += [tuple(a + (bb - a) * j / k) for j in range(k)]
    poly = dense
    t = math.radians(hc.THETA)
    out = ["solid shell"]

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
    n = len(poly)
    area = sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1]
               for i in range(n))
    if area < 0:
        poly = poly[::-1]
    for i in range(n):
        a, bb = poly[i], poly[(i + 1) % n]
        tri(P(a, -1), P(bb, -1), P(bb, 1))
        tri(P(a, -1), P(bb, 1), P(a, 1))
    # caps: the shell is not convex, so fan each cap by ear clipping
    idx = list(range(n))
    tris = []
    guard = 0
    while len(idx) > 3 and guard < 100000:
        guard += 1
        for k in range(len(idx)):
            i0, i1, i2 = idx[k - 1], idx[k], idx[(k + 1) % len(idx)]
            p0, p1, p2 = (np.array(poly[i]) for i in (i0, i1, i2))
            cr = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])
            if cr <= 1e-14:
                continue
            ok = True
            for j in idx:
                if j in (i0, i1, i2):
                    continue
                q = np.array(poly[j])
                d = [((b2[0] - a2[0]) * (q[1] - a2[1]) - (b2[1] - a2[1]) * (q[0] - a2[0]))
                     for a2, b2 in ((p0, p1), (p1, p2), (p2, p0))]
                if min(d) >= -1e-14:
                    ok = False
                    break
            if ok:
                tris.append((i0, i1, i2))
                idx.pop(k)
                break
    tris.append(tuple(idx))
    for i0, i1, i2 in tris:
        tri(P(poly[i0], 1), P(poly[i1], 1), P(poly[i2], 1))
        tri(P(poly[i0], -1), P(poly[i2], -1), P(poly[i1], -1))
    out.append("endsolid shell")
    return "\n".join(out) + "\n"


def buffered_case(work, name, model="mc", check=None, c=None, nProcs=1, **mesh):
    c = dict(hc.COUTO, **(c or {}))
    b, ht = hc.hopper_geometry(c)
    H = ht + c["Hc"]
    check = dict(check or {})
    check.setdefault("fill", [-1.0, Y0, 1.0, H])
    check.setdefault("U0", [0.0, 0.0, 0.0])
    # material only inside the wall: left of the middle of the shell
    prof = inner_profile(c)
    ct = math.cos(math.radians(hc.THETA))
    half = 0.4 * TS * math.hypot(1.0, math.tan(hc.hopper_geometry(c)[0]))
    check.setdefault("fillRight", [[Y0 - 1.0, (c["Rout"] + 0.5 * TS) * ct]]
                     + [[float(y), float((r + half) * ct)] for r, y in prof]
                     + [[H + 1.0, (c["R"] + 0.5 * TS) * ct]])
    p = vc.make_case(f"{work}/{name}", 1, 1, hc.material(c, model), check, nProcs=nProcs)
    buffered_mesh(p, c, **mesh)
    cs.write(f"{p}/constant/triSurface/shell.stl", shell_stl(c))
    dw = math.degrees(math.atan(c["mu"]))
    cs.write(f"{p}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n    shell\n    {\n        surface shell.stl;\n"
             f"        friction {dw!r};\n    }}\n}}\n")
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read()
    i = txt.index("walls\n{")
    txt = txt[:i] + ("walls\n{\n    right { friction 0; }\n    bottom { friction 0; }\n"
                     "    top { friction 0; }\n}\n")
    txt = txt.replace(f"rho {vc.RHO};", f"rho {c['rho']!r};")
    t_open = check.get("tOpen", 1e9)
    txt += f"outlets\n{{\n    outlet {{ friction {dw!r}; openTime {t_open!r}; }}\n}}\n"
    open(mp, "w").write(txt)
    return p, c


def inside_silo(c, r, y):
    """True where (r, y) lies inside the silo's inner wall."""
    prof = inner_profile(c)
    rin = np.interp(y, prof[:, 1], prof[:, 0])
    return (r < rin) & (y >= 0)


def leak_report(c, cells):
    """Mass per radian... returned as (wedge mass outside the shell, its
    lowest and highest y, its r range) from a cells array
    (x y alpha rho ... open), wedge volume weights from the cell size."""
    t = math.radians(hc.THETA)
    r = cells[:, 0] / math.cos(t)
    y = cells[:, 1]
    prof = inner_profile(c)
    rin = np.interp(y, prof[:, 1], prof[:, 0], left=c["Rout"], right=c["R"])
    tb = math.tan(hc.hopper_geometry(c)[0])
    out = (r > rin + TS * np.where(y < hc.hopper_geometry(c)[1], math.hypot(1.0, tb), 1.0)) & (cells[:, 2] > 1e-6)
    return out, r, y


def picture(p, c, png, heading=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    sys.path.insert(0, os.path.join(HERE, "..", "angled"))
    import angled_checks as ac
    m = f"{p}/constant/polyMesh"
    pts = ac.foam_list(f"{m}/points", "point")
    faces = ac.foam_list(f"{m}/faces", "face")
    t = math.radians(hc.THETA)
    segs = []
    for f in faces:
        q = pts[f]
        if np.all(q[:, 2] <= 1e-12):          # front wedge face (z <= 0 side)
            rr = q[:, 0] / math.cos(t)
            for i in range(len(f)):
                j = (i + 1) % len(f)
                segs.append([(rr[i], q[i, 1]), (rr[j], q[j, 1])])
    poly = shell_polygon(c)
    b, ht = hc.hopper_geometry(c)
    fig, axs = plt.subplots(1, 3, figsize=(15, 9), gridspec_kw=dict(width_ratios=[1.1, 1, 1]))
    views = [((-0.02, 0.6), (-0.05, 2.65), "Whole silo"),
             ((0.10, 0.32), (-0.03, 0.16), "Outlet edge (fillet 15 mm, 20 mm spout)"),
             ((0.36, 0.60), (0.36, 0.62), "Transition (fillet 50 mm)")]
    for ax, (xl, yl, title) in zip(axs, views):
        ax.add_collection(LineCollection(segs, colors="#9aa4b1", linewidths=0.35))
        pp = np.vstack([poly, poly[:1]])
        ax.fill(pp[:, 0], pp[:, 1], color="#c0392b", alpha=0.35, lw=0)
        ax.plot(pp[:, 0], pp[:, 1], color="#c0392b", lw=1.2, label="STL shell (20 mm)")
        prof = np.vstack([[c["Rout"], Y0], inner_profile(c)])
        ax.fill_betweenx(prof[:, 1], 0, prof[:, 0], color="#e6c27a", alpha=0.35,
                         label="material (filled to the top)")
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("r (m)")
    axs[0].set_ylabel("y (m)")
    axs[0].legend(loc="upper left", fontsize=8)
    ncell = ac.foam_list(f'{m}/owner', 'label').max() + 1
    fig.suptitle(heading.format(ncell=ncell) if heading else
                 f"Couto silo, STL shell inside a graded wedge mesh ({ncell} cells; "
                 f"outlet 4 mm, walls 7 mm, {BUFFER} buffer cells)")
    fig.tight_layout()
    fig.savefig(png, dpi=130)


if __name__ == "__main__":
    work = sys.argv[1]
    p, c = buffered_case(work, "couto_mc_buffered", "mc")
    import subprocess
    subprocess.run(["blockMesh"], cwd=p, stdout=open(f"{p}/log.blockMesh", "w"),
                   stderr=subprocess.STDOUT, check=True)
    subprocess.run(["checkMesh"], cwd=p, stdout=open(f"{p}/log.checkMesh", "w"),
                   stderr=subprocess.STDOUT)
    print(open(f"{p}/log.checkMesh").read().split("Checking geometry")[0][-600:])
    if len(sys.argv) > 2:
        picture(p, c, sys.argv[2])


# --- Wall-aligned version (Sia, 2026-11-20: "We need the aligned mesh") ----

def offset_path(pts, d):
    """Offset an (r, y) path (going up, material on its left = smaller r)
    by d to its right (out of the material), along the vertex normals."""
    pts = np.asarray(pts, float)
    tng = np.gradient(pts, axis=0)
    tng /= np.linalg.norm(tng, axis=1)[:, None]
    nrm = np.column_stack([tng[:, 1], -tng[:, 0]])      # right of an upward path
    return pts + d * nrm


def wall_path(c, y_top, n_extra=0):
    """Inner wall from the mesh floor Y0 up to y_top: the lip r = Rout,
    the outlet fillet, the cone, the transition fillet, the cylinder."""
    prof = inner_profile(c)
    pts = [(c["Rout"], Y0)] + [tuple(p) for p in prof[:-1] if p[1] < y_top]
    pts.append((c["R"], y_top))
    # densify straight runs so the offset follows them
    dense = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = np.array(a), np.array(b)
        k = max(1, int(math.ceil(np.linalg.norm(b - a) / 0.005)))
        dense += [a + (b - a) * j / k for j in range(k)]
    dense.append(np.array(pts[-1]))
    return np.array(dense)


def outer_path(c, D, y_top):
    """The domain's outer side, D outside the inner wall: the lip and the
    cone offsets meeting at a corner, the transition fillet (radius
    R_TRANS + D), the cylinder r = R + D."""
    b, ht = hc.hopper_geometry(c)
    R, Ro = c["R"], c["Rout"]
    tb = math.tan(b)
    s = math.hypot(1.0, tb)
    y_c = D * (1 - s) / tb                       # lip offset meets cone offset
    pts = [(Ro + D, Y0)]
    if y_c > Y0:
        pts.append((Ro + D, y_c))
    yt = (R - R_TRANS - Ro + R_TRANS * s) / tb
    c2 = np.array([R - R_TRANS, yt])
    nrm = np.array([1.0, -tb]) / s
    for a in np.linspace(math.atan2(-tb, 1.0), 0.0, 41):
        pts.append(tuple(c2 + (R_TRANS + D) * np.array([math.cos(a), math.sin(a)])))
    pts.append((R + D, y_top))
    dense = []
    for a, b2 in zip(pts[:-1], pts[1:]):
        a, b2 = np.array(a), np.array(b2)
        k = max(1, int(math.ceil(np.linalg.norm(b2 - a) / 0.005)))
        dense += [a + (b2 - a) * j / k for j in range(k)]
    dense.append(np.array(pts[-1]))
    return np.array(dense)


def aligned_buffered_mesh(path, c, ncols=44, wall_ratio=0.6, h_out=0.004, h_hop=0.012,
                          h_top=0.025, y_t=0.55, n_buf=4, head_rows=3,
                          outlet_wall_ratio=None):
    """Four curved blocks in the wedge: inside the shell's mid-line M (the
    inner wall offset by TS/2) the hopper block (axis to M, its right edge
    follows the lip, both fillets and the cone, so the columns lean with the
    wall) and the cylinder block; outside M, n_buf columns out to the
    domain side O (TS + 3 x 7 mm outside the inner wall): the shell's outer
    half and about 3 cells of buffer, aligned with the wall. Rows: h_out at
    the outlet growing to h_hop at y_t, then to h_top in the cylinder;
    head_rows rows of h_top above the silo top (room for a heap when the
    silo is filled through an inlet). outlet_wall_ratio (default
    wall_ratio): the column grading (wall column / axis column) along the
    hopper block's bottom edge, the outlet; the top edge keeps wall_ratio,
    so the columns match the cylinder block and blockMesh blends the two
    down the cone (edgeGrading). 1 keeps the outlet cells at the block
    width over ncols (decision 0117: >= 3 mm with 58 columns)."""
    b, ht = hc.hopper_geometry(c)
    R, Htop = c["R"], ht + c["Hc"]
    Hm = Htop + head_rows * h_top
    D = TS + 3 * 0.007
    t = math.radians(hc.THETA)
    co, si = math.cos(t), math.sin(t)
    W = wall_path(c, y_t)
    M = offset_path(W, 0.5 * TS)
    M[0, 1], M[-1, 1] = Y0, y_t
    O = outer_path(c, D, y_t)

    def P(r, y, sg):
        r, y = float(r), float(y)
        return (r * co, y, sg * r * si)
    a = [(0.0, Y0), (0.0, y_t), (0.0, Hm)]
    m = [tuple(M[0]), tuple(M[-1]), (R + 0.5 * TS, Hm)]
    o = [tuple(O[0]), tuple(O[-1]), (R + D, Hm)]
    v = [P(r, y, 0) for r, y in a]
    for sg in (-1, 1):
        v += [P(r, y, sg) for r, y in m] + [P(r, y, sg) for r, y in o]
    # 0-2 axis; front 3-5 m, 6-8 o; back 9-11 m, 12-14 o
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    n1, g1 = graded(y_t - Y0, [seg(y_t - Y0, h_out, h_hop)])
    n2, g2 = graded(Hm - y_t, [seg(0.5, h_hop, h_top), seg(Hm - y_t - 0.5, h_top, h_top)])
    wr = f"{wall_ratio!r}"
    owr = f"{(wall_ratio if outlet_wall_ratio is None else outlet_wall_ratio)!r}"
    # x edges of hex (0 3 4 1 0 9 10 1): bottom front, top front, top back, bottom back
    hop = (f"edgeGrading ({owr} {wr} {wr} {owr} {g1} {g1} {g1} {g1} 1 1 1 1)"
           if owr != wr else f"simpleGrading ({wr} {g1} 1)")
    blocks = (f"    hex (0 3 4 1 0 9 10 1) ({ncols} {n1} 1) {hop}\n"
              f"    hex (1 4 5 2 1 10 11 2) ({ncols} {n2} 1) simpleGrading ({wr} {g2} 1)\n"
              f"    hex (3 6 7 4 9 12 13 10) ({n_buf} {n1} 1) simpleGrading (1 {g1} 1)\n"
              f"    hex (4 7 8 5 10 13 14 11) ({n_buf} {n2} 1) simpleGrading (1 {g2} 1)\n")

    def pl(i, j, pts, sg):
        inner = " ".join("({!r} {!r} {!r})".format(*P(r, y, sg)) for r, y in pts[1:-1])
        return f"    polyLine {i} {j} ({inner})\n"
    edges = (pl(3, 4, M, -1) + pl(9, 10, M, 1) + pl(6, 7, O, -1) + pl(12, 13, O, 1))
    patches = {
        "outlet": ("patch", "(0 3 9 0)"),
        "bottom": ("patch", "(3 6 12 9)"),
        "outer": ("patch", "(6 7 13 12) (7 8 14 13)"),
        "top": ("patch", "(2 5 11 2) (5 8 14 11)"),
        "front": ("wedge", "(0 1 4 3) (1 2 5 4) (3 4 7 6) (4 5 8 7)"),
        "back": ("wedge", "(0 9 10 1) (1 10 11 2) (9 12 13 10) (10 13 14 11)"),
    }
    pt = [f"    {k}\n    {{\n        type {ty};\n        faces ({f});\n    }}"
          for k, (ty, f) in patches.items()]
    pt.append("    axis\n    {\n        type empty;\n        faces ((0 1 1 0) (1 2 2 1));\n    }")
    cs.write(f"{path}/system/blockMeshDict",
             cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n{blocks});\n"
             f"edges\n(\n{edges});\n"
             "boundary\n(\n" + "\n".join(pt) + "\n);\n")
    return Hm


def aligned_buffered_case(work, name, model="mc", check=None, c=None, nProcs=1, **mesh):
    c = dict(hc.COUTO, **(c or {}))
    b, ht = hc.hopper_geometry(c)
    H = ht + c["Hc"]
    check = dict(check or {})
    check.setdefault("fill", [-1.0, Y0, 1.0, H])
    check.setdefault("U0", [0.0, 0.0, 0.0])
    prof = inner_profile(c)
    ct = math.cos(math.radians(hc.THETA))
    half = 0.4 * TS * math.hypot(1.0, math.tan(b))
    check.setdefault("fillRight", [[Y0 - 1.0, (c["Rout"] + 0.4 * TS) * ct]]
                     + [[float(y), float((r + half) * ct)] for r, y in prof]
                     + [[H + 1.0, (c["R"] + 0.4 * TS) * ct]])
    p = vc.make_case(f"{work}/{name}", 1, 1, hc.material(c, model), check, nProcs=nProcs)
    Hm = aligned_buffered_mesh(p, c, **mesh)
    cs.write(f"{p}/constant/triSurface/shell.stl", shell_stl(c, y_top=Hm + 0.05))
    dw = math.degrees(math.atan(c["mu"]))
    cs.write(f"{p}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n    shell\n    {\n        surface shell.stl;\n"
             f"        friction {dw!r};\n    }}\n}}\n")
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read()
    i = txt.index("walls\n{")
    txt = txt[:i] + ("walls\n{\n    outer { friction 0; }\n    bottom { friction 0; }\n"
                     "    top { friction 0; }\n}\n")
    txt = txt.replace(f"rho {vc.RHO};", f"rho {c['rho']!r};")
    t_open = check.get("tOpen", 1e9)
    txt += f"outlets\n{{\n    outlet {{ friction {dw!r}; openTime {t_open!r}; }}\n}}\n"
    open(mp, "w").write(txt)
    return p, c
