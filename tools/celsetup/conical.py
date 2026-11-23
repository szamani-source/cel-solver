"""Template conicalHopper: a cylinder on a conical hopper, on an
axisymmetric wedge, with the setup that worked best for the Couto silo
(decisions 0105, 0107, 0115; manual, Meshing):

- a wall-aligned blockMesh: four curved blocks whose column lines follow the
  wall (lip, outlet fillet, cone, transition fillet, cylinder), rows
  horizontal, fine rows at the outlet growing upwards;
- the silo wall as a closed STL shell of thickness `walls.thickness` inside
  the mesh, with rounded corners and a short spout below the outlet, and
  about three buffer cells outside it (a leak shows up there);
- the outlet is the mesh patch under the spout.

The geometry code is the one of tests/hopper/buffered_case.py
(aligned_buffered_mesh), with the Couto numbers turned into inputs. With
the Couto inputs and the same mesh numbers it writes the same mesh
(tests/test_celsetup.sh checks this).
"""

import math

import numpy as np

from . import foam, schema
from .schema import E

R_ = schema.REQUIRED
SCHEMA = {
    "geometry": [
        E("diameter", R_, "m", "inner diameter of the cylinder"),
        E("cylinderHeight", R_, "m", "height of the cylinder above the hopper"),
        E("hopperAngle", R_, "deg", "hopper wall angle from the vertical (0 to 89)"),
        E("outletDiameter", R_, "m", "outlet diameter"),
        E("fillHeight", "top", "m", "material fill height above the outlet, or top"),
        E("wedgeAngle", 5.0, "deg", "whole angle of the axisymmetric wedge (5 is tested)",
          level="advanced"),
    ],
    "walls": [
        E("frictionCoefficient", None, "", "wall friction coefficient mu_w (or give friction)"),
        E("friction", None, "deg", "wall friction angle, instead of frictionCoefficient"),
        E("thickness", 0.02, "m", "thickness of the STL shell that forms the wall"),
        E("filletTransition", 0.05, "m", "radius of the rounded cylinder-hopper corner "
          "(rounded corners raised the rate without changing wall pressures)"),
        E("filletOutlet", 0.015, "m", "radius of the rounded outlet edge"),
        E("spout", 0.02, "m", "short vertical lip below the outlet; the outlet patch is its end"),
    ],
    "mesh": [
        E("preset", "standard", "", "standard: the validated 7,584-cell Couto mesh numbers; fine: "
          "cells near the nonlocal length (decision 0115); entries below set to a value override it",
          ["standard", "fine"]),
        E("mesher", "blockMesh", "", "blockMesh (OpenFOAM) or gmsh (Gmsh, GPL, transfinite "
          "quads: the same block layout and cell counts, interior nodes placed a little "
          "differently)", ["blockMesh", "gmsh"]),
        E("outletCell", "preset", "m", "row height at the outlet: 3 to 5 mm sets the discharge rate; preset: from the mesh preset"),
        E("hopperTopCell", "preset", "m", "row height at the top of the hopper block; preset: from the mesh preset"),
        E("topCell", "preset", "m", "row height in the upper cylinder; preset: from the mesh preset"),
        E("columns", "preset", "", "columns from the axis to the wall; auto: outlet radius / outletCell; preset: from the mesh preset"),
        E("wallRatio", "preset", "", "wall column width / axis column width (smaller: finer at the wall); preset: from the mesh preset"),
        E("growthLength", 0.5, "m", "rows grow from hopperTopCell to topCell over this height",
          level="advanced"),
        E("bufferCell", 0.007, "m", "size of the buffer cells outside the shell", level="advanced"),
        E("bufferCells", 3, "", "buffer width outside the shell, in bufferCell (leak monitor)",
          level="advanced"),
        E("bufferColumns", 4, "", "columns from the shell middle to the domain side",
          level="advanced"),
        E("headroom", None, "m", "empty space above the fill; default 3 top rows",
          level="advanced"),
        E("transitionMargin", 0.05, "m", "the hopper block ends this far above the transition "
          "fillet", level="advanced"),
    ],
}

# Mesh presets (overridden by any entry given in the YAML file)
PRESETS = {
    # the 7,584-cell Couto mesh: MC 124 kg/s, Beverloo 110-131 (decision 0107)
    "standard": dict(outletCell=0.004, hopperTopCell=0.012, topCell=0.025, columns="auto",
                     wallRatio=0.6),
    # decision 0115 addendum: cells near the nonlocal length (about 12 mm)
    "fine": dict(outletCell=0.004, hopperTopCell=0.004, topCell=0.012, columns=58,
                 wallRatio=0.25),
}


class Geometry:
    """Profile of the inner wall, the shell polygon and the domain side, in
    (r, y) with y = 0 at the outlet and the mesh floor at -spout."""

    def __init__(self, g, w):
        self.R = 0.5 * float(g["diameter"])
        self.Ro = 0.5 * float(g["outletDiameter"])
        self.Hc = float(g["cylinderHeight"])
        self.beta = math.radians(float(g["hopperAngle"]))
        self.TS = float(w["thickness"])
        self.R_TRANS = float(w["filletTransition"])
        self.R_OUT = float(w["filletOutlet"])
        self.Y0 = -float(w["spout"])
        self.theta = math.radians(0.5 * float(g["wedgeAngle"]))
        if not 0 < self.Ro < self.R:
            raise ValueError("geometry: need 0 < outletDiameter < diameter")
        if not 0 < self.beta < math.radians(89):
            raise ValueError("geometry: hopperAngle is measured from the vertical, "
                             "0 to 89 deg (use the flatBottomSilo template for 90)")
        self.ht = (self.R - self.Ro) / math.tan(self.beta)
        self.H = self.ht + self.Hc
        tb = math.tan(self.beta)
        s = math.hypot(1.0, tb)
        # transition fillet centre height and its top (tangent to r = R)
        self.yt_fillet = (self.R - self.R_TRANS - self.Ro + self.R_TRANS * s) / tb
        if self.R_TRANS > 0 and self.yt_fillet - self.R_TRANS * (1 - 1 / s) < 0:
            raise ValueError("walls: filletTransition too large for this hopper")

    def inner_profile(self):
        """(r, y) of the inner wall from just below the outlet to the top:
        the lip r = Ro rounded into the cone (radius R_OUT), the cone, and
        the cone rounded into the cylinder (radius R_TRANS)."""
        R, Ro, Htop = self.R, self.Ro, self.H
        tb = math.tan(self.beta)
        s = math.hypot(1.0, tb)
        yv = -self.R_OUT * (s - 1) / tb
        c1 = np.array([Ro + self.R_OUT, yv])
        a1 = math.atan2(tb, -1.0)
        c2 = np.array([R - self.R_TRANS, self.yt_fillet])
        a2 = math.atan2(-tb, 1.0)
        pts = [c1 + self.R_OUT * np.array([math.cos(a), math.sin(a)])
               for a in np.linspace(math.pi, a1, 31)]
        pts += [c2 + self.R_TRANS * np.array([math.cos(a), math.sin(a)])
                for a in np.linspace(a2, 0.0, 41)]
        pts.append(np.array([R, Htop]))
        return np.array(pts)

    def shell_polygon(self, y_top=None, yb=-0.05):
        R, Ro = self.R, self.Ro
        Htop = self.H if y_top is None else y_top
        tb = math.tan(self.beta)
        inner = self.inner_profile()
        inner[-1, 1] = Htop
        off = self.TS * math.hypot(1.0, tb)
        y_join = (R + self.TS - Ro - off) / tb
        y_low = -0.01
        r_low = Ro + off + y_low * tb
        yb = min(yb, self.Y0 - 0.03)
        poly = [(Ro, yb)] + [tuple(p) for p in inner] + \
            [(R + self.TS, Htop), (R + self.TS, y_join), (r_low, y_low), (r_low, yb)]
        return np.array(poly)

    def wall_path(self, y_top):
        prof = self.inner_profile()
        pts = [(self.Ro, self.Y0)] + [tuple(p) for p in prof[:-1] if p[1] < y_top]
        pts.append((self.R, y_top))
        return _densify(pts)

    def outer_path(self, D, y_top):
        R, Ro = self.R, self.Ro
        tb = math.tan(self.beta)
        s = math.hypot(1.0, tb)
        y_c = D * (1 - s) / tb
        pts = [(Ro + D, self.Y0)]
        if y_c > self.Y0:
            pts.append((Ro + D, y_c))
        c2 = np.array([R - self.R_TRANS, self.yt_fillet])
        for a in np.linspace(math.atan2(-tb, 1.0), 0.0, 41):
            pts.append(tuple(c2 + (self.R_TRANS + D) * np.array([math.cos(a), math.sin(a)])))
        pts.append((R + D, y_top))
        return _densify(pts)

    def fill_polygon(self, y_fill):
        """Material region at the start: inside the inner wall, from the
        closed outlet at the mesh floor up to y_fill."""
        prof = self.inner_profile()
        r_at = float(np.interp(y_fill, prof[:, 1], prof[:, 0]))
        pts = [(0.0, self.Y0), (self.Ro, self.Y0)]
        pts += [tuple(p) for p in prof if p[1] < y_fill]
        pts += [(r_at, y_fill), (0.0, y_fill)]
        return np.array(pts)

    def volume_below(self, y):
        """Volume inside the inner wall from the mesh floor up to y (m^3)."""
        prof = np.vstack([[self.Ro, self.Y0], self.inner_profile()])
        ys = np.linspace(self.Y0, y, 4001)
        r = np.interp(ys, prof[:, 1], prof[:, 0])
        f = math.pi * r ** 2
        return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(ys)))


def _densify(pts, step=0.005):
    dense = []
    for a, b in zip(pts[:-1], pts[1:]):
        a, b = np.array(a, float), np.array(b, float)
        k = max(1, int(math.ceil(np.linalg.norm(b - a) / step)))
        dense += [a + (b - a) * j / k for j in range(k)]
    dense.append(np.array(pts[-1], float))
    return np.array(dense)


def offset_path(pts, d):
    """Offset an upward (r, y) path by d to its right (out of the material)."""
    pts = np.asarray(pts, float)
    tng = np.gradient(pts, axis=0)
    tng /= np.linalg.norm(tng, axis=1)[:, None]
    nrm = np.column_stack([tng[:, 1], -tng[:, 0]])
    return pts + d * nrm


def graded(segs):
    """blockMesh multi-grading for one direction: segs = [(length, n, ratio)]."""
    tot = sum(s[0] for s in segs)
    ntot = sum(s[1] for s in segs)
    inner = " ".join(f"({s[0] / tot!r} {s[1] / ntot!r} {s[2]!r})" for s in segs)
    return ntot, f"({inner})"


def seg(length, h0, h1):
    """(length, n, ratio) for cells growing from h0 to h1 over length."""
    if abs(h1 - h0) < 1e-12:
        return (length, max(1, int(round(length / h0))), 1.0)
    return (length, max(1, int(round(length / (0.5 * (h0 + h1))))), h1 / h0)


def layout(G, m, y_fill):
    """The block layout shared by the blockMesh and Gmsh meshers: block
    corners, the shell mid-line M and the domain side O (r, y paths), column
    and row counts with their grading segments [(length, n, ratio)]."""
    R = G.R
    h_out, h_hop, h_top = m["outletCell"], m["hopperTopCell"], m["topCell"]
    headroom = m["headroom"] if m["headroom"] is not None else 3 * h_top
    Hm = y_fill + headroom
    D = G.TS + m["bufferCells"] * m["bufferCell"]
    y_t = math.ceil((G.yt_fillet + m["transitionMargin"]) / 0.05 - 1e-9) * 0.05
    if y_t >= Hm - h_top:
        raise ValueError("mesh: the fill top is below the end of the hopper block")
    ncols = m["columns"]
    if ncols in (None, "auto"):
        ncols = int(round((G.Ro + 0.5 * G.TS) / h_out))
    ncols = int(ncols)
    W = G.wall_path(y_t)
    M = offset_path(W, 0.5 * G.TS)
    M[0, 1], M[-1, 1] = G.Y0, y_t
    O = G.outer_path(D, y_t)
    rows1 = [seg(y_t - G.Y0, h_out, h_hop)]
    gl = m["growthLength"]
    if Hm - y_t > gl + 1e-9:
        rows2 = [seg(gl, h_hop, h_top), seg(Hm - y_t - gl, h_top, h_top)]
    else:
        rows2 = [seg(Hm - y_t, h_hop, h_top)]
    return dict(R=R, Hm=Hm, D=D, y_t=y_t, ncols=ncols, nb=int(m["bufferColumns"]), M=M, O=O,
                rows1=rows1, rows2=rows2, wallRatio=float(m["wallRatio"]))


def block_mesh_dict(G, m, y_fill, src, depth=None):
    """The wall-aligned four-block wedge mesh (aligned_buffered_mesh).
    depth: a plane-strain mesh of that thickness instead (x = r, symmetry
    plane at x = 0)."""
    L = layout(G, m, y_fill)
    R, Hm, D, y_t, ncols, M, O = L["R"], L["Hm"], L["D"], L["y_t"], L["ncols"], L["M"], L["O"]
    h_out, h_hop, h_top = m["outletCell"], m["hopperTopCell"], m["topCell"]
    co, si = math.cos(G.theta), math.sin(G.theta)

    def P(r, y, sg):
        r, y = float(r), float(y)
        if depth:
            return (r, y, sg * 0.5 * depth)
        return (r * co, y, sg * r * si)
    a = [(0.0, G.Y0), (0.0, y_t), (0.0, Hm)]
    mm = [tuple(M[0]), tuple(M[-1]), (R + 0.5 * G.TS, Hm)]
    o = [tuple(O[0]), tuple(O[-1]), (R + D, Hm)]
    if depth:
        v = [P(r, y, -1) for r, y in a]
    else:
        v = [P(r, y, 0) for r, y in a]
    for sg in (-1, 1):
        v += [P(r, y, sg) for r, y in mm] + [P(r, y, sg) for r, y in o]
    if depth:
        v += [P(r, y, 1) for r, y in a]          # 15-17: back axis
    A0, A1, A2 = (15, 16, 17) if depth else (0, 1, 2)
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    n1, g1 = graded(L["rows1"])
    n2, g2 = graded(L["rows2"])
    wr = f"{float(m['wallRatio'])!r}"
    nb = int(m["bufferColumns"])
    # outlet columns: shared (the cylinder's grading, as validated) or even
    # (uniform across the outlet, optionally a wider axis column: WP S.6)
    oc, fa = m.get("outletColumns", "shared"), float(m.get("axisColumn", 1.0))
    if oc == "shared" and fa == 1.0:
        hop_grading = f"simpleGrading ({wr} {g1} 1)"
    else:
        L0 = float(M[0][0])
        if oc == "shared":
            raise ValueError("mesh: axisColumn needs outletColumns: even")
        w0 = L0 / ncols
        if fa != 1.0:
            nx, xb = graded([(fa * w0, 1, 1.0), (L0 - fa * w0, ncols - 1, 1.0)])
        else:
            xb = "1"
        hop_grading = (f"edgeGrading ({xb} {wr} {wr} {xb}  {g1} {g1} {g1} {g1}  1 1 1 1)")
    blocks = (f"    hex (0 3 4 1 {A0} 9 10 {A1}) ({ncols} {n1} 1) {hop_grading}\n"
              f"    hex (1 4 5 2 {A1} 10 11 {A2}) ({ncols} {n2} 1) simpleGrading ({wr} {g2} 1)\n"
              f"    hex (3 6 7 4 9 12 13 10) ({nb} {n1} 1) simpleGrading (1 {g1} 1)\n"
              f"    hex (4 7 8 5 10 13 14 11) ({nb} {n2} 1) simpleGrading (1 {g2} 1)\n")

    def pl(i, j, pts, sg):
        inner = " ".join("({!r} {!r} {!r})".format(*P(r, y, sg)) for r, y in pts[1:-1])
        return f"    polyLine {i} {j} ({inner})\n"
    edges = pl(3, 4, M, -1) + pl(9, 10, M, 1) + pl(6, 7, O, -1) + pl(12, 13, O, 1)
    side = "empty" if depth else "wedge"
    patches = {
        "outlet": ("patch", f"(0 3 9 {A0})"),
        "bottom": ("patch", "(3 6 12 9)"),
        "outer": ("patch", "(6 7 13 12) (7 8 14 13)"),
        "top": ("patch", f"(2 5 11 {A2}) (5 8 14 11)"),
        "front": (side, "(0 1 4 3) (1 2 5 4) (3 4 7 6) (4 5 8 7)"),
        "back": (side, f"({A0} 9 10 {A1}) ({A1} 10 11 {A2}) (9 12 13 10) (10 13 14 11)"),
        "axis": ("symmetryPlane" if depth else "empty", f"(0 1 {A1} {A0}) (1 2 {A2} {A1})"),
    }
    pt = [f"    {k}\n    {{\n        type {ty};\n        faces ({f});\n    }}"
          for k, (ty, f) in patches.items()]
    comment = (f"// conicalHopper: R {R:g} m, outlet radius {G.Ro:g} m, hopper wall "
               f"{math.degrees(G.beta):g} deg from the vertical,\n"
               f"// cylinder {G.Hc:g} m; wedge half angle {math.degrees(G.theta):g} deg about y.\n"
               f"// Columns: {ncols} from the axis to the shell middle (wall/axis {wr}),\n"
               f"// {nb} outside it. Rows: {h_out * 1000:g} mm at the outlet growing to "
               f"{h_hop * 1000:g} mm at y = {y_t:g} m, then to {h_top * 1000:g} mm.\n")
    txt = (foam.header("dictionary", "blockMeshDict", src) + comment
           + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n{blocks});\n"
           f"edges\n(\n{edges});\n"
           "boundary\n(\n" + "\n".join(pt) + "\n);\n")
    info = dict(ncols=ncols, y_t=y_t, Hm=Hm, nRows=n1 + n2, nCells=(ncols + nb) * (n1 + n2))
    return txt, info


def shell_stl(G, y_top, ext=6.0, seg_len=0.005):
    """The shell polygon extruded across the wedge (beyond both wedge sides,
    decision 0093), caps by ear clipping."""
    poly = G.shell_polygon(y_top=y_top)
    dense = []
    n = len(poly)
    for i in range(n):
        a, bb = poly[i], poly[(i + 1) % n]
        k = max(1, int(math.ceil(np.linalg.norm(bb - a) / seg_len)))
        dense += [tuple(a + (bb - a) * j / k) for j in range(k)]
    return extruded_stl("shell", dense, G.theta, ext)


def extruded_stl(name, poly, theta, ext=6.0, wedge=True, depth=None):
    """Closed STL of an (r, y) polygon extruded across a wedge (z = +-ext r
    sin theta) or a plane mesh (z = +-depth)."""
    out = [f"solid {name}"]

    def P(rv, sg):
        r, y = rv
        if wedge:
            return np.array([r * math.cos(theta), y, sg * ext * r * math.sin(theta)])
        return np.array([r, y, sg * depth])

    def tri(p, q, w):
        nn = np.cross(q - p, w - p)
        nn /= np.linalg.norm(nn)
        out.append(f"facet normal {nn[0]:.9e} {nn[1]:.9e} {nn[2]:.9e}\nouter loop")
        for x in (p, q, w):
            out.append(f"vertex {x[0]:.12e} {x[1]:.12e} {x[2]:.12e}")
        out.append("endloop\nendfacet")
    poly = [tuple(p) for p in poly]
    n = len(poly)
    area = sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1]
               for i in range(n))
    if area < 0:
        poly = poly[::-1]
    for i in range(n):
        a, bb = poly[i], poly[(i + 1) % n]
        tri(P(a, -1), P(bb, -1), P(bb, 1))
        tri(P(a, -1), P(bb, 1), P(a, 1))
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
    out.append(f"endsolid {name}")
    return "\n".join(out) + "\n"


def build(cfg, src):
    """Everything template-specific: returns a dict with the files to write
    (relative path -> text), wall/outlet friction, patches, the fill and
    body polygons for the alpha fill, and notes for the checker."""
    g, w, m = cfg["geometry"], cfg["walls"], cfg["mesh"]
    for k in ("diameter", "cylinderHeight", "hopperAngle", "outletDiameter"):
        if g.get(k) is None:
            raise ValueError(f"geometry: '{k}' is required for conicalHopper")
    G = Geometry(g, w)
    y_fill = G.H if g["fillHeight"] in ("top", None) else float(g["fillHeight"])
    if not G.Y0 < y_fill <= G.H + 1e-9:
        raise ValueError(f"geometry: fillHeight must be between 0 and the silo top {G.H:.4f} m")
    bmd, info = block_mesh_dict(G, m, y_fill, src)
    files = {"constant/triSurface/shell.stl": shell_stl(G, info["Hm"] + 0.05)}
    if m.get("mesher", "blockMesh") == "gmsh":
        from . import gmsh_mesh
        if m.get("outletColumns", "shared") != "shared":
            raise ValueError("mesh: outletColumns needs mesher: blockMesh")
        files["system/mesh.geo"] = gmsh_mesh.geo_script(layout(G, m, y_fill), G)
    else:
        files["system/blockMeshDict"] = bmd
    return dict(
        files=files,
        wedge=True,
        theta=G.theta,
        mesher=m.get("mesher", "blockMesh"),
        patchTypes={"front": "wedge", "back": "wedge"},
        bodies={"shell": dict(surface="shell.stl", friction=float(w["friction"]))},
        walls={"outer": 0.0, "bottom": 0.0, "top": 0.0},
        wallsNote="domain sides outside the STL wall: frictionless",
        outlets={"outlet": float(w["friction"])},
        fill=G.fill_polygon(y_fill).tolist(),
        solids=[G.shell_polygon(y_top=info["Hm"] + 0.05).tolist()],
        volume=G.volume_below(y_fill),
        decompose="(1 {n} 1)",
        notes=dict(
            template="conicalHopper",
            outlet=dict(r0=0.0, r1=G.Ro, y=0.0),
            wall=G.wall_path(G.H).tolist(),
            siloTop=G.H, fillTop=y_fill, hopperHeight=G.ht,
            mesh=info,
        ),
        summary=(f"conical hopper: D {2 * G.R:g} m, outlet {2 * G.Ro:g} m, hopper "
                 f"{math.degrees(G.beta):g} deg from vertical ({G.ht:.3f} m tall), cylinder "
                 f"{G.Hc:g} m; fill to {y_fill:.3f} m; mesh {info['ncols']}+"
                 f"{int(m['bufferColumns'])} columns x {info['nRows']} rows "
                 f"= {info['nCells']} cells"),
    )
