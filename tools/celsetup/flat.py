"""Template flatBottom (WP G.5): a flat-bottomed silo or tank with a
central outlet, axisymmetric (a wedge) or plane (half of it, a symmetry
plane at x = 0). Sia's first target problem was the plane form: a
rectangular tank with a small bottom opening.

The walls are straight, so the mesh follows them without an STL shell: two
blocks, the outlet column (uniform cells of outletCell across the outlet)
and the rest (columns growing from outletCell at the outlet edge to
sideCell at the side wall). Rows grow from outletCell at the floor to
topCell over growthLength. Patches: outlet (floor under the outlet),
bottom (the rest of the floor, a wall), wall (the side), top, front, back,
axis. blockMesh or Gmsh.
"""

import math

from . import foam, schema
from .conical import graded, seg
from .schema import E

R_ = schema.REQUIRED
SCHEMA = {
    "geometry": [
        E("shape", "axisymmetric", "", "axisymmetric: a cylindrical silo on a wedge; plane: a "
          "tank in plane strain (diameters read as widths), half of it meshed with a symmetry "
          "plane", ["axisymmetric", "plane"]),
        E("diameter", R_, "m", "inner diameter (plane: width) of the silo"),
        E("height", R_, "m", "height of the silo (the mesh top)"),
        E("outletDiameter", R_, "m", "outlet diameter (plane: width), in the middle of the floor"),
        E("fillHeight", "top", "m", "material fill height above the floor, or top"),
        E("wedgeAngle", 5.0, "deg", "whole angle of the axisymmetric wedge (5 is tested)",
          level="advanced"),
        E("depth", 0.01, "m", "thickness of a plane mesh (one cell); results are per metre",
          level="advanced"),
    ],
    "walls": [
        E("frictionCoefficient", None, "", "wall friction coefficient mu_w (or give friction)"),
        E("friction", None, "deg", "wall friction angle, instead of frictionCoefficient"),
    ],
    "mesh": [
        E("preset", "standard", "", "standard: 4 mm at the outlet growing to 20 mm; fine: "
          "half of that", ["standard", "fine"]),
        E("mesher", "blockMesh", "", "blockMesh (OpenFOAM) or gmsh (Gmsh, GPL): the same "
          "blocks and cell counts", ["blockMesh", "gmsh"]),
        E("outletCell", "preset", "m", "cell size across the outlet and in the floor rows "
          "(3 to 5 mm on the 0.35 m Couto outlet set its rate); preset: from the mesh preset"),
        E("sideCell", "preset", "m", "column width at the side wall; preset: from the mesh "
          "preset"),
        E("topCell", "preset", "m", "row height in the upper part; preset: from the mesh preset"),
        E("growthLength", "preset", "m", "rows grow from outletCell to topCell over this "
          "height; preset: from the mesh preset", level="advanced"),
    ],
}
PRESETS = {
    "standard": dict(outletCell=0.004, sideCell=0.02, topCell=0.02, growthLength=0.3),
    "fine": dict(outletCell=0.002, sideCell=0.01, topCell=0.01, growthLength=0.3),
}


class Geometry:
    def __init__(self, g):
        self.plane = g.get("shape", "axisymmetric") == "plane"
        self.R = 0.5 * float(g["diameter"])
        self.Ro = 0.5 * float(g["outletDiameter"])
        self.H = float(g["height"])
        self.theta = math.radians(0.5 * float(g["wedgeAngle"]))
        self.depth = float(g.get("depth") or 0.01)
        if not 0 < self.Ro < self.R:
            raise ValueError("geometry: need 0 < outletDiameter < diameter")
        if self.H <= 0:
            raise ValueError("geometry: height must be positive")


def layout(G, m):
    """Column and row segments [(length, n, ratio)] of the two blocks."""
    h_o, h_s, h_t = float(m["outletCell"]), float(m["sideCell"]), float(m["topCell"])
    gl = float(m["growthLength"])
    nA = max(1, int(round(G.Ro / h_o)))
    colsA = [(G.Ro, nA, 1.0)]
    colsB = [seg(G.R - G.Ro, h_o, h_s)]
    rows = ([seg(gl, h_o, h_t), seg(G.H - gl, h_t, h_t)] if G.H > gl + 1e-9
            else [seg(G.H, h_o, h_t)])
    return dict(colsA=colsA, colsB=colsB, rows=rows)


def block_mesh_dict(G, L, src):
    co, si = math.cos(G.theta), math.sin(G.theta)

    def P(r, y, sg):
        if G.plane:
            return (r, y, sg * 0.5 * G.depth)
        return (r * co, y, sg * r * si)
    rz = [(0.0, 0.0), (G.Ro, 0.0), (G.R, 0.0), (0.0, G.H), (G.Ro, G.H), (G.R, G.H)]
    v = [P(r, y, -1) for r, y in rz] + [P(r, y, 1) for r, y in rz]
    # wedge: the back axis points are the front ones (collapsed)
    b0, b3 = (6, 9) if G.plane else (0, 3)
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    nA, gA = graded(L["colsA"])
    nB, gB = graded(L["colsB"])
    nR, gR = graded(L["rows"])
    blocks = (f"    hex (0 1 4 3 {b0} 7 10 {b3}) ({nA} {nR} 1) simpleGrading ({gA} {gR} 1)\n"
              f"    hex (1 2 5 4 7 8 11 10) ({nB} {nR} 1) simpleGrading ({gB} {gR} 1)\n")
    side = "empty" if G.plane else "wedge"
    patches = {
        "outlet": ("patch", f"(0 1 7 {b0})"),
        "bottom": ("patch", "(1 2 8 7)"),
        "wall": ("patch", "(2 5 11 8)"),
        "top": ("patch", f"(3 4 10 {b3}) (4 5 11 10)"),
        "front": (side, "(0 1 4 3) (1 2 5 4)"),
        "back": (side, f"({b0} {b3} 10 7) (7 10 11 8)"),
        "axis": ("symmetryPlane" if G.plane else "empty", f"(0 3 {b3} {b0})"),
    }
    pt = [f"    {k}\n    {{\n        type {ty};\n        faces ({f});\n    }}"
          for k, (ty, f) in patches.items()]
    txt = (foam.header("dictionary", "blockMeshDict", src)
           + f"// flatBottom ({'plane, half' if G.plane else 'axisymmetric wedge'}): "
           f"R {G.R:g} m, outlet {G.Ro:g} m, height {G.H:g} m\n"
           + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n{blocks});\nedges\n(\n);\n"
           "boundary\n(\n" + "\n".join(pt) + "\n);\n")
    return txt, dict(nA=nA, nB=nB, nRows=nR, nCells=(nA + nB) * nR)


def geo_script(G, L):
    """The same two blocks for Gmsh (built-in kernel): transfinite quads,
    one layer revolved (wedge) or extruded (plane)."""
    out = ["// written by celSetup (flatBottom, WP G.5)", 'SetFactory("Built-in");',
           "Mesh.MshFileVersion = 2.2;"]
    xs = [0.0, G.Ro, G.R]
    ys = [0.0]
    for ln, _, _ in L["rows"]:
        ys.append(ys[-1] + ln)
    ys[-1] = G.H
    cols = [L["colsA"][0], L["colsB"][0]]
    pid = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            pid[i, j] = len(pid) + 1
            out.append(f"Point({pid[i, j]}) = {{{x!r}, {y!r}, 0}};")
    cid = [0]

    def line(a, b, n, ratio):
        cid[0] += 1
        coef = 1.0 if n <= 1 or ratio == 1 else ratio ** (1.0 / (n - 1))
        out.append(f"Line({cid[0]}) = {{{a}, {b}}};")
        out.append(f"Transfinite Curve {{{cid[0]}}} = {n + 1} Using Progression {coef!r};")
        return cid[0]
    h = {(i, j): line(pid[i, j], pid[i + 1, j], cols[i][1], cols[i][2])
         for j in range(len(ys)) for i in range(2)}
    v = {(i, j): line(pid[i, j], pid[i, j + 1], L["rows"][j][1], L["rows"][j][2])
         for j in range(len(ys) - 1) for i in range(3)}
    surfs = []
    for j in range(len(ys) - 1):
        for i in range(2):
            s = len(surfs) + 1
            out.append(f"Curve Loop({s}) = {{{h[i, j]}, {v[i + 1, j]}, -{h[i, j + 1]}, "
                       f"-{v[i, j]}}};")
            out.append(f"Plane Surface({s}) = {{{s}}};")
            out.append(f"Transfinite Surface {{{s}}};\nRecombine Surface {{{s}}};")
            surfs.append((s, i, j))
    sl = ", ".join(str(s) for s, _, _ in surfs)
    if G.plane:
        out.append(f"Translate {{0, 0, {-0.5 * G.depth!r}}} {{ Surface{{{sl}}}; }}")
        out.append(f"ex[] = Extrude {{0, 0, {G.depth!r}}} {{ Surface{{{sl}}}; Layers{{1}}; "
                   "Recombine; };")
    else:
        out.append(f"Rotate {{{{0, 1, 0}}, {{0, 0, 0}}, {G.theta!r}}} {{ Surface{{{sl}}}; }}")
        out.append(f"ex[] = Extrude {{{{0, 1, 0}}, {{0, 0, 0}}, {-2 * G.theta!r}}} "
                   f"{{ Surface{{{sl}}}; Layers{{1}}; Recombine; }};")
    # per surface: top, volume, then sides lo, right, hi, left; revolving drops
    # the left side of the axis column (it collapses)
    start, k = {}, 0
    for s, i, j in surfs:
        start[s] = k
        k += 5 if (i == 0 and not G.plane) else 6
    groups = {"outlet": [], "bottom": [], "wall": [], "top": [], "axis": []}
    nj = len(ys) - 1
    for s, i, j in surfs:
        side = lambda q: f"ex[{start[s] + 2 + q}]"  # noqa: E731
        if j == 0:
            groups["outlet" if i == 0 else "bottom"].append(side(0))
        if j == nj - 1:
            groups["top"].append(side(2))
        if i == 1:
            groups["wall"].append(side(1))
        if i == 0 and G.plane:
            groups["axis"].append(side(3))
    for name, lst in groups.items():
        if lst:
            out.append(f'Physical Surface("{name}") = {{{", ".join(lst)}}};')
    out.append(f'Physical Surface("front") = {{{sl}}};')
    out.append(f'Physical Surface("back") = {{{", ".join(f"ex[{start[s]}]" for s, _, _ in surfs)}}};')
    out.append(f'Physical Volume("internal") = '
               f'{{{", ".join(f"ex[{start[s] + 1}]" for s, _, _ in surfs)}}};')
    return "\n".join(out) + "\n"


def build(cfg, src):
    g, w, m = cfg["geometry"], cfg["walls"], cfg["mesh"]
    G = Geometry(g)
    y_fill = G.H if g["fillHeight"] in ("top", None) else float(g["fillHeight"])
    if not 0 < y_fill <= G.H + 1e-9:
        raise ValueError(f"geometry: fillHeight must be between 0 and the height {G.H:g} m")
    if w.get("friction") is None:
        raise ValueError("walls: give friction (deg) or frictionCoefficient")
    fr = float(w["friction"])
    L = layout(G, m)
    bmd, info = block_mesh_dict(G, L, src)
    mesher = m.get("mesher", "blockMesh")
    if mesher == "gmsh":
        files = {"system/mesh.geo": geo_script(G, L)}
    else:
        files = {"system/blockMeshDict": bmd}
    vol = (2 * G.R * y_fill) if G.plane else math.pi * G.R ** 2 * y_fill
    return dict(
        files=files,
        wedge=not G.plane, theta=None if G.plane else G.theta,
        plane=dict(depth=G.depth, mirror=2.0) if G.plane else None,
        massUnit="kg/m" if G.plane else "kg",
        mesher=mesher,
        patchTypes=({"front": "empty", "back": "empty", "axis": "symmetryPlane"} if G.plane
                    else {"front": "wedge", "back": "wedge"}),
        bodies={},
        walls={"bottom": fr, "wall": fr, "top": 0.0},
        wallsNote="the floor and side wall with the wall friction; the top frictionless",
        outlets={"outlet": fr},
        fill=[[0.0, 0.0], [G.R, 0.0], [G.R, y_fill], [0.0, y_fill]],
        solids=[], volume=vol,
        decompose="(1 {n} 1)",
        notes=dict(template="flatBottom", outlet=dict(r0=0.0, r1=G.Ro, y=0.0),
                   wall=[[G.R, 0.0], [G.R, G.H]], siloTop=G.H, fillTop=y_fill, mesh=info),
        summary=(f"flat-bottomed {'tank (plane, per metre, half meshed): width' if G.plane else 'silo: D'} "
                 f"{2 * G.R:g} m, outlet {2 * G.Ro:g} m, height {G.H:g} m; fill to {y_fill:g} m; "
                 f"mesh {info['nA']}+{info['nB']} columns x {info['nRows']} rows "
                 f"= {info['nCells']} cells"),
    )
