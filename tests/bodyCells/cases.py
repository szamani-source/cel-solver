"""WP 4.2: build the celBodyCells test cases and check their output.

    cases.py make <dir>     write the cases plane, sphere10, sphere20, cube
    cases.py check <dir>    check the fields celBodyCells wrote
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "prototype"))
import contact2d as c2d  # noqa: E402

HDR = """FoamFile {{ version 2.0; format ascii; class dictionary; object {}; }}
"""


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def stl(path, tris):
    out = ["solid body"]
    for t in tris:
        n = np.cross(t[1] - t[0], t[2] - t[0])
        n = n / np.linalg.norm(n)
        out += [f"facet normal {n[0]:.17g} {n[1]:.17g} {n[2]:.17g}", "outer loop"]
        out += [f"vertex {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}" for p in t]
        out += ["endloop", "endfacet"]
    out.append("endsolid body")
    write(path, "\n".join(out) + "\n")


def box_tris(lo, hi, R=np.eye(3), shift=np.zeros(3)):
    """Closed box [lo, hi] (outward normals), rotated by R then shifted."""
    V = np.array([[x, y, z] for z in (lo[2], hi[2]) for y in (lo[1], hi[1])
                  for x in (lo[0], hi[0])], float)
    V = V @ R.T + shift
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    t = []
    for a, b, c, d in quads:
        t += [(V[a], V[b], V[c]), (V[a], V[c], V[d])]
    return t


def seam_box_tris(lo, hi, xs):
    """Closed box [lo, hi] whose faces normal to y and z are split into strips
    at the x values xs (a structured surface, as a hex8 body's): the strip
    edges are triangle edges along z and y at each x in xs."""
    X = [lo[0]] + list(xs) + [hi[0]]
    t = []
    for a, b in zip(X[:-1], X[1:]):
        for y, sy in ((lo[1], -1), (hi[1], 1)):        # faces normal to y
            q = [np.array(v, float) for v in ((a, y, lo[2]), (b, y, lo[2]),
                                              (b, y, hi[2]), (a, y, hi[2]))]
            if sy > 0:
                q = q[::-1]
            t += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
        for z, sz in ((lo[2], -1), (hi[2], 1)):        # faces normal to z
            q = [np.array(v, float) for v in ((a, lo[1], z), (b, lo[1], z),
                                              (b, hi[1], z), (a, hi[1], z))]
            if sz < 0:
                q = q[::-1]
            t += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
    for x, sx in ((lo[0], -1), (hi[0], 1)):            # the two ends
        q = [np.array(v, float) for v in ((x, lo[1], lo[2]), (x, hi[1], lo[2]),
                                          (x, hi[1], hi[2]), (x, lo[1], hi[2]))]
        if sx < 0:
            q = q[::-1]
        t += [(q[0], q[1], q[2]), (q[0], q[2], q[3])]
    for tri in t:                                      # outward normals
        c = np.mean(tri, axis=0)
        n = np.cross(tri[1] - tri[0], tri[2] - tri[0])
        mid = (np.array(lo) + np.array(hi)) / 2
        assert np.dot(n, c - mid) > 0
    return t


def icosphere(r, levels, centre):
    p = (1 + 5 ** 0.5) / 2
    V = [(-1, p, 0), (1, p, 0), (-1, -p, 0), (1, -p, 0), (0, -1, p), (0, 1, p),
         (0, -1, -p), (0, 1, -p), (p, 0, -1), (p, 0, 1), (-p, 0, -1), (-p, 0, 1)]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    F = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4),
         (11, 10, 2), (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8),
         (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(levels):
        cache, F2 = {}, []

        def mid(a, b):
            k = (min(a, b), max(a, b))
            if k not in cache:
                m = V[a] + V[b]
                V.append(m / np.linalg.norm(m))
                cache[k] = len(V) - 1
            return cache[k]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            F2 += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        F = F2
    P = [r * v + centre for v in V]
    return [(P[a], P[b], P[c]) for a, b, c in F]


def stl_volume(path):
    """Volume of a closed STL as OpenFOAM reads it (single precision)."""
    v = [list(map(float, l.split()[1:4])) for l in open(path) if l.strip().startswith("vertex")]
    P = np.array(v, dtype=np.float32).astype(float).reshape(-1, 3, 3)
    return np.einsum("ij,ij->i", P[:, 0], np.cross(P[:, 1], P[:, 2])).sum() / 6.0


def case(d, n, L, bodies, times="0"):
    """Box [0, L]^3 with n^3 cells; bodies: list of (name, tris, extra dict text)."""
    write(f"{d}/system/blockMeshDict", HDR.format("blockMeshDict") + f"""
vertices ( (0 0 0) ({L} 0 0) ({L} {L} 0) (0 {L} 0) (0 0 {L}) ({L} 0 {L}) ({L} {L} {L}) (0 {L} {L}) );
blocks ( hex (0 1 2 3 4 5 6 7) ({n[0]} {n[1]} {n[2]}) simpleGrading (1 1 1) );
boundary ( walls {{ type wall; faces ( (0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7) ); }} );
""")
    write(f"{d}/system/controlDict", HDR.format("controlDict") + """
application celBodyCells; startFrom startTime; startTime 0; stopAt endTime; endTime 1;
deltaT 1; writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 17;
timeFormat general; timePrecision 8;
""")
    write(f"{d}/system/fvSchemes", HDR.format("fvSchemes") + """
ddtSchemes { default Euler; } gradSchemes { default Gauss linear; }
divSchemes { default none; } laplacianSchemes { default none; }
interpolationSchemes { default linear; } snGradSchemes { default corrected; }
""")
    write(f"{d}/system/fvSolution", HDR.format("fvSolution"))
    write(f"{d}/system/decomposeParDict", HDR.format("decomposeParDict") +
          "numberOfSubdomains 2; method simple; coeffs { n (2 1 1); }\n")
    entries = []
    for name, tris, extra in bodies:
        stl(f"{d}/constant/triSurface/{name}.stl", tris)
        entries.append(f"    {name} {{ surface {name}.stl; {extra} }}")
    write(f"{d}/constant/rigidBodies", HDR.format("rigidBodies") +
          "bodies\n{\n" + "\n".join(entries) + "\n}\n")
    for t in times.split():
        os.makedirs(f"{d}/{t}", exist_ok=True)


# ---- the cases ------------------------------------------------------------------

PLANE_TH = 30.0
PLANE_C = 0.4321


def make(root):
    L = 1.0
    # plane: half-space n.x <= c, n = (-sin 30, cos 30, 0), as a big rotated box
    th = math.radians(PLANE_TH)
    n = np.array([-math.sin(th), math.cos(th), 0.0])
    t = np.array([math.cos(th), math.sin(th), 0.0])
    R = np.column_stack([t, n, [0, 0, 1]])        # local (u, v, w) -> world
    shift = PLANE_C * n
    case(f"{root}/plane", (10, 10, 2), L,
         [("plane", box_tris((-5, -5, -1), (5, 0, 2), R, shift), "")])
    ctr = np.array([0.5123, 0.4871, 0.5037])
    case(f"{root}/sphere10", (17, 17, 17), 1.02,
         [("sphere", icosphere(0.3, 3, ctr), "")])
    case(f"{root}/sphere20", (34, 34, 34), 1.02,
         [("sphere", icosphere(0.3, 3, ctr), "")])
    # cube 0.25 m, moving one cell (0.1 m) along x per second and turning about z
    case(f"{root}/cube", (10, 10, 4), 1.0,
         [("cube", box_tris((0.2117, 0.3089, 0.1043), (0.4617, 0.5589, 0.3543)),
           "centre (0.3367 0.4339 0.2293); velocity constant (0.1 0 0);")],
         times="0 1")
    case(f"{root}/cubeTurn", (10, 10, 4), 1.0,
         [("cube", box_tris((0.2117, 0.3089, 0.1043), (0.4617, 0.5589, 0.3543)),
           "centre (0.3367 0.4339 0.2293); omega constant (0 0 0.7);")],
         times="0 1")
    # seam (decision 0093): a box whose surface strips meet on the mesh face
    # x = 0.5, as the nodes of a hex8 body can
    case(f"{root}/seam", (10, 10, 4), 1.0,
         [("seam", seam_box_tris((0.2117, 0.3089, 0.1043), (0.7883, 0.5589, 0.3543),
                                 [0.35, 0.5, 0.65]), "")])


# ---- reading OpenFOAM ascii fields ------------------------------------------------

def read_internal(path, vector=False):
    s = open(path).read()
    i = s.index("internalField")
    tail = s[i:]
    if tail.split()[1] == "uniform":
        return None
    j = tail.index("(")
    n = int(tail[:j].split()[-1])
    body = tail[j + 1:]
    if vector:
        vals = []
        k = 0
        for _ in range(n):
            a = body.index("(", k)
            b = body.index(")", a)
            vals.append([float(x) for x in body[a + 1:b].split()])
            k = b + 1
        return np.array(vals)
    end = body.index(")")
    return np.array([float(x) for x in body[:end].split()])


fails = []


def check(name, err, tol):
    ok = err <= tol
    print(f"{'ok  ' if ok else 'FAIL'} {name}: {err:.3e} (limit {tol:g})")
    if not ok:
        fails.append(name)


def checks(root):
    # plane: every cell fraction exact (2D clip; cells span the body in z)
    d = f"{root}/plane"
    a = read_internal(f"{d}/0/alphaBody")
    W = read_internal(f"{d}/0/wallArea", vector=True)
    h, nz = 0.1, 2
    th = math.radians(PLANE_TH)
    n2 = (-math.sin(th), math.cos(th))
    hp = c2d.HalfPlane(n2, PLANE_C)
    g = c2d.Grid(10, 10, h)
    bex, _, _, _ = hp.fractions(g, 0.0)
    ex = np.tile(bex.ravel(), nz)             # blockMesh order: x fastest, then y, z
    # the STL holds single-precision vertices (about 3e-8 m here), so the plane
    # read is not exactly n.x = c: limits 1e-6 (the sphere and cube checks below
    # compare against the STL's own geometry and are exact to round-off)
    check("plane: cell fractions (single-precision STL)", float(np.abs(a - ex).max()), 1e-6)
    # wall face of cut cells: |W| = cut length x depth, direction into the body (-n)
    cut = (ex > 1e-12) & (ex < 1 - 1e-12)
    dz = 1.0 / nz
    errW = 0.0
    for idx in np.nonzero(cut)[0]:
        i, j = idx % 10, (idx // 10) % 10
        x0, y0 = i * h, j * h
        sq = [(x0, y0), (x0 + h, y0), (x0 + h, y0 + h), (x0, y0 + h)]
        # length of the line n.x = c inside the square
        pts = []
        for k in range(4):
            p, q = sq[k], sq[(k + 1) % 4]
            fp = n2[0] * p[0] + n2[1] * p[1] - PLANE_C
            fq = n2[0] * q[0] + n2[1] * q[1] - PLANE_C
            if fp * fq < 0:
                tt = fp / (fp - fq)
                pts.append((p[0] + tt * (q[0] - p[0]), p[1] + tt * (q[1] - p[1])))
        Lc = math.dist(pts[0], pts[1]) if len(pts) == 2 else 0.0
        exW = -Lc * dz * np.array([n2[0], n2[1], 0.0])
        errW = max(errW, float(np.abs(W[idx] - exW).max()))
    check("plane: wall face area vector (single-precision STL)", errW, 1e-6 * h * dz)

    # spheres: covered volume = volume of the faceted sphere
    for nm in ("sphere10", "sphere20"):
        d = f"{root}/{nm}"
        n = 17 if nm == "sphere10" else 34
        a = read_internal(f"{d}/0/alphaBody")
        V = (1.02 / n) ** 3
        Vs = stl_volume(f"{d}/constant/triSurface/sphere.stl")
        check(f"{nm}: covered volume against the faceted sphere (relative)",
              abs(a.sum() * V / Vs - 1), 1e-10)
        check(f"{nm}: faceted against true sphere (relative, information)",
              abs(Vs / (4 / 3 * math.pi * 0.3 ** 3) - 1), 0.01)
        if os.path.isdir(f"{d}/parallel"):
            ap = read_internal(f"{d}/parallel/0/alphaBody")
            check(f"{nm}: 2 processors against serial", float(np.abs(ap - a).max()), 1e-12)

    # moving cube: one cell along x in 1 s; fractions shift by exactly one cell
    d = f"{root}/cube"
    a0 = read_internal(f"{d}/0/alphaBody").reshape(4, 10, 10)
    a1 = read_internal(f"{d}/1/alphaBody").reshape(4, 10, 10)
    check("cube: fractions at t = 1 are those at t = 0 one cell along x",
          float(np.abs(a1[:, :, 1:] - a0[:, :, :-1]).max()), 1e-12)
    Vc = stl_volume(f"{d}/constant/triSurface/cube.stl")
    check("cube: covered volume at t = 0", abs(a0.sum() * 0.1 * 0.1 * 0.25 / Vc - 1), 1e-10)
    # turning cube: volume kept while it turns 0.7 rad
    d = f"{root}/cubeTurn"
    a1 = read_internal(f"{d}/1/alphaBody")
    check("turning cube: covered volume at t = 1", abs(a1.sum() * 0.1 * 0.1 * 0.25 / Vc - 1),
          1e-10)

    # seam: cells 3 to 6 along x see the same body; the strip edges on the
    # mesh face x = 0.5 must not change their fractions (decision 0093)
    d = f"{root}/seam"
    a = read_internal(f"{d}/0/alphaBody").reshape(4, 10, 10)
    check("seam: fractions either side of the strip edge on a mesh face",
          float(np.abs(a[:, :, 3:7] - a[:, :, 3:4]).max()), 1e-12)
    Vs = stl_volume(f"{d}/constant/triSurface/seam.stl")
    check("seam: covered volume", abs(a.sum() * 0.1 * 0.1 * 0.25 / Vs - 1), 1e-10)

    if fails:
        print("FAILED:", ", ".join(fails))
        sys.exit(1)
    print("all body cell checks passed")


if __name__ == "__main__":
    {"make": make, "check": checks}[sys.argv[1]](sys.argv[2])
