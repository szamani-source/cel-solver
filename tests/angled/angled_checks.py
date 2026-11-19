"""WP G.1 done-when checks: the solver on angled (non-orthogonal) hex cells
(docs/spec/angled-cells-audit.md, decision 0103). About 1 minute.

[1] P-wave on a 45 degree sheared periodic mesh (celCycleCheck, fixed step):
    L1 error within 10 % of the Cartesian mesh at n 25 and 50, and second
    order (error ratio >= 3.5 from n 25 to 50).
[2] Elastic column under gravity in a box of vertical walls, meshed as two
    blocks split by an inclined line, so every column of cells leans and the
    cells differ in shape (celFoam, no free surface, solver's own step):
    - the first step equals maxCo times the smallest perpendicular height
      between opposite faces (V over the largest face, from the mesh) over
      the P-wave speed;
    - it settles (largest speed below 1e-3 m/s at the end);
    - frictionless walls: s_yy = -rho g (H - y) and s_xx = K s_yy, K =
      nu/(1 - nu), within 3 % of rho g H in the cells from 10 to 90 % height.

    python3 angled_checks.py WORKDIR
"""

import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "eulerian"))
import cycle_checks as cc  # noqa: E402

fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


# [1] sheared wave -----------------------------------------------------------

def sheared_wave(work, n, s):
    L, h = 1.0, 1.0 / n
    cells, size = (n, n, 1), (L, L, h)
    X, idx = cc.centres(cells, size)
    X[:, 0] += s * X[:, 1]
    k = (1, 0, 0)
    U, S = cc.wave_state(X, idx, n, "P", k, "sine", np.zeros(3))
    t_end = 1.0 / cc.CP
    nst = cc.steps_for(t_end, h / math.sqrt(1 + s * s), 0.4)
    cyc = dict(deltaT=t_end / nst, nSteps=nst, frameNormal=[1.0, 0.0, 0.0],
               streaming=[0.0, 0.0, 0.0], writeSteps=[nst])
    p = cc.make_case(os.path.join(work, f"wave_s{s:g}_n{n}"), cells, size, U, S, cyc)
    bm = open(f"{p}/system/blockMeshDict").read()
    for z in ("0", repr(h)):
        bm = bm.replace(f"({L} {L} {z}) (0 {L} {z})", f"({L + s} {L} {z}) ({s} {L} {z})")
    open(f"{p}/system/blockMeshDict", "w").write(bm)
    cc.run(p)
    cl = cc.cells_at(p, nst)
    t = cc.history(p)[-1, 1]
    phi = 2 * np.pi * cl[:, 0]
    exact = -cc.CP * cc.A * np.sin(phi - 2 * np.pi * cc.CP * t)
    return np.abs(cl[:, 3] - exact).mean() / (cc.CP * cc.A)


def check_wave(work):
    print("[1] P-wave on a 45 degree sheared mesh", flush=True)
    e = {(s, n): sheared_wave(work, n, s) for s in (0.0, 1.0) for n in (25, 50)}
    for n in (25, 50):
        r = e[1.0, n] / e[0.0, n]
        report(f"n {n} error against Cartesian", abs(r - 1) <= 0.10,
               f"sheared {e[1.0, n]:.3e}, Cartesian {e[0.0, n]:.3e}, ratio {r:.3f}")
    q = e[1.0, 25] / e[1.0, 50]
    report("second order on the sheared mesh", q >= 3.5, f"error ratio {q:.2f} (>= 3.5)")


# [2] leaning-cell column ----------------------------------------------------

W, H, NX, NY = 0.5, 1.0, 10, 40
X0, X1 = 0.1, 0.4              # split line from (X0, 0) to (X1, H)
RHO, E, NU, G, CO = 1500.0, 1.0e7, 0.3, 9.81, 0.4
HEADER = """FoamFile
{{
    version 2.0;
    format ascii;
    class {cls};
    object {obj};
}}
"""


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(text)


def column_case(p, t_end):
    d = W / NX
    v = [(0, 0, 0), (X0, 0, 0), (X1, H, 0), (0, H, 0), (W, 0, 0), (W, H, 0)]
    v += [(x, y, d) for x, y, _ in v]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    write(f"{p}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict")
          + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n"
          f"    hex (0 1 2 3 6 7 8 9) ({NX // 2} {NY} 1) simpleGrading (1 1 1)\n"
          f"    hex (1 4 5 2 7 10 11 8) ({NX // 2} {NY} 1) simpleGrading (1 1 1)\n);\n"
          "boundary\n(\n"
          "    left { type patch; faces ((0 6 9 3)); }\n"
          "    right { type patch; faces ((4 5 11 10)); }\n"
          "    bottom { type patch; faces ((0 1 7 6) (1 4 10 7)); }\n"
          "    top { type patch; faces ((3 9 8 2) (2 8 11 5)); }\n"
          "    frontAndBack { type empty; faces ((0 3 2 1) (1 2 5 4) (6 7 8 9) "
          "(7 10 11 8)); }\n);\n")
    write(f"{p}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict")
          + f"application celFoam;\nstartFrom startTime; startTime 0; stopAt endTime; "
          f"endTime {t_end!r}; deltaT 1e-5;\nwriteControl adjustableRunTime; "
          f"writeInterval {t_end!r}; writeFormat ascii; writePrecision 12;\n"
          f"timeFormat general; timePrecision 6;\nrunTimeModifiable false;\n"
          f"maxCo {CO!r};\ngravityRamp 0.05;\ndamping 200;\n")
    write(f"{p}/system/fvSchemes", HEADER.format(cls="dictionary", obj="fvSchemes")
          + "ddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
          "divSchemes { default none; }\nlaplacianSchemes { default none; }\n"
          "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    write(f"{p}/system/fvSolution", HEADER.format(cls="dictionary", obj="fvSolution"))
    walls = "".join(f"    {w} {{ friction 0; }}\n" for w in ("left", "right", "bottom", "top"))
    write(f"{p}/constant/materialProperties",
          HEADER.format(cls="dictionary", obj="materialProperties")
          + f"rho {RHO!r};\ntype linearElastic;\nE {E!r};\nnu {NU!r};\n"
          + "walls\n{\n" + walls + "}\n")
    write(f"{p}/constant/g", HEADER.format(cls="uniformDimensionedVectorField", obj="g")
          + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G!r} 0);\n")


def sh(p, cmd, log):
    with open(f"{p}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=p, stdout=f, stderr=subprocess.STDOUT,
                           env=dict(os.environ, PWD=p))
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed:\n" + open(f"{p}/log.{log}").read()[-3000:])


def foam_list(path, kind):
    """points, faces, owner or neighbour from an ascii polyMesh file."""
    s = open(path).read()
    s = s[s.index("}", s.index("FoamFile")) + 1:]
    s = re.sub(r"//.*", "", s)
    n = int(re.search(r"^\s*(\d+)\s*\(", s, re.M).group(1))
    body = s[s.index("(") + 1:]
    if kind == "label":
        return np.array(body.split(")")[0].split(), int)[:n]
    if kind == "point":
        return np.array(re.findall(r"\(([^()]*)\)", body)[:n] and
                        [list(map(float, m.split())) for m in
                         re.findall(r"\(([^()]*)\)", body)[:n]])
    return [list(map(int, m.split())) for m in re.findall(r"\d+\(([^()]*)\)", body)[:n]]


def smallest_height(p):
    """min over cells of V / (largest |S| of the cell's faces, empty sides left out)."""
    m = f"{p}/constant/polyMesh"
    pts = foam_list(f"{m}/points", "point")
    faces = foam_list(f"{m}/faces", "face")
    own = foam_list(f"{m}/owner", "label")
    nei = foam_list(f"{m}/neighbour", "label")
    nc = own.max() + 1
    V = np.zeros(nc)
    Smax = np.zeros(nc)
    empty = set()
    b = open(f"{m}/boundary").read()
    for name, body in re.findall(r"(\w+)\s*\{([^{}]*)\}", b):
        if re.search(r"type\s+empty", body):
            n0 = int(re.search(r"nFaces\s+(\d+)", body).group(1))
            s0 = int(re.search(r"startFace\s+(\d+)", body).group(1))
            empty.update(range(s0, s0 + n0))
    for f, vs in enumerate(faces):
        q = pts[vs]
        c = q.mean(axis=0)
        S = 0.5 * sum(np.cross(q[i] - c, q[(i + 1) % len(q)] - c) for i in range(len(q)))
        v = np.dot(S, c) / 3.0
        cells = [(own[f], v)] + ([(nei[f], -v)] if f < len(nei) else [])
        for ci, vv in cells:
            V[ci] += vv
            if f not in empty:
                Smax[ci] = max(Smax[ci], np.linalg.norm(S))
    return (V / Smax).min()


def read_field(path):
    """Internal field of an ascii vector or symmTensor field, one row per cell."""
    s = open(path).read()
    s = s[s.index("internalField"):]
    n = int(re.search(r"List<\w+>\s*(\d+)", s).group(1))
    rows = re.findall(r"\(([^()]*)\)", s[s.index("(", s.index("List")):])[:n]
    return np.array([list(map(float, r.split())) for r in rows])


def check_column(work):
    print("[2] Elastic column on leaning cells (celFoam, solver's own step)", flush=True)
    p = os.path.join(work, "column")
    t_end = 0.4
    column_case(p, t_end)
    sh(p, ["blockMesh"], "blockMesh")
    sh(p, ["celFoam"], "celFoam")
    sh(p, ["postProcess", "-func", "writeCellCentres", "-time", f"{t_end:g}"], "centres")
    log = open(f"{p}/log.celFoam").read()
    dt0 = float(re.search(r"Stable time step\s+deltaT = (\S+)", log).group(1))
    M = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    h = smallest_height(p)
    want = CO * h / math.sqrt(M / RHO)
    report("first step from the perpendicular cell height", abs(dt0 / want - 1) < 1e-6,
           f"deltaT {dt0:.6e}, maxCo h/c {want:.6e}, h {h * 1e3:.2f} mm")
    tdir = f"{p}/{t_end:g}"
    U = read_field(f"{tdir}/U")
    S = read_field(f"{tdir}/sigma")
    C = read_field(f"{tdir}/C")
    um = float(np.abs(U).max())
    report("settled", np.isfinite(um) and um < 1e-3, f"largest speed {um:.2e} m/s")
    y = C[:, 1]
    band = (y > 0.1 * H) & (y < 0.9 * H)
    syy = -RHO * G * (H - y[band])
    K = NU / (1 - NU)
    ref = RHO * G * H
    eyy = np.abs(S[band, 3] - syy).max() / ref
    exx = np.abs(S[band, 0] - K * syy).max() / ref
    report("s_yy hydrostatic", eyy <= 0.03, f"largest error {eyy * 100:.2f} % of rho g H")
    report("s_xx = K s_yy", exx <= 0.03, f"largest error {exx * 100:.2f} % of rho g H")


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    check_wave(work)
    check_column(work)
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
