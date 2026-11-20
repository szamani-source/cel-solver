"""WP G.3 done-when checks: walls and free surface on angled faces
(docs/spec/angled-cells-audit.md, decision 0106). About 2 minutes.

An infinite 30 degree slope: a strip of cells along the slope, periodic
(cyclic) at its two ends, with a frictional wall underneath and void above
an elastic layer 60 mm thick (6 cells). Gravity (0, -g, 0) with a cosine
ramp of T = 0.05 s. With wall friction 15 deg a rigid slide accelerates at
a = g (sin 30 - tan 15 cos 30) = 2.63 m/s^2, speed a (t - T/2) after the
ramp; with friction 60 deg it stays at rest.

Three meshes of 10 mm cells:
  cartesian  the strip untilted, gravity turned instead (reference);
  turned     the strip turned 30 deg: no face along an axis, cells square;
  sheared    lines along the slope and vertical lines (parallelograms).

[walls] celFoam without a free surface: the strip full, a frictionless lid
on top: every cell's speed down the slope within 1 % of the closed form,
normal speed below 1 % of it (wall friction on angled faces).
[free surface] celVoidCheck, layer and void:
for each mesh: mean speed down the slope within 2 % of the closed form at
0.3 s, mean speed normal to the slope below 2 % of it, mass to 1e-10; at
friction 60 deg, mean speed below 1 % of the free slide.

    python3 slope_checks.py WORKDIR
"""

import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
import void_checks as vc  # noqa: E402
cs = vc.cs

THETA, FRICTION, G, RAMP, T_END = 30.0, 15.0, 9.81, 0.05, 0.3
LS, HN, NS, NN = 0.4, 0.1, 40, 10           # strip along and normal to the slope
A0, A1 = 0.15, 0.25                         # split line (A0, 0) to (A1, HN)
th = math.radians(THETA)
ACC = G * (math.sin(th) - math.tan(math.radians(FRICTION)) * math.cos(th))
V_EXACT = ACC * (T_END - RAMP / 2)
fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def strip_case(p, fill, tilted):
    """Strip of NS x NN cells; tilted: turned by -THETA, split by a leaning
    line, gravity -y; "sheared": vertical lines and lines along the slope,
    gravity -y; otherwise one Cartesian block with gravity turned."""
    d = LS / NS
    if tilted == "sheared":
        tn = math.tan(th)
        lx, hv = LS * math.cos(th), HN / math.cos(th)
        loc = [(0, 0), (lx / 2, 0), (lx / 2, hv), (0, hv), (lx, 0), (lx, hv)]
        v = [(x, y - x * tn) for x, y in loc]
        n1 = NS // 2
        blocks = (f"    hex (0 1 2 3 6 7 8 9) ({n1} {NN} 1) simpleGrading (1 1 1)\n"
                  f"    hex (1 4 5 2 7 10 11 8) ({NS - n1} {NN} 1) simpleGrading (1 1 1)\n")
        bottom, top = "(0 1 7 6) (1 4 10 7)", "(3 9 8 2) (2 8 11 5)"
        empty = "(0 3 2 1) (1 2 5 4) (6 7 8 9) (7 10 11 8)"
        gvec = (0.0, -G)
    elif tilted:
        loc = [(0, 0), (A0, 0), (A1, HN), (0, HN), (LS, 0), (LS, HN)]
        ca, sa = math.cos(-th), math.sin(-th)
        v = [(a * ca - b * sa, a * sa + b * ca) for a, b in loc]
        n1 = NS // 2
        blocks = (f"    hex (0 1 2 3 6 7 8 9) ({n1} {NN} 1) simpleGrading (1 1 1)\n"
                  f"    hex (1 4 5 2 7 10 11 8) ({NS - n1} {NN} 1) simpleGrading (1 1 1)\n")
        bottom, top = "(0 1 7 6) (1 4 10 7)", "(3 9 8 2) (2 8 11 5)"
        empty = "(0 3 2 1) (1 2 5 4) (6 7 8 9) (7 10 11 8)"
        gvec = (0.0, -G)
    else:
        v = [(0, 0), (LS / 2, 0), (LS / 2, HN), (0, HN), (LS, 0), (LS, HN)]
        n1 = NS // 2
        blocks = (f"    hex (0 1 2 3 6 7 8 9) ({n1} {NN} 1) simpleGrading (1 1 1)\n"
                  f"    hex (1 4 5 2 7 10 11 8) ({NS - n1} {NN} 1) simpleGrading (1 1 1)\n")
        bottom, top = "(0 1 7 6) (1 4 10 7)", "(3 9 8 2) (2 8 11 5)"
        empty = "(0 3 2 1) (1 2 5 4) (6 7 8 9) (7 10 11 8)"
        gvec = (G * math.sin(th), -G * math.cos(th))
    pts = [(x, y, 0.0) for x, y in v] + [(x, y, d) for x, y in v]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in pts)
    H_ = cs.HEADER
    cs.write(f"{p}/system/blockMeshDict", H_.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n{blocks});\n"
             "boundary\n(\n"
             "    upstream { type cyclic; neighbourPatch downstream; faces ((0 6 9 3)); }\n"
             "    downstream { type cyclic; neighbourPatch upstream; faces ((4 5 11 10)); }\n"
             f"    bottom {{ type patch; faces ({bottom}); }}\n"
             f"    top {{ type patch; faces ({top}); }}\n"
             f"    frontAndBack {{ type empty; faces ({empty}); }}\n);\n")
    cs.write(f"{p}/system/controlDict", H_.format(cls="dictionary", obj="controlDict")
             + "application celVoidCheck;\nstartFrom startTime; startTime 0; stopAt endTime;"
             " endTime 100; deltaT 1;\nwriteControl timeStep; writeInterval 100000000;"
             f" writeFormat ascii;\ngravityRamp {RAMP!r};\n")
    cs.write(f"{p}/system/fvSchemes", H_.format(cls="dictionary", obj="fvSchemes")
             + "ddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
             "divSchemes { default none; }\nlaplacianSchemes { default none; }\n"
             "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    cs.write(f"{p}/system/fvSolution", H_.format(cls="dictionary", obj="fvSolution"))
    cs.write(f"{p}/constant/materialProperties",
             H_.format(cls="dictionary", obj="materialProperties")
             + f"rho {vc.RHO};\nfreeSurface yes;\n" + vc.ELASTIC
             + f"walls\n{{\n    bottom {{ friction {FRICTION!r}; }}\n}}\n")
    cs.write(f"{p}/constant/g", H_.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue ({gvec[0]!r} {gvec[1]!r} 0);\n")
    check = dict(fill=list(fill), U0=[0, 0, 0], tEnd=T_END, maxCo=0.4)
    if tilted == "sheared":
        check["fillSlope"] = THETA
        check["fill"] = [fill[0], fill[1], fill[2], fill[3] / math.cos(th)] \
            if fill[3] < HN else list(fill)
    elif tilted:
        check["fillAngle"] = -THETA
    body = "\n".join(f"{k} {cs.foam_value(x)};" for k, x in check.items())
    cs.write(f"{p}/system/voidCheckDict",
             H_.format(cls="dictionary", obj="voidCheckDict") + body + "\n")
    return p


def slide(work, name, fill, tilted):
    p = strip_case(os.path.join(work, name), fill, tilted)
    vc.run(p)
    h = vc.history(p)
    s = np.array([math.cos(th), -math.sin(th)]) if tilted else np.array([1.0, 0.0])
    n = np.array([-s[1], s[0]])
    m = h[-1, vc.MASS]
    P = h[-1, [vc.PX, vc.PY]] / m
    dm = abs(m / h[0, vc.MASS] - 1)
    return float(P @ s), float(P @ n), dm, h[-1, vc.T]


def lid_slide(work, name, kind):
    """Full strip under a frictionless lid, celFoam, no free surface."""
    p = strip_case(os.path.join(work, name), (-1, -1, LS + 1, HN + 1), kind)
    mp = open(f"{p}/constant/materialProperties").read()
    mp = mp.replace("freeSurface yes;", "freeSurface no;").replace(
        "walls\n{\n", "walls\n{\n    top { friction 0; }\n")
    open(f"{p}/constant/materialProperties", "w").write(mp)
    cs.write(f"{p}/system/controlDict", cs.HEADER.format(cls="dictionary", obj="controlDict")
             + f"application celFoam;\nstartFrom startTime; startTime 0; stopAt endTime; "
             f"endTime {T_END!r}; deltaT 1e-5;\nwriteControl adjustableRunTime; "
             f"writeInterval {T_END!r}; writeFormat ascii; writePrecision 12;\n"
             "timeFormat general; timePrecision 6;\nrunTimeModifiable false;\n"
             f"maxCo 0.4;\ngravityRamp {RAMP!r};\n")
    for cmd, log in ((["blockMesh"], "blockMesh"), (["celFoam"], "celFoam")):
        with open(f"{p}/log.{log}", "w") as f:
            r = subprocess.run(cmd, cwd=p, stdout=f, stderr=subprocess.STDOUT,
                               env=dict(os.environ, PWD=p))
        if r.returncode != 0:
            raise RuntimeError(open(f"{p}/log.{log}").read()[-3000:])
    U = read_vectors(f"{p}/{T_END:g}/U")
    sd = np.array([math.cos(th), -math.sin(th), 0.0]) if kind else np.array([1.0, 0, 0])
    nd = np.array([-sd[1], sd[0], 0.0])
    want = ACC * (T_END - RAMP / 2)
    es = np.abs(U @ sd / want - 1).max()
    en = np.abs(U @ nd).max() / want
    report(f"{name}: every cell's speed", es <= 0.01 and en <= 0.01,
           f"largest error down the slope {es * 100:.2f} %, normal {en * 100:.2f} % "
           f"(closed form {want:.4f} m/s)")


def read_vectors(path):
    s = open(path).read()
    s = s[s.index("internalField"):]
    n = int(re.search(r"List<\w+>\s*(\d+)", s).group(1))
    rows = re.findall(r"\(([^()]*)\)", s[s.index("(", s.index("List")):])[:n]
    return np.array([list(map(float, r.split())) for r in rows])


def check(name, vs, vn, dm, t, tol):
    want = ACC * (t - RAMP / 2)
    e = abs(vs / want - 1)
    report(f"{name}: speed down the slope", e <= tol,
           f"{vs:.4f} m/s at {t:.3f} s, closed form {want:.4f}, error {e * 100:.2f} %")
    report(f"{name}: speed normal to the slope", abs(vn) <= 0.02 * abs(want),
           f"{vn:.2e} m/s")
    report(f"{name}: mass", dm <= 1e-10, f"relative change {dm:.1e}")


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    global FRICTION
    print(f"Slope {THETA:g} deg, wall friction {FRICTION:g} deg: a = {ACC:.4f} m/s^2")
    layer = (-1, -1, LS + 1, 0.06)
    global A0, A1
    A0 = A1 = LS / 2                        # turned strip: square cells
    print("[walls] full strip under a frictionless lid, no free surface", flush=True)
    FRICTION = 15.0
    for name, kind in (("turned", True), ("sheared", "sheared")):
        lid_slide(work, f"{name}_lid", kind)
    for name, kind in (("cartesian", False), ("turned", True), ("sheared", "sheared")):
        print(f"[{name}]", flush=True)
        FRICTION = 15.0
        check(name, *slide(work, f"{name}_slide", layer, kind), 0.02)
        FRICTION = 60.0
        vs, _, dm, t = slide(work, f"{name}_stick", layer, kind)
        free = ACC * (t - RAMP / 2)
        report(f"{name}: at rest at friction 60 deg", abs(vs) <= 0.01 * free and dm <= 1e-10,
               f"speed {vs:.2e} m/s, mass change {dm:.1e}")
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
