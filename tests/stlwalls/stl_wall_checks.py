"""WP 4.3a done-when checks: fixed STL walls with friction as cut-cell walls
(decision 0042, 0047). Run by tests/test_stl_walls.sh after celSiloCheck is
built:

    python3 stl_wall_checks.py WORKDIR

The WP 3.1 frictional tank (decision 0030, tests/silo/silo_checks.py):
0.4 x 1.2 m, T0 material (cohesionless Mohr-Coulomb, rho 1800, E 10 MPa,
nu 0.3, phi 35), wall friction 20 deg, 10 mm cells, gravity ramp 0.1 s and
damping 250 1/s until KE <= 1e-6 strain energy, then 10 % more steps
undamped.

1. Reference: side walls and bottom are mesh patches (as WP 3.1).
2. Aligned: the side walls are STL boxes aligned with the mesh; their inner
   faces cut cells at a quarter (left: open 0.75, a wall cell) and three
   quarters (right: open 0.25, a small cell mixed with its neighbour); the
   bottom stays a mesh wall. Mean vertical stress within 2 % of the
   reference at every depth >= W.
3. Tilted: the whole tank (walls, bottom, a smooth lid at the top that the
   settling material leaves, as a free top) is STL boxes turned 10 deg
   against the mesh, gravity turned with it. Mean stress along the tank axis
   within 5 % of the reference at every depth >= W.
In 2 and 3: material that entered through the STL walls (not the lid) below
1e-6 of the total; weight = wall forces along gravity (bodies plus mesh
walls) within 1 %; settled and steady as in WP 3.1.

Prints PASS/FAIL per check, exits non-zero on any failure.
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "silo"))
sys.path.insert(0, os.path.join(HERE, "..", "lagrange"))
import cases as cs  # noqa: E402
import silo_checks as sc  # noqa: E402

W, H = 0.4, 1.2
CELL = sc.CELL
DELTA_W = sc.DELTA_W
TILT = 10.0
SETTLE = dict(mode="settle", keTol=sc.KE_TOL, checkEvery=20, tMax=1.0,
              continueFraction=0.1)


def box_stl(lo, hi, rot=0.0):
    """ASCII STL of the box lo..hi (outward normals), turned by rot deg about z
    through the origin."""
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    p = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6),
             (0, 4, 7, 3), (1, 2, 6, 5)]
    out = ["solid box"]
    for q in quads:
        for t in ((q[0], q[1], q[2]), (q[0], q[2], q[3])):
            a, b, d = (R @ p[i] for i in t)
            n = np.cross(b - a, d - a)
            n /= np.linalg.norm(n)
            out.append(f"facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}")
            out.append("outer loop")
            for v in (a, b, d):
                out.append(f"vertex {v[0]:.12e} {v[1]:.12e} {v[2]:.12e}")
            out.append("endloop")
            out.append("endfacet")
    out.append("endsolid box")
    return "\n".join(out) + "\n"


def make_stl_case(path, x0, nx, y0, ny, bodies, walls, g):
    """Box mesh x0 + CELL*(0..nx), y0 + CELL*(0..ny); bodies: name ->
    (lo, hi, rot, friction); walls: mesh patches that are walls (20 deg)."""
    H_ = cs.HEADER
    x1, y1 = x0 + CELL * nx, y0 + CELL * ny
    sides = {"left": "xMin", "right": "xMax", "bottom": "yMin", "top": "yMax"}
    patches = [f"    {n}\n    {{\n        type patch;\n        faces ({cs.FACES[f]});\n    }}"
               for n, f in sides.items()]
    patches.append(f"    frontAndBack\n    {{\n        type empty;\n        faces "
                   f"({cs.FACES['zMin']} {cs.FACES['zMax']});\n    }}")
    cs.write(f"{path}/system/blockMeshDict", H_.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices
(
    ({x0!r} {y0!r} 0) ({x1!r} {y0!r} 0) ({x1!r} {y1!r} 0) ({x0!r} {y1!r} 0)
    ({x0!r} {y0!r} {CELL}) ({x1!r} {y0!r} {CELL}) ({x1!r} {y1!r} {CELL}) ({x0!r} {y1!r} {CELL})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1) );
boundary
(
{chr(10).join(patches)}
);
""")
    cs.write(f"{path}/system/controlDict", H_.format(cls="dictionary", obj="controlDict") + f"""
application celSiloCheck;
startFrom startTime; startTime 0; stopAt endTime; endTime 100; deltaT 1;
writeControl timeStep; writeInterval 100000000; writeFormat ascii;
gravityRamp {sc.T_RAMP!r};
damping {sc.ALPHA!r};
""")
    for name in ("fvSchemes", "fvSolution"):
        src = open(f"{path}/../ref/system/{name}").read()
        cs.write(f"{path}/system/{name}", src)
    wl = "".join(f"    {w} {{ friction {DELTA_W!r}; }}\n" for w in walls)
    cs.write(f"{path}/constant/materialProperties",
             H_.format(cls="dictionary", obj="materialProperties")
             + f"rho {sc.RHO};\n" + sc.MC + "walls\n{\n" + wl + "}\n")
    cs.write(f"{path}/constant/g", H_.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue ({g[0]!r} {g[1]!r} 0);\n")
    body = ""
    for name, (lo, hi, rot, fr) in bodies.items():
        cs.write(f"{path}/constant/triSurface/{name}.stl", box_stl(lo, hi, rot))
        body += f"    {name}\n    {{\n        surface {name}.stl;\n        friction {fr!r};\n    }}\n"
    cs.write(f"{path}/constant/rigidBodies", H_.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")
    chk = "\n".join(f"{k} {cs.foam_value(v)};" for k, v in SETTLE.items())
    cs.write(f"{path}/system/siloCheckDict",
             H_.format(cls="dictionary", obj="siloCheckDict") + chk + "\n")
    return path


def profile(path, rot):
    """Mean stress along the tank axis per 10 mm row of the tank frame:
    depth below the top, -s_aa averaged over the material volume."""
    a = np.array(sc.load(path, "cells.dat"), dtype=float)
    o = a[:, 10]
    a, o = a[o > 0], o[o > 0]
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    # tank frame coordinates of the centres and stress s' = R^T s R
    xp = c * a[:, 0] + s * a[:, 1]
    yp = -s * a[:, 0] + c * a[:, 1]
    sxx, syy, sxy = a[:, 4], a[:, 5], a[:, 7]
    syy_p = s * s * sxx - 2 * s * c * sxy + c * c * syy
    row = np.floor(yp / CELL + 1e-9).astype(int)
    rows = np.arange(int(round(H / CELL)))
    sv = np.array([-(o[row == r] * syy_p[row == r]).sum() / o[row == r].sum()
                   for r in rows])
    z = H - (rows + 0.5) * CELL
    return z, sv, (xp, yp)


def settled(s):
    ok, det = sc.settle_report(s, sc.RHO * sc.G_ACC * H)
    return ok, det


def compare(name, path, rot, ref, tol, mesh_walls):
    s = sc.summary(path)
    good, det = settled(s)
    z, sv, _ = profile(path, rot)
    zr, svr = ref
    deep = z >= W - 1e-9
    e = np.abs(sv[deep] / svr[deep] - 1.0)
    det.append(f"mean stress along the axis against the mesh-wall tank, depth >= W: "
               f"largest difference {100 * e.max():.2f} % (<= {100 * tol:g} %)")
    for frac in (0.5, 0.75, 1.0):
        i = int(np.argmin(np.abs(z - frac * H)))
        det.append(f"  depth {z[i]:.3f} m: {sv[i]:8.1f} Pa, mesh walls {svr[i]:8.1f} Pa")
    # Leak: mass in through the STL walls (lid excluded: a free top)
    mass = s["mass"]
    bodies = sorted(k[len("massIn_"):] for k in s if k.startswith("massIn_"))
    walls_in = sum(abs(s[f"massIn_{b}"]) for b in bodies if b != "lid")
    leak = walls_in / mass
    det.append(f"material through the STL walls {leak:.1e} of the total (<= 1e-6)"
               + (f"; through the lid (free top) {s['massIn_lid'] / mass:.1e}"
                  if "lid" in bodies else ""))
    # Force balance along gravity: weight per depth against the walls' push
    gx, gy = s["gx"], s["gy"]
    gm = math.hypot(gx, gy)
    up = (-gx / gm, -gy / gm)
    fb = {b: s[f"bforce_{b}_x"] * up[0] + s[f"bforce_{b}_y"] * up[1] for b in bodies}
    total = sum(fb.values()) + sum(s[f"force_{w}"] for w in mesh_walls)
    ef = abs(total - s["weight"]) / s["weight"]
    det.append("weight {:.1f} N/m = ".format(s["weight"])
               + " + ".join(f"{b} {v:.1f}" for b, v in fb.items())
               + "".join(f" + mesh {w} {s['force_' + w]:.1f}" for w in mesh_walls)
               + f": imbalance {100 * ef:.3f} % (<= 1 %)")
    ok = good and e.max() <= tol and leak <= 1e-6 and ef <= 0.01
    return sc.report(name, ok, det)


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    gm = sc.G_ACC

    # 1. Reference: WP 3.1 frictional tank with mesh walls
    ref = sc.make_case(f"{work}/ref", W, H, True, SETTLE, sc.T_RAMP, sc.ALPHA)
    sec = sc.run(ref)
    s = sc.summary(ref)
    good, det = settled(s)
    z, sv, _ = profile_ref(ref)
    ok = sc.report("Reference: mesh-wall tank settled", good,
                   det + [f"wall clock {sec:.0f} s"])

    # 2. Aligned STL side walls cutting cells at 1/4 and 3/4
    T = 0.1
    al = make_stl_case(
        f"{work}/aligned", -0.0525, 51, 0.0, 120,
        {"leftWall": ((-T, -0.05, -0.01), (0.0, H + 0.05, 0.02), 0.0, DELTA_W),
         "rightWall": ((W, -0.05, -0.01), (W + T, H + 0.05, 0.02), 0.0, DELTA_W)},
        ["bottom"], (0.0, -gm))
    sec = sc.run(al)
    ok &= compare("Aligned STL side walls (cut at 1/4 and 3/4 of a cell)",
                  al, 0.0, (z, sv), 0.02, ["bottom"])
    print(f"    wall clock {sec:.0f} s")

    # 3. The tank turned 10 deg against the mesh, gravity with it
    T = 0.3
    th = math.radians(TILT)
    c, s_ = math.cos(th), math.sin(th)
    corners = np.array([[0, 0], [W, 0], [W, H], [0, H]]) @ np.array([[c, s_], [-s_, c]])
    lo, hi = corners.min(axis=0) - 2.3 * CELL, corners.max(axis=0) + 2.3 * CELL
    nx, ny = int(math.ceil((hi[0] - lo[0]) / CELL)), int(math.ceil((hi[1] - lo[1]) / CELL))
    z0 = (-0.01, 0.02)
    ti = make_stl_case(
        f"{work}/tilted", float(lo[0]), nx, float(lo[1]), ny,
        {"leftWall": ((-T, -T, z0[0]), (0.0, H + T, z0[1]), TILT, DELTA_W),
         "rightWall": ((W, -T, z0[0]), (W + T, H + T, z0[1]), TILT, DELTA_W),
         "bottomWall": ((0.0, -T, z0[0]), (W, 0.0, z0[1]), TILT, DELTA_W),
         "lid": ((0.0, H, z0[0]), (W, H + T, z0[1]), TILT, 0.0)},
        [], (gm * s_, -gm * c))
    sec = sc.run(ti)
    ok &= compare(f"Tank turned {TILT:g} deg against the mesh (gravity turned with it)",
                  ti, TILT, (z, sv), 0.05, [])
    print(f"    wall clock {sec:.0f} s")

    print(f"\nWP 4.3a checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


def profile_ref(path):
    """Reference profile from the mesh-wall tank (cells.dat has no open
    column there: all cells open)."""
    a = np.array(sc.load(path, "cells.dat"), dtype=float)
    row = np.floor(a[:, 1] / CELL + 1e-9).astype(int)
    rows = np.arange(int(round(H / CELL)))
    sv = np.array([-a[row == r, 5].mean() for r in rows])
    return H - (rows + 0.5) * CELL, sv, None


if __name__ == "__main__":
    main()
