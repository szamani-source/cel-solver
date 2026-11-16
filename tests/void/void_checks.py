"""WP 3.2a done-when checks: one material plus void in celFoam's cycle
(materialProperties freeSurface yes, decision 0031) against the prototype
cases of WP 3.2p (prototype/void2d.py --verify, decision 0021):

1. Lake at rest: K0 column between smooth walls in a closed 0.1 x 0.6 m box,
   filled to 0.455 m (top cell half full), 10 mm cells, full gravity, no
   damping, 1000 steps: speeds <= 1e-8 m/s, no stress in empty cells, mass
   to 1e-10.
2. Elastic block of 10 x 10 cells translating through void at 5 m/s for 50
   cells, no gravity: mass to 1e-10, |stress| <= 1e-6 Pa, material at most 3
   cells beyond the exact block, centroid travel within 0.1 %.
3. Elastic block (10 x 10 cells) falling 0.1 m in void onto a rigid floor:
   acceleration within 0.5 % of g, |stress| <= 1e-6 Pa while falling; lowest
   centroid >= 0.04 m (half the block less a cell), peak floor force above
   the weight, mass to 1e-10.
4. Cohesionless block (T0) pulled apart at +-0.5 m/s, no gravity: void
   opened inside within 5 % of 2 v t, at least one empty column, no tension
   (normal stress <= 1e-6 Pa), momentum and mass to 1e-10.
5. Two processors match serial to 1e-10 (stresses scaled by rho c v): case 2
   with 25 cells of travel across the split, and case 4 with the split at
   the crack.

Run by tests/test_void.sh after celVoidCheck is built:

    python3 void_checks.py WORKDIR
"""

import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lagrange"))
import cases as cs  # noqa: E402  (HEADER, FACES, write, foam_value)

RHO, E, NU, PHI = 1800.0, 10.0e6, 0.3, 35.0
G_ACC = 9.81
CELL = 0.01
MC = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
      "sigmaT 0;\n")
ELASTIC = f"type linearElastic;\nE {E!r};\nnu {NU!r};\n"


def make_case(path, nx, ny, material, check, gravity=True, nProcs=1):
    W, H = nx * CELL, ny * CELL
    H_ = cs.HEADER
    sides = {"left": "xMin", "right": "xMax", "bottom": "yMin", "top": "yMax"}
    patches = [f"    {n}\n    {{\n        type patch;\n        faces ({cs.FACES[f]});\n    }}"
               for n, f in sides.items()]
    patches.append(f"    frontAndBack\n    {{\n        type empty;\n        faces "
                   f"({cs.FACES['zMin']} {cs.FACES['zMax']});\n    }}")
    cs.write(f"{path}/system/blockMeshDict", H_.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices
(
    (0 0 0) ({W} 0 0) ({W} {H} 0) (0 {H} 0)
    (0 0 {CELL}) ({W} 0 {CELL}) ({W} {H} {CELL}) (0 {H} {CELL})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1) );
boundary
(
{chr(10).join(patches)}
);
""")
    cs.write(f"{path}/system/controlDict", H_.format(cls="dictionary", obj="controlDict") + """
application celVoidCheck;
startFrom startTime; startTime 0; stopAt endTime; endTime 100; deltaT 1;
writeControl timeStep; writeInterval 100000000; writeFormat ascii;
""")
    cs.write(f"{path}/system/fvSchemes", H_.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default none; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    cs.write(f"{path}/system/fvSolution", H_.format(cls="dictionary", obj="fvSolution"))
    cs.write(f"{path}/system/decomposeParDict", H_.format(cls="dictionary", obj="decomposeParDict")
             + f"numberOfSubdomains {nProcs};\nmethod simple;\ncoeffs {{ n ({nProcs} 1 1); }}\n")
    walls = "".join(f"    {w} {{ friction 0; }}\n" for w in sides)
    cs.write(f"{path}/constant/materialProperties",
             H_.format(cls="dictionary", obj="materialProperties")
             + f"rho {RHO};\nfreeSurface yes;\n" + material
             + os.environ.get("CEL_EXTRA_MATERIAL", "") + "walls\n{\n" + walls + "}\n")
    if gravity:
        cs.write(f"{path}/constant/g", H_.format(cls="uniformDimensionedVectorField", obj="g")
                 + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G_ACC!r} 0);\n")
    body = "\n".join(f"{k} {cs.foam_value(v)};" for k, v in check.items())
    cs.write(f"{path}/system/voidCheckDict",
             H_.format(cls="dictionary", obj="voidCheckDict") + body + "\n")
    return path


def run(path, nProcs=1):
    def sh(cmd, log, env=None):
        env = dict(env or os.environ, PWD=os.path.abspath(path))
        with open(f"{path}/log.{log}", "w") as f:
            r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n"
                               + open(f"{path}/log.{log}").read()[-3000:])
    t0 = time.perf_counter()
    sh(["blockMesh"], "blockMesh")
    if nProcs > 1:
        sh(["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        sh(["mpirun", "-np", str(nProcs), "--oversubscribe", "celVoidCheck", "-parallel"],
           "driver", env)
    else:
        sh(["celVoidCheck"], "driver")
    return time.perf_counter() - t0


def history(path):
    return np.loadtxt(f"{path}/history.dat")


def cells(path, nProcs=1):
    files = [f"{path}/cells.dat"] if nProcs == 1 else \
        [f"{path}/cells.dat.proc{p}" for p in range(nProcs)]
    a = np.vstack([np.loadtxt(f) for f in files])
    return a[np.lexsort((np.floor(a[:, 0] / CELL), np.floor(a[:, 1] / CELL)))]


def report(name, ok, detail):
    print(f"  {detail}: {'PASS' if ok else 'FAIL'}")
    return ok


# history columns
STEP, T, MASS, PX, PY, CX, CY, SMAX, SEMPTY, UMAX, FY, SNMAX = range(12)


def check_lake(work):
    print("[1] Lake at rest: K0 column between smooth walls, 0.1 x 0.6 m box filled to")
    print("    0.455 m (top cell half full), 10 mm cells, full gravity, 1000 steps")
    p = make_case(f"{work}/lake", 10, 60, MC,
                  {"fill": [-1.0, -1.0, 1.0, 0.455], "U0": [0, 0, 0],
                   "initK0": NU / (1 - NU), "nSteps": 1000})
    t = run(p)
    h = history(p)
    c = cells(p)
    vmax = h[1:, UMAX].max()
    dm = abs(h[-1, MASS] / h[0, MASS] - 1)
    part = c[(c[:, 2] > 1e-12) & (c[:, 2] < 1 - 1e-12), 2]
    ok = vmax <= 1e-8 and h[:, SEMPTY].max() == 0.0 and dm <= 1e-10
    return report("lake", ok,
                  f"largest speed {vmax:.1e} m/s (<= 1e-8; prototype 2.5e-16); stress in "
                  f"empty cells {h[:, SEMPTY].max():.1e} Pa (0); partly filled cells al = "
                  f"{part.mean() if part.size else float('nan'):.3f}; mass change {dm:.0e} "
                  f"(<= 1e-10); {t:.0f} s")


def translating(work, name, travel=50, nProcs=1):
    v = 5.0
    p = make_case(f"{work}/{name}", 80, 20, ELASTIC,
                  {"fill": [0.05, 0.05, 0.15, 0.15], "U0": [v, 0, 0],
                   "tEnd": travel * CELL / v}, gravity=False, nProcs=nProcs)
    t = run(p, nProcs)
    return p, t


def check_translate(work):
    print("\n[2] Block of 10 x 10 cells translating through void at 5 m/s for 50 cells")
    print("    (elastic, no gravity, 10 mm cells)")
    p, t = translating(work, "translate")
    h = history(p)
    c = cells(p)
    v, tend = 5.0, h[-1, T]
    dm = abs(h[-1, MASS] / h[0, MASS] - 1)
    smax = h[:, SMAX].max()
    xa, xb, y0, y1 = 0.05 + v * tend, 0.15 + v * tend, 0.05, 0.15
    x, y, al = c[:, 0], c[:, 1], c[:, 2]
    has = al > 1e-12
    dist = np.maximum(np.maximum(np.maximum(xa - x, x - xb), 0),
                      np.maximum(np.maximum(y0 - y, y - y1), 0)) / CELL
    spread = dist[has].max() + 0.5
    rows = np.floor(y / CELL).astype(int)
    partial = max(((al > 1e-6) & (al < 1 - 1e-6) & (rows == r)).sum() for r in set(rows))
    dcx = (h[-1, CX] - h[0, CX]) / (v * tend) - 1
    ok = dm <= 1e-10 and smax <= 1e-6 and spread <= 3.0 and abs(dcx) <= 1e-3
    return report("translate", ok,
                  f"{int(h[-1, STEP])} steps; mass change {dm:.0e} (<= 1e-10); largest "
                  f"|stress| {smax:.1e} Pa (<= 1e-6); material at most {spread:.2f} cells "
                  f"beyond the exact block (<= 3; prototype 0.5); partly filled cells per row "
                  f"{partial}; centroid travel error {100 * dcx:.3f} % (<= 0.1 %); "
                  f"vx {c[has, 4].min():.6f} to {c[has, 4].max():.6f} m/s; {t:.0f} s")


def check_fall(work):
    print("\n[3] Elastic block (10 x 10 cells) falling 0.1 m in void onto a rigid floor")
    y0 = 0.10
    t_land = math.sqrt(2 * y0 / G_ACC)
    p = make_case(f"{work}/fall", 20, 30, ELASTIC,
                  {"fill": [0.05, y0, 0.15, y0 + 0.1], "U0": [0, 0, 0],
                   "tEnd": t_land + 0.03})
    t = run(p)
    h = history(p)
    fall = (h[:, T] > 0) & (h[:, T] < 0.9 * t_land)
    acc = -2 * np.polyfit(h[fall, T], h[fall, CY], 2)[0]
    e_acc = abs(acc / G_ACC - 1)
    s_fall = h[fall, SMAX].max()
    after = h[:, T] >= 0.9 * t_land
    ymin = h[after, CY].min()
    weight = h[0, MASS] * G_ACC / CELL     # per m depth (mesh 10 mm thick)
    fb = h[after, FY].max() / CELL
    dm = abs(h[-1, MASS] / h[0, MASS] - 1)
    ok1 = report("fall", e_acc <= 0.005 and s_fall <= 1e-6,
                 f"free fall: centroid acceleration {acc:.5f} m/s^2, {100 * e_acc:.3f} % from "
                 f"g (<= 0.5 %); largest |stress| while falling {s_fall:.1e} Pa (<= 1e-6)")
    ok2 = report("land", dm <= 1e-10 and ymin >= 0.04 and fb > weight,
                 f"landing at {math.sqrt(2 * G_ACC * y0):.2f} m/s: lowest centroid "
                 f"{ymin:.4f} m (>= 0.04); peak floor force {fb:.0f} N/m (weight "
                 f"{weight:.0f} N/m); mass change {dm:.0e} (<= 1e-10); {t:.0f} s")
    return ok1 and ok2


def separating(work, name, nProcs=1):
    v, tend = 0.5, 0.02
    p = make_case(f"{work}/{name}", 60, 10, MC,
                  {"fill": [0.2, -1.0, 0.4, 1.0], "U0": [v, 0, 0], "xSplit": 0.3,
                   "Uleft": [-v, 0, 0], "tEnd": tend}, gravity=False, nProcs=nProcs)
    return p, run(p, nProcs)


def check_separation(work):
    print("\n[4] Separation: cohesionless block (T0) pulled apart at +-0.5 m/s, no gravity")
    v = 0.5
    p, t = separating(work, "separate")
    h = history(p)
    c = cells(p)
    x, al = c[:, 0], c[:, 2]
    m0 = h[0, MASS]
    dm = abs(h[-1, MASS] / m0 - 1)
    dp = abs(h[-1, PX] - h[0, PX]) / (m0 * v)
    sn = h[:, SNMAX].max()
    cols = np.unique(np.floor(x / CELL).astype(int))
    gap = sum(1 for k in cols if 0.2 < (k + 0.5) * CELL < 0.4
              and (al[np.floor(x / CELL).astype(int) == k] < 1e-6).all())
    inner = (x > 0.25) & (x < 0.35)
    void_cells = (1 - al[inner]).sum() / 10.0
    expect = 2 * v * h[-1, T] / CELL
    e_gap = abs(void_cells / expect - 1)
    ok = dm <= 1e-10 and sn <= 1e-6 and gap >= 1 and dp <= 1e-10 and e_gap <= 0.05
    return report("separate", ok,
                  f"void opened inside the block {void_cells:.3f} cells per row (2 v t = "
                  f"{expect:.3f}, {100 * e_gap:.1f} % off, <= 5 %), {gap} empty column(s); "
                  f"largest normal stress {sn:.1e} Pa (<= 1e-6); momentum change {dp:.0e}; "
                  f"mass change {dm:.0e}; {t:.0f} s")


def compare(ps, pp, v):
    a, b = cells(ps), cells(pp, 2)
    cols = {"alpha": 2, "rho": 3, "Ux": 4, "Uy": 5, "s_xx": 6, "s_yy": 7, "s_xy": 8}
    # stresses against the impact stress rho c v of the block's speed
    M = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    sv = RHO * math.sqrt(M / RHO) * v
    scale = {"alpha": 1.0, "rho": RHO, "Ux": v, "Uy": v, "s_xx": sv, "s_yy": sv,
             "s_xy": sv}
    if a.shape != b.shape:
        return False, "cell counts differ"
    diffs = {k: np.abs(a[:, i] - b[:, i]).max() / scale[k] for k, i in cols.items()}
    return max(diffs.values()) <= 1e-10, \
        "largest difference " + ", ".join(f"{k} {d:.0e}" for k, d in diffs.items())


def check_parallel(work):
    print("\n[5] Two processors against serial (scaled differences <= 1e-10): translating")
    print("    block, 25 cells of travel across the split at x = 0.4 m; cohesionless block")
    print("    pulled apart exactly at the split (x = 0.3 m)")
    ok1, d1 = compare(translating(work, "ser", travel=25)[0],
                      translating(work, "par", travel=25, nProcs=2)[0], 5.0)
    ok2, d2 = compare(separating(work, "sepSer")[0], separating(work, "sepPar", 2)[0], 0.5)
    return report("parallel", ok1 and ok2, f"translating: {d1}; separating: {d2}")


def main():
    work = sys.argv[1]
    only = sys.argv[2:] or ["lake", "translate", "fall", "separate", "parallel"]
    checks = {"lake": check_lake, "translate": check_translate, "fall": check_fall,
              "separate": check_separation, "parallel": check_parallel}
    t0 = time.perf_counter()
    ok = True
    for k in only:
        ok &= checks[k](work)
    print(f"\nWP 3.2a checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
