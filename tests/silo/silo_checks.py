"""WP 3.1 done-when checks: celFoam's cycle settles a tank under gravity with
smooth and frictional walls and matches the 2D prototype (WP 3.1p,
prototype/silo2d.py, decision 0020). The C++ scheme is decision 0030.

Run by tests/test_silo_settle.sh after celSiloCheck is built:

    python3 silo_checks.py WORKDIR

Material T0 (cohesionless Mohr-Coulomb, rho 1800, E 10 MPa, nu 0.3, phi 35,
psi 0), g = 9.81 along -y, 10 mm cells, default maxCo:

1. Lake at rest: a K0 column (K0 = nu/(1 - nu)) in discrete equilibrium
   between smooth walls, 0.1 x 0.6 m, keeps speeds below 1e-10 m/s for 1000
   steps.
2. Smooth walls, 0.1 x 1.2 m, settled from zero stress (gravity ramp 0.1 s,
   damping 250 1/s until KE <= 1e-6 strain energy): bottom stress within 2 %
   of the weight per width, s_xx / s_yy within 5 % of K0 below 10 % depth;
   then 10 % more steps undamped: stress change <= 2 % of rho g H.
3. Frictional walls and bottom (20 deg), 0.4 x 1.2 m: mean vertical stress
   within 10 % of Janssen with the K measured at the walls (depth >= W);
   measured and best-fitting K in 0.27 to 0.45 (decision 0019); wall shear /
   normal stress (cells extrapolated to the wall, as the prototype) within
   5 % of tan 20 deg where the wall slips in the settled state (at least 10
   faces per wall); weight =
   bottom + wall forces within 1 %; settled and steady as in 2.
4. Two processors (split across the height) match serial to 1e-10 of
   rho g H over 300 steps of case 3 during the ramp.

Prints PASS/FAIL per check and the prototype's values beside, exits
non-zero on any failure.
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
DELTA_W = 20.0
K0 = NU / (1.0 - NU)
T_RAMP, ALPHA, KE_TOL = 0.1, 250.0, 1e-6

MC = (f"type mohrCoulomb;\nE {E!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\npsi 0;\n"
      "sigmaT 0;\n")


def make_case(path, W, H, friction, check, ramp=0.0, damping=0.0, nProcs=1):
    nx, ny = int(round(W / CELL)), int(round(H / CELL))
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
    cs.write(f"{path}/system/controlDict", H_.format(cls="dictionary", obj="controlDict") + f"""
application celSiloCheck;
startFrom startTime; startTime 0; stopAt endTime; endTime 100; deltaT 1;
writeControl timeStep; writeInterval 100000000; writeFormat ascii;
gravityRamp {ramp!r};
damping {damping!r};
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
             + f"numberOfSubdomains {nProcs};\nmethod simple;\ncoeffs {{ n (1 {nProcs} 1); }}\n")
    fr = DELTA_W if friction else 0.0
    walls = "".join(f"    {w} {{ friction {fr!r}; }}\n" for w in ("left", "right", "bottom"))
    cs.write(f"{path}/constant/materialProperties",
             H_.format(cls="dictionary", obj="materialProperties")
             + f"rho {RHO};\n" + MC + "walls\n{\n" + walls + "}\n")
    cs.write(f"{path}/constant/g", H_.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G_ACC!r} 0);\n")
    body = "\n".join(f"{k} {cs.foam_value(v)};" for k, v in check.items())
    cs.write(f"{path}/system/siloCheckDict",
             H_.format(cls="dictionary", obj="siloCheckDict") + body + "\n")
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
        sh(["mpirun", "-np", str(nProcs), "--oversubscribe", "celSiloCheck", "-parallel"],
           "driver", env)
    else:
        sh(["celSiloCheck"], "driver")
    return time.perf_counter() - t0


def summary(path):
    out = {}
    for line in open(f"{path}/summary.dat"):
        k, v = line.split()
        out[k] = float(v)
    return out


def load(path, name, nProcs=1):
    files = [f"{path}/{name}"] if nProcs == 1 else \
        [f"{path}/{name}.proc{p}" for p in range(nProcs)]
    rows = []
    for fn in files:
        for line in open(fn):
            if not line.startswith("#"):
                rows.append(line.split())
    return rows


def cells(path, nProcs=1):
    a = np.array(load(path, "cells.dat", nProcs), dtype=float)
    # sorted by row, then x (centres carry round-off, so sort on cell indices)
    return a[np.lexsort((np.floor(a[:, 0] / CELL), np.floor(a[:, 1] / CELL)))]


def walls(path):
    rows = load(path, "walls.dat")
    out = {}
    for r in rows:
        out.setdefault(r[0], []).append([float(x) for x in r[1:]])
    return {k: np.array(v) for k, v in out.items()}


def janssen(z, W, K, mu):
    a = 2.0 * mu * K / W
    return RHO * G_ACC / a * (1.0 - np.exp(-a * z))


def grid(c, W, H):
    nx, ny = int(round(W / CELL)), int(round(H / CELL))
    return {k: c[:, i].reshape(ny, nx) for i, k in
            enumerate(("x", "y", "ux", "uy", "sxx", "syy", "szz", "sxy", "rho", "ep"))}


def report(name, ok, detail):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    for d in detail:
        print(f"    {d}")
    return ok


def settle_report(s, scale):
    ds = s["contMaxStressChange"] / scale
    good = s["settled"] == 1 and ds <= 0.02 and s["contKEratio"] <= 10 * KE_TOL
    return good, [
        f"settled after {int(s['steps'])} steps of dt = {1e6 * s['dt']:.1f} us to t = "
        f"{s['t']:.3f} s: {bool(s['settled'])}",
        f"{int(s['contSteps'])} more steps undamped: stress change {100 * ds:.3f} % of "
        f"rho g H (<= 2 %), KE / strain energy {s['contKEratio']:.1e} (<= {10 * KE_TOL:g}), "
        f"largest speed {s['contVmax']:.1e} m/s",
        f"mass change {s['massChange']:.1e} (free top: material enters with the top "
        f"cells' values as the column settles)"]


def check_lake(work):
    W, H = 0.1, 0.6
    p = make_case(f"{work}/lake", W, H, False, dict(mode="fixed", nSteps=1000, initK0=K0))
    sec = run(p)
    s = summary(p)
    ok = s["vmax"] <= 1e-10
    return report("Lake at rest: K0 column between smooth walls, 1000 steps", ok, [
        f"largest speed {s['vmax']:.1e} m/s (<= 1e-10; prototype 3e-18), dt = "
        f"{1e6 * s['dt']:.1f} us, {sec:.1f} s"])


def check_smooth(work):
    W, H = 0.1, 1.2
    p = make_case(f"{work}/smooth", W, H, False,
                  dict(mode="settle", keTol=KE_TOL, checkEvery=20, tMax=1.0,
                       continueFraction=0.1), T_RAMP, ALPHA)
    sec = run(p)
    s = summary(p)
    good, det = settle_report(s, RHO * G_ACC * H)
    g = grid(cells(p), W, H)
    eb = abs(s["force_bottom"] / s["weight"] - 1.0)
    below = (H - g["y"]) >= 0.1 * H
    kx = g["sxx"][below] / g["syy"][below]
    ek = np.abs(kx / K0 - 1.0).max()
    ok = good and eb <= 0.02 and ek <= 0.05
    return report("Smooth walls: settled column", ok, det + [
        f"bottom force {s['force_bottom']:.2f} N/m against weight {s['weight']:.2f} N/m: "
        f"{100 * eb:.3f} % (<= 2 %; prototype 0.13 %)",
        f"s_xx / s_yy below 10 % depth: {kx.min():.4f} to {kx.max():.4f}, "
        f"{100 * ek:.2f} % from K0 = {K0:.4f} (<= 5 %)",
        f"wall clock {sec:.0f} s"])


def check_friction(work):
    W, H = 0.4, 1.2
    p = make_case(f"{work}/friction", W, H, True,
                  dict(mode="settle", keTol=KE_TOL, checkEvery=20, tMax=1.0,
                       continueFraction=0.1), T_RAMP, ALPHA)
    sec = run(p)
    s = summary(p)
    good, det = settle_report(s, RHO * G_ACC * H)
    mu = math.tan(math.radians(DELTA_W))
    g = grid(cells(p), W, H)
    wl = walls(p)
    z = H - g["y"][:, 0]
    sv = -g["syy"].mean(axis=1)
    deep = z >= W
    # Horizontal stress at the walls: cells extrapolated to the wall (prototype)
    sh = 0.5 * (-(1.5 * g["sxx"][:, 0] - 0.5 * g["sxx"][:, 1])
                - (1.5 * g["sxx"][:, -1] - 0.5 * g["sxx"][:, -2]))
    K_wall = sh[deep].sum() / sv[deep].sum()
    err = lambda K: np.abs(janssen(z, W, K, mu)[deep] / sv[deep] - 1.0)  # noqa: E731
    Ks = np.linspace(0.2, 0.8, 601)
    K_fit = Ks[np.argmin([err(K).max() for K in Ks])]
    e_wall = err(K_wall).max()
    ok1 = e_wall <= 0.10 and 0.27 <= K_wall <= 0.45 and 0.27 <= K_fit <= 0.45
    det += [f"K at the walls (wall s_xx / mean s_v, depth >= W) {K_wall:.4f}; "
            f"best fit {K_fit:.4f} (error {100 * err(K_fit).max():.2f} %); prototype 0.425",
            f"mean s_v against Janssen with the measured K, depth >= W: largest error "
            f"{100 * e_wall:.2f} % (<= 10 %; prototype 1.7 %)"]
    for frac in (0.25, 0.5, 1.0):
        i = int(np.argmin(np.abs(z - frac * H)))
        det.append(f"  depth {z[i]:.3f} m: s_v {sv[i]:8.1f} Pa, Janssen "
                   f"{janssen(z[i], W, K_wall, mu):8.1f}, hydrostatic {RHO * G_ACC * z[i]:8.1f}")
    # Cell stresses extrapolated linearly to the wall (as the prototype), on
    # the faces that slip in the settled state (tangential face velocity not
    # zero), depth >= W. The prototype's zone (slipped at least 10 % of the
    # largest slip over the run) also takes in faces near the bottom that
    # slipped during the ramp and stick again: printed for information
    # (decision 0030).
    worst, nslip = 0.0, 0
    for side, (i0, i1) in (("left", (0, 1)), ("right", (-1, -2))):
        a = wl[side][np.argsort(wl[side][:, 1])]
        ext = lambda q: 1.5 * q[:, i0] - 0.5 * q[:, i1]  # noqa: E731
        ratio = np.abs(ext(g["sxy"]) / ext(g["sxx"])) / mu
        sl = a[:, 6]
        now = (np.abs(a[:, 5]) > 0) & deep
        old = (sl >= 0.1 * sl.max()) & deep
        e = np.abs(ratio[now] - 1.0).max()
        worst = max(worst, e)
        nslip = min(nslip, now.sum()) if nslip else now.sum()
        det.append(f"{side} wall slips on {int(now.sum())} faces at depth >= W in the settled "
                   f"state: shear / normal stress (cells extrapolated to the wall) / "
                   f"tan(delta_w) {ratio[now].min():.4f} to {ratio[now].max():.4f}; on the "
                   f"{int(old.sum())} faces that slipped >= 10 % of the largest slip "
                   f"(prototype's zone) {ratio[old].min():.4f} to {ratio[old].max():.4f}")
    worst = worst if nslip >= 10 else 1.0
    ok2 = worst <= 0.05
    det.append(f"wall ratio within {100 * worst:.2f} % of tan(delta_w) where it slips (<= 5 %; "
               f"prototype 1.8 % on its zone)")
    fb, fl, fr, wt = (s["force_bottom"], s["force_left"], s["force_right"], s["weight"])
    ef = abs(fb + fl + fr - wt) / wt
    ok3 = ef <= 0.01
    det.append(f"weight {wt:.1f} N/m = bottom {fb:.1f} + walls {fl:.1f} + {fr:.1f}: "
               f"imbalance {100 * ef:.3f} % (<= 1 %; prototype 0.08 %)")
    det.append(f"wall clock {sec:.0f} s")
    return report("Frictional walls and bottom: Janssen", good and ok1 and ok2 and ok3, det)


def check_parallel(work):
    W, H = 0.4, 1.2
    chk = dict(mode="fixed", nSteps=300)
    a = make_case(f"{work}/par1", W, H, True, chk, T_RAMP, ALPHA)
    b = make_case(f"{work}/par2", W, H, True, chk, T_RAMP, ALPHA, nProcs=2)
    run(a)
    run(b, 2)
    ca, cb = cells(a), cells(b, 2)
    scale = RHO * G_ACC * H
    d = np.abs(ca[:, 4:8] - cb[:, 4:8]).max() / scale
    dv = np.abs(ca[:, 2:4] - cb[:, 2:4]).max()
    ok = d <= 1e-10 and dv <= 1e-10
    return report("Two processors against serial (300 steps, frictional tank)", ok, [
        f"largest stress difference {d:.1e} of rho g H, velocity {dv:.1e} m/s (<= 1e-10)"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    ok = True
    for f in (check_lake, check_smooth, check_friction, check_parallel):
        ok &= f(work)
    print(f"\nWP 3.1 checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
