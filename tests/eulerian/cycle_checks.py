#!/usr/bin/env python3
"""WP 2.4b done-when checks: the full C++ cycle (Lagrange step, then remap).

Run by tests/test_eulerian_cycle.sh after celCycleCheck is built:

    python3 cycle_checks.py WORKDIR [--study]

Each check writes a box case (initial U, sigma in 0/), runs celCycleCheck
(tests/eulerian/cycleCheck) and compares with the 2D prototype checks of
WP 1.2 and 1.3b (prototype/lr1d_plastic.py, prototype/lr2d_plastic.py):

1. Elastic P and S sine waves along x and at 45 degrees in material
   streaming diagonally at 0.25 c_p: speed within 1 %, L1 order >= 1.7 over
   25, 50, 100 cells; step waves: no new extrema above 1 % after 10 cells of
   travel; mass and momentum to 1e-12.
2. Pre-stressed elastic disc in prescribed rigid rotation, one turn, 100
   cells across: principal stresses within 2 % inside, below 1 % outside the
   6-cell smeared edge.
3. von Mises step loading (WP 1.2) along x on a 1000 x 2 x 2 3D mesh,
   streaming at 0.25 c_p: precursor and plastic wave speeds and the HEL
   within 2 % of the analytic values.
4. 45 degree impact (von Mises, Drucker-Prager, elastic) against the same
   impact along x: beyond the along-x range by at most 1 % after the
   start-up; mass and momentum to 1e-12; 2 processors match serial.
5. Stability at the default maxCo: 45 degree (2D) and body-diagonal (3D)
   step waves with streaming, no growth over 4 crossings.

--study also scans maxCo for the stable Courant number (decision 0021).
Prints one PASS/FAIL line per check, exits non-zero on any failure.
"""

import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "lagrange"))
sys.path.insert(0, os.path.join(ROOT, "prototype"))

import cases as cs  # noqa: E402  (HEADER, FACES, MATERIALS, write, foam_value)
from lr1d import front_position  # noqa: E402

RHO, E, NU = cs.RHO, cs.E, cs.NU
K, G, M = cs.moduli()
CP = math.sqrt(M / RHO)
CS = math.sqrt(G / RHO)
A = 1e-6                    # elastic wave strain amplitude
STREAM = 0.25               # streaming speed as a fraction of c_p
SIGMA_Y = 10.0e3            # von Mises yield stress of cases.MATERIALS

failures = []


def report(name, ok, detail):
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}", flush=True)
    if not ok:
        failures.append(name)


# ----------------------------------------------------------------------------
# Cases
# ----------------------------------------------------------------------------

def centres(cells, size):
    """Cell centres in blockMesh order (x fastest, then y, then z) and the
    integer indices (i, j, k)."""
    nx, ny, nz = cells
    h = [size[d] / cells[d] for d in range(3)]
    k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    idx = np.stack([i.ravel(), j.ravel(), k.ravel()], axis=1)
    return (idx + 0.5) * np.array(h), idx


def field_list(a):
    rows = "\n".join("(" + " ".join(repr(float(x)) for x in r) + ")" for r in a)
    return f"nonuniform List<{'vector' if a.shape[1] == 3 else 'symmTensor'}>\n{len(a)}\n(\n{rows}\n)"


def make_case(path, cells, size, U, sigma, cycle, sides=None, material="linearElastic",
              maxCo=None, solveMomentum=True, nProcs=1):
    """Box with cells (nx, ny, nz) and size (Lx, Ly, Lz). sides: side name ->
    cyclic | empty | fixed; default cyclic, and empty for z when nz == 1.
    U (n, 3) and sigma (n, 6, OpenFOAM order xx xy xz yy yz zz) are the
    initial fields; cycle is cycleCheckDict."""
    sides = sides or {}
    kinds = {}
    for s in cs.FACES:
        default = "empty" if (cells[2] == 1 and s.startswith("z")) else "cyclic"
        kinds[s] = sides.get(s, default)
    patches = []
    for s, k in kinds.items():
        if k == "cyclic":
            other = s[0] + ("Max" if s.endswith("Min") else "Min")
            body = f"type cyclic;\n        neighbourPatch {other};"
        elif k == "empty":
            body = "type empty;"
        else:
            body = "type patch;"
        patches.append(f"    {s}\n    {{\n        {body}\n        faces ({cs.FACES[s]});\n    }}")
    Lx, Ly, Lz = size
    H = cs.HEADER
    cs.write(f"{path}/system/blockMeshDict", H.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices
(
    (0 0 0) ({Lx} 0 0) ({Lx} {Ly} 0) (0 {Ly} 0)
    (0 0 {Lz}) ({Lx} 0 {Lz}) ({Lx} {Ly} {Lz}) (0 {Ly} {Lz})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({cells[0]} {cells[1]} {cells[2]}) simpleGrading (1 1 1) );
boundary
(
{chr(10).join(patches)}
);
""")
    extra = f"maxCo {maxCo!r};\n" if maxCo is not None else ""
    extra += f"solveMomentum {'yes' if solveMomentum else 'no'};\n"
    cs.write(f"{path}/system/controlDict", H.format(cls="dictionary", obj="controlDict") + f"""
application celCycleCheck;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 17;
{extra}""")
    cs.write(f"{path}/system/fvSchemes", H.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default none; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    cs.write(f"{path}/system/fvSolution", H.format(cls="dictionary", obj="fvSolution"))
    cs.write(f"{path}/system/decomposeParDict", H.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {nProcs};
method simple;
coeffs {{ n ({nProcs} 1 1); }}
""")
    cs.write(f"{path}/constant/materialProperties",
             H.format(cls="dictionary", obj="materialProperties")
             + f"rho {RHO};\n" + cs.MATERIALS[material].format(E=E, NU=NU))
    bU, bS = [], []
    for s, k in kinds.items():
        tU = {"cyclic": "cyclic", "empty": "empty",
              "fixed": "fixedValue;\n        value uniform (0 0 0)"}[k]
        tS = {"cyclic": "cyclic", "empty": "empty", "fixed": "zeroGradient"}[k]
        bU.append(f"    {s}\n    {{\n        type {tU};\n    }}")
        bS.append(f"    {s}\n    {{\n        type {tS};\n    }}")
    cs.write(f"{path}/0/U", H.format(cls="volVectorField", obj="U")
             + f"dimensions [0 1 -1 0 0 0 0];\ninternalField {field_list(U)};\n"
             + "boundaryField\n{\n" + "\n".join(bU) + "\n}\n")
    cs.write(f"{path}/0/sigma", H.format(cls="volSymmTensorField", obj="sigma")
             + f"dimensions [1 -1 -2 0 0 0 0];\ninternalField {field_list(sigma)};\n"
             + "boundaryField\n{\n" + "\n".join(bS) + "\n}\n")
    body = "\n".join(f"{k} {cs.foam_value(v)};" for k, v in cycle.items())
    cs.write(f"{path}/system/cycleCheckDict",
             H.format(cls="dictionary", obj="cycleCheckDict") + body + "\n")
    return path


def run(path, nProcs=1):
    def sh(cmd, log, env=None):
        env = dict(env or os.environ, PWD=os.path.abspath(path))   # quiet cwd warnings
        with open(f"{path}/log.{log}", "w") as f:
            r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n"
                               + open(f"{path}/log.{log}").read()[-3000:])

    sh(["blockMesh"], "blockMesh")
    if nProcs > 1:
        sh(["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        sh(["mpirun", "-np", str(nProcs), "--oversubscribe", "celCycleCheck", "-parallel"],
           "driver", env)
    else:
        sh(["celCycleCheck"], "driver")
    return path


# history.dat columns
MASS, PX, PSCALE, EP = 2, 3, 6, 7
FIELDS = ("vn", "vt", "snn", "stt", "szz", "snt", "tau", "mean")


def history(path):
    return np.loadtxt(f"{path}/history.dat")


def frange(h, key):
    """(min, max) columns of a frame field over the rows of h."""
    c = 8 + 2 * FIELDS.index(key)
    return h[:, c], h[:, c + 1]


def cells_at(path, step, nProcs=1):
    if nProcs == 1:
        return np.loadtxt(f"{path}/cells_{step}.dat")
    return np.vstack([np.loadtxt(f"{path}/cells_{step}.dat.proc{p}") for p in range(nProcs)])


def cost(path):
    out = {}
    for line in open(f"{path}/log.driver"):
        m = re.match(r"(Lagrange step|Remap|Cycle) per cell per step = (\S+)", line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def conservation(h):
    dm = np.abs(h[:, MASS] - h[0, MASS]).max() / h[0, MASS]
    dp = np.abs(h[:, PX:PX + 3] - h[0, PX:PX + 3]).max() / h[:, PSCALE].max()
    return dm, dp


def steps_for(t_end, h, maxCo=0.4, c=CP, multiple=1):
    """Number of steps (a multiple of `multiple`) so that dt <= maxCo h/c."""
    n = math.ceil(t_end / (maxCo * h / c))
    return multiple * math.ceil(n / multiple)


# ----------------------------------------------------------------------------
# 1. Elastic waves in streaming material
# ----------------------------------------------------------------------------

def polarisation(mode, k):
    kv = np.array(k, float)
    kh = kv / np.linalg.norm(kv)
    if mode == "P":
        return kh
    p = np.array([0.0, 1.0, 0.0]) if k[1] == 0 else np.array([-1.0, 1.0, 0.0])
    p -= p.dot(kh) * kh
    return p / np.linalg.norm(p)


def wave_state(X, idx, n, mode, k, profile, v0):
    """Right-going simple wave (as celLagrangeStepCheck) on streaming v0."""
    kv = np.array(k, float)
    kh = kv / np.linalg.norm(kv)
    p = polarisation(mode, k)
    c = CP if mode == "P" else CS
    if profile == "sine":
        e = A * np.sin(2 * np.pi * (X @ kv))
    else:
        s = (idx @ np.array(k)) % n         # integer indices (decision 0013)
        e = np.where(2 * s < n, A, 0.0)
    U = v0[None, :] - c * e[:, None] * p[None, :]
    lam = K - 2 * G / 3
    if mode == "P":
        S = lam * np.eye(3) + 2 * G * np.outer(kh, kh)
    else:
        S = G * (np.outer(kh, p) + np.outer(p, kh))
    s6 = np.array([S[0, 0], S[0, 1], S[0, 2], S[1, 1], S[1, 2], S[2, 2]])
    return U, e[:, None] * s6[None, :]


def wave_case(work, name, n, mode, k, profile="sine", crossings=1.0, dim=2,
              maxCo=None, writeEnd=True):
    L = 1.0
    h = L / n
    cells = (n, n, n if dim == 3 else 1)
    size = (L, L, L if dim == 3 else h)
    X, idx = centres(cells, size)
    d = np.array([1.0, 1.0, 1.0 if dim == 3 else 0.0])
    v0 = STREAM * CP * d / np.linalg.norm(d)
    U, S = wave_state(X, idx, n, mode, k, profile, v0)
    c = CP if mode == "P" else CS
    t_end = crossings * L / c
    nst = steps_for(t_end, h, maxCo or 0.4)
    kh = np.array(k, float) / np.linalg.norm(k)
    cyc = dict(deltaT=t_end / nst, nSteps=nst, frameNormal=list(kh),
               streaming=list(v0), writeSteps=[nst] if writeEnd else [])
    path = make_case(os.path.join(work, name), cells, size, U, S, cyc, maxCo=maxCo)
    return run(path), v0, nst


def sine_error(path, nst, mode, k, v0):
    """L1 error of the polarised velocity against the exact sine wave carried
    by the streaming, and the phase speed of the fundamental relative to c."""
    h = history(path)
    c = cells_at(path, nst)
    t = h[-1, 1]
    kv = np.array(k, float)
    kn = np.linalg.norm(kv)
    cw = CP if mode == "P" else CS
    v = (c[:, 3:6] - v0) @ polarisation(mode, k)
    phi = 2 * np.pi * ((c[:, :3] - v0 * t) @ kv)      # material phase
    wt = 2 * np.pi * kn * cw * t
    exact = -cw * A * np.sin(phi - wt)
    err = np.abs(v - exact).mean() / (cw * A)
    a = np.sum(v * np.exp(-1j * phi))
    ph = -np.angle(a / (0.5j * cw * A))
    ph = wt + (ph - wt + np.pi) % (2 * np.pi) - np.pi
    return err, ph / (2 * np.pi * kn * t) / cw - 1


def check_waves(work):
    worst = [0.0, 0.0]
    costs = []
    for mode, k in [("P", (1, 0, 0)), ("S", (1, 0, 0)), ("P", (1, 1, 0)), ("S", (1, 1, 0))]:
        angle = "45 deg" if k[1] else "along x"
        errs, speeds = [], []
        for n in (25, 50, 100):
            p, v0, nst = wave_case(work, f"wave{mode}{k[1]}_{n}", n, mode, k)
            e, s = sine_error(p, nst, mode, k, v0)
            errs.append(e)
            speeds.append(s)
            dm, dp = conservation(history(p))
            worst = [max(worst[0], dm), max(worst[1], dp)]
            if n == 100:
                costs.append(cost(p))
        orders = [math.log2(errs[0] / errs[1]), math.log2(errs[1] / errs[2])]
        report(f"{mode}-wave {angle}, streaming: speed", max(abs(s) for s in speeds) <= 0.01,
               "speed error " + ", ".join(f"{100 * s:+.3f} %" for s in speeds)
               + " on 25, 50, 100 cells (<= 1 %)")
        report(f"{mode}-wave {angle}, streaming: order", min(orders) >= 1.7,
               f"L1 errors {errs[0]:.3e} {errs[1]:.3e} {errs[2]:.3e}, orders "
               f"{orders[0]:.2f} {orders[1]:.2f} (>= 1.7)")

        # Step wave: range of the wave-frame fields after 10 cells of travel
        n = 100
        p, v0, nst = wave_case(work, f"step{mode}{k[1]}", n, mode, k, "step", writeEnd=False)
        h = history(p)
        cw = CP if mode == "P" else CS
        hn = (1.0 / n) / np.linalg.norm(k)        # cell size along the normal
        late = cw * h[:, 1] >= 10 * hn
        key = "snn" if mode == "P" else "snt"
        jump = (M if mode == "P" else G) * A
        lo, hi = frange(h[late], key)
        over = max(hi.max() - jump, -lo.min(), 0.0) / jump
        vlo, vhi = frange(h[late], "vn" if mode == "P" else "vt")
        vj = cw * A
        vover = max(vhi.max(), -vj - vlo.min(), 0.0) / vj
        lo0, hi0 = frange(h[~late], key)
        start = max(hi0.max() - jump, -lo0.min(), 0.0) / jump
        report(f"{mode} step {angle}, streaming: extrema", max(over, vover) <= 0.01,
               f"beyond the step after 10 cells of travel: stress {100 * over:.3f} %, "
               f"velocity {100 * vover:.3f} % (<= 1 %); start-up {100 * start:.3f} %")
        dm, dp = conservation(h)
        worst = [max(worst[0], dm), max(worst[1], dp)]
    report("mass and momentum, streaming waves", max(worst) <= 1e-12,
           f"largest change: mass {worst[0]:.1e}, momentum {worst[1]:.1e} of sum(m|U|) (<= 1e-12)")
    c = {k: np.mean([x[k] for x in costs]) for k in costs[0]}
    print(f"      cost (2D 100 x 100, elastic, serial): Lagrange {c['Lagrange step']:.2e}, "
          f"remap {c['Remap']:.2e}, cycle {c['Cycle']:.2e} s per cell per step")


# ----------------------------------------------------------------------------
# 2. Pre-stressed disc in rigid rotation
# ----------------------------------------------------------------------------

def rotated(s0, angle):
    """OpenFOAM-order stress s0 (6,) rotated in plane by angle (rad)."""
    c, sn = math.cos(angle), math.sin(angle)
    Q = np.array([[c, -sn, 0], [sn, c, 0], [0, 0, 1]])
    T = np.array([[s0[0], s0[1], s0[2]], [s0[1], s0[3], s0[4]], [s0[2], s0[4], s0[5]]])
    R = Q @ T @ Q.T
    return np.array([R[0, 0], R[0, 1], R[0, 2], R[1, 1], R[1, 2], R[2, 2]])


def principal_2d(s):
    sxx, sxy, syy = s[:, 0], s[:, 1], s[:, 3]
    m = 0.5 * (sxx + syy)
    r = np.hypot(0.5 * (sxx - syy), sxy)
    return m + r, m - r, 0.5 * np.arctan2(2 * sxy, sxx - syy)


def check_disc(work, cells=100, margin=8, band=6, steps=2000):
    """WP 1.3b disc (lr2d_plastic.disc_case) in the full cycle: radius 0.5 m,
    omega 1 rad/s, prescribed motion (solveMomentum no, U held), fixed-velocity
    box walls, `steps` cycles per turn."""
    R = 0.5
    dx = 2 * R / cells
    n = cells + 2 * margin
    L = n * dx
    xc = 0.5 * L
    X, _ = centres((n, n, 1), (L, L, dx))
    sub = 8
    frac = np.zeros(len(X))
    for a in range(sub):
        for b in range(sub):
            x = X[:, 0] - 0.5 * dx + (a + 0.5) / sub * dx
            y = X[:, 1] - 0.5 * dx + (b + 0.5) / sub * dx
            frac += (x - xc) ** 2 + (y - xc) ** 2 < R * R
    frac /= sub * sub
    s1, s2, th = -40.0e3, -100.0e3, math.radians(30.0)
    s0 = rotated(np.array([s1, 0.0, 0.0, s2, 0.0, -60.0e3]), th)
    S = frac[:, None] * s0[None, :]
    U = np.zeros((len(X), 3))
    w = 1.0
    cyc = dict(deltaT=2 * math.pi / w / steps, nSteps=steps, holdU="yes",
               prescribe="rotation", omega=[0.0, 0.0, w], centre=[xc, xc, 0.0],
               writeSteps=[steps // 4, steps // 2, steps])
    sides = {s: "fixed" for s in ("xMin", "xMax", "yMin", "yMax")}
    p = run(make_case(os.path.join(work, "disc"), (n, n, 1), (L, L, dx), U, S, cyc,
                      sides=sides, solveMomentum=False))
    r = np.hypot(X[:, 0] - xc, X[:, 1] - xc)
    inner = r < R - band * dx
    outer = r > R + band * dx
    scale = abs(s2)
    worst = dict(princ=0.0, tensor=0.0, ang=0.0, out=0.0)
    for q, st in ((0.25, steps // 4), (0.5, steps // 2), (1.0, steps)):
        c = cells_at(p, st)
        s = c[:, 6:12]
        ex = rotated(s0, 2 * math.pi * q)
        p1, p2, ang = principal_2d(s[inner])
        e_p = max(np.abs(p1 - s1).max(), np.abs(p2 - s2).max()) / scale
        e_t = np.abs(s[inner] - ex).max() / scale
        da = (ang - (th + 2 * math.pi * q) + 0.5 * math.pi) % math.pi - 0.5 * math.pi
        e_a = math.degrees(np.abs(da).max())
        e_o = np.abs(s[outer]).max() / scale
        row = (X[:, 1] > xc) & (X[:, 1] < xc + dx)
        prof = np.abs(c[row, 11]) / 60.0e3
        width = int(np.sum((prof > 0.1) & (prof < 0.9)))
        print(f"      {q:.2f} turn: principal stresses {100 * e_p:.3f} %, tensor "
              f"{100 * e_t:.3f} %, direction {e_a:.3f} deg; outside {100 * e_o:.3f} %; "
              f"edge cells (10-90 %, both sides) {width}")
        worst = dict(princ=max(worst["princ"], e_p), tensor=max(worst["tensor"], e_t),
                     ang=max(worst["ang"], e_a), out=max(worst["out"], e_o))
    report("disc in rigid rotation, principal stresses", worst["princ"] <= 0.02,
           f"largest error inside {100 * worst['princ']:.3f} % of |s2| (<= 2 %; "
           f"prototype 0.03 %), tensor {100 * worst['tensor']:.3f} %, direction "
           f"{worst['ang']:.3f} deg")
    report("disc in rigid rotation, nothing beyond the edge", worst["out"] <= 0.01,
           f"largest |sigma| more than {band} cells outside {100 * worst['out']:.3f} % (<= 1 %)")
    dm, _ = conservation(history(p))
    print(f"      mass change {dm:.1e} (open walls: reported only)")


# ----------------------------------------------------------------------------
# 3. 1D plastic step wave on a 3D mesh
# ----------------------------------------------------------------------------

def check_plastic_1d(work, n=1000):
    """WP 1.2 step loading (lr1d_plastic.step_loading) along x on n x 2 x 2."""
    h = 1.0 / n
    hel = SIGMA_Y * M / (2 * G)
    c_p = math.sqrt(K / RHO)
    v0 = STREAM * CP
    v1 = hel / (RHO * CP)
    V = 3 * v1
    sig_final = -(hel + RHO * c_p * (V - v1))
    cells, size = (n, 2, 2), (1.0, 2 * h, 2 * h)
    X, _ = centres(cells, size)
    U = np.zeros((len(X), 3))
    U[:, 0] = np.where(X[:, 0] < 0.5, v0 + V, v0 - V)
    S = np.zeros((len(X), 6))
    t1, t2 = 0.08 / CP, 0.22 / CP
    nst = steps_for(t2, h, multiple=11)          # t1 / t2 = 4 / 11
    n1 = nst * 4 // 11
    cyc = dict(deltaT=t2 / nst, nSteps=nst, writeSteps=[n1, nst],
               streaming=[v0, 0.0, 0.0])
    p = run(make_case(os.path.join(work, "plastic1d"), cells, size, U, S, cyc,
                      material="vonMises"))

    def profile(st):
        c = cells_at(p, st)
        o = np.lexsort((c[:, 2], c[:, 1], c[:, 0]))
        c = c[o].reshape(n, 4, -1)
        spread = np.abs(c[:, :, 6] - c[:, :1, 6]).max()
        return c[:, 0, 0], -c[:, :, 6].mean(axis=1), spread

    def fronts(st):
        t = st * t2 / nst
        xc, sig, sp = profile(st)
        x_imp = 0.5 + v0 * t
        xe = front_position(xc, sig, 0.5 * hel, x_imp + (CP - V) * t, 1.0)
        xp = front_position(xc, sig, hel + 0.5 * (-sig_final - hel), x_imp + c_p * t, 1.0)
        return xe, xp, xc, sig, sp

    xe1, xp1, _, _, sp1 = fronts(n1)
    xe2, xp2, xc, sig, sp2 = fronts(nst)
    u_e = (xe2 - xe1) / (t2 - t1) - (v0 - V)
    u_p = (xp2 - xp1) / (t2 - t1) - (v0 - V + v1)
    rel = (xc - xp2) % 1.0
    gap = (xe2 - xp2) % 1.0
    amp = np.median(sig[(rel > gap / 3) & (rel < 2 * gap / 3)])
    e_e, e_p, e_a = u_e / CP - 1, u_p / c_p - 1, amp / hel - 1
    report("plastic step (3D mesh), precursor speed", abs(e_e) <= 0.02,
           f"{u_e:.3f} m/s vs c_p {CP:.3f} ({100 * e_e:+.3f} %, <= 2 %)")
    report("plastic step (3D mesh), plastic wave speed", abs(e_p) <= 0.02,
           f"{u_p:.3f} m/s vs sqrt(K/rho) {c_p:.3f} ({100 * e_p:+.3f} %, <= 2 %)")
    report("plastic step (3D mesh), HEL", abs(e_a) <= 0.02,
           f"precursor {amp:.1f} Pa vs Y M/(2G) {hel:.1f} Pa ({100 * e_a:+.3f} %, <= 2 %)")
    h_ = history(p)
    dm, dp = conservation(h_)
    report("plastic step (3D mesh), stays 1D and conserves",
           max(sp1, sp2) <= 1e-9 * abs(sig_final) and max(dm, dp) <= 1e-12,
           f"spread across y, z {max(sp1, sp2):.1e} Pa; mass {dm:.1e}, momentum {dp:.1e} "
           f"(<= 1e-12)")


# ----------------------------------------------------------------------------
# 4. 45 degree elastic-plastic impact
# ----------------------------------------------------------------------------

def impact_case(work, name, angle, n, material, nProcs=1):
    """lr2d_plastic.impact: s < period/2 moves at +V along the normal, the
    rest at -V, plus diagonal streaming at 0.25 c_p; travel 0.22 period."""
    hel = SIGMA_Y * M / (2 * G)
    V = 3.0 * hel / (RHO * CP)
    if angle == 45:
        cells, size = (n, n, 1), (1.0, 1.0, 1.0 / n)
        nvec = np.array([1.0, 1.0, 0.0]) / math.sqrt(2)
        period = math.sqrt(0.5)
        X, idx = centres(cells, size)
        s = ((idx[:, 0] + idx[:, 1] + 1) % n) * (1.0 / n) * period
        hn = 1.0 / n / math.sqrt(2)
    else:
        cells, size = (n, 4, 1), (1.0, 4.0 / n, 1.0 / n)
        nvec = np.array([1.0, 0.0, 0.0])
        period = 1.0
        X, idx = centres(cells, size)
        s = X[:, 0]
        hn = 1.0 / n
    sign = np.where(s < 0.5 * period, 1.0, -1.0)
    vd = STREAM * CP / math.sqrt(2)
    stream = np.array([vd, vd, 0.0])
    U = stream[None, :] + V * sign[:, None] * nvec[None, :]
    S = np.zeros((len(X), 6))
    t_end = 0.22 * period / CP
    nst = steps_for(t_end, 1.0 / n)
    cyc = dict(deltaT=t_end / nst, nSteps=nst, frameNormal=list(nvec),
               streaming=list(stream), writeSteps=[nst])
    p = run(make_case(os.path.join(work, name), cells, size, U, S, cyc,
                      material=material, nProcs=nProcs), nProcs)
    return p, V, hn, nst


def check_impact(work, n=100):
    for material in ("vonMises", "druckerPrager", "linearElastic"):
        p0, V, hn0, _ = impact_case(work, f"impact0_{material}", 0, n, material)
        p, _, hn, nst = impact_case(work, f"impact45_{material}", 45, n, material)
        h0, h = history(p0), history(p)
        rng0 = {k: (frange(h0, k)[0].min(), frange(h0, k)[1].max()) for k in FIELDS}
        sscale = rng0["snn"][1] - rng0["snn"][0]
        late = CP * h[:, 1] >= 10 * hn
        e = {"start": 0.0, "main": 0.0}
        worst_field = ""
        for win, sel in (("start", ~late), ("main", late)):
            for k in FIELDS:
                lo, hi = frange(h[sel], k)
                scale = 2 * V if k.startswith("v") else sscale
                x = max(hi.max() - rng0[k][1], rng0[k][0] - lo.min(), 0.0) / scale
                if win == "main" and x > e["main"]:
                    worst_field = k
                e[win] = max(e[win], x)
        dm, dp = conservation(h)
        report(f"45 deg impact, {material}", e["main"] <= 0.01 and max(dm, dp) <= 1e-12,
               f"beyond the along-x range after the start-up {100 * e['main']:.3f} % "
               f"({worst_field or '-'}, <= 1 %), start-up {100 * e['start']:.3f} % "
               f"(reported); mass {dm:.0e}, momentum {dp:.0e} (<= 1e-12)")
        if material == "vonMises":
            pp, _, _, _ = impact_case(work, "impact45_vonMises_par", 45, n, material, nProcs=2)
            a = cells_at(p, nst)
            b = cells_at(pp, nst, nProcs=2)
            def order(c):       # by cell index (centres differ by round-off)
                return c[np.lexsort((np.floor(c[:, 1] * n), np.floor(c[:, 0] * n)))]
            a, b = order(a), order(b)
            dU = np.abs(a[:, 3:6] - b[:, 3:6]).max() / V
            dS = np.abs(a[:, 6:12] - b[:, 6:12]).max() / sscale
            report("45 deg impact, 2 processors vs serial", max(dU, dS) <= 1e-12,
                   f"largest difference U {dU:.1e} of V, sigma {dS:.1e} (<= 1e-12)")


# ----------------------------------------------------------------------------
# 5. Stability at the default maxCo
# ----------------------------------------------------------------------------

def growth(work, name, dim, maxCo, crossings=4.0, n=None):
    """45 degree (2D, n 50) or body-diagonal (3D, n 20) step P-wave with
    streaming: largest |U - v0| / (c_p e) in each crossing. Stable: finite and
    never above 1.05 (the stable scheme overshoots 0.5 % in 2D and 2.2 % in
    3D, and the maximum moves by about 0.5 % from crossing to crossing as the
    streaming shifts the front against the mesh, so a strict "no growth"
    comparison of two crossings is noise; an unstable mode passes 1.05 within
    a few crossings)."""
    n = n or (50 if dim == 2 else 20)
    k = (1, 1, 0) if dim == 2 else (1, 1, 1)
    p, _, nst = wave_case(work, name, n, "P", k, "step", crossings=crossings, dim=dim,
                          maxCo=maxCo, writeEnd=False)
    h = history(p)
    lo, hi = frange(h, "vn")
    vmax = np.maximum(np.abs(lo), np.abs(hi)) / (CP * A)
    per = (len(h) - 1) // int(crossings)
    each = [float(vmax[1 + i * per:1 + (i + 1) * per].max()) for i in range(int(crossings))]
    ok = bool(np.all(np.isfinite(h[:, 2:]))) and max(each) <= 1.05
    return ok, each


def check_stability(work, maxCo):
    for dim in (2, 3):
        ok, each = growth(work, f"stab{dim}d", dim, maxCo)
        what = "45 degree (2D)" if dim == 2 else "body-diagonal (3D)"
        report(f"stable at maxCo {maxCo}, {what} step with streaming, 4 crossings", ok,
               "largest |U - v0| / (c_p e) per crossing "
               + " ".join(f"{x:.4f}" for x in each) + " (<= 1.05)")


def study(work):
    """Stable Courant number: unstable if a step P-wave with streaming passes
    1.5 times its amplitude (or overflows) within 16 crossings. The bounded
    slow rise of the overshoot over many crossings (2D: 0.5 % to 5.7 % over
    40 crossings at maxCo 0.4, also without streaming) is not counted."""
    print("Stable Courant number of the full cycle (step P-wave, streaming 0.25 c_p,")
    print("16 crossings; largest |U - v0| / (c_p e) per crossing, first and last):")
    for dim, cos in ((2, (0.6, 0.7, 0.75, 0.8, 0.9)), (3, (0.5, 0.55, 0.6, 0.65, 0.7))):
        for co in cos:
            try:
                _, each = growth(work, f"study{dim}d_{co}", dim, co, crossings=16.0)
                ok = max(each) < 1.5
                print(f"  {dim}D maxCo {co}: first {each[0]:.4f}, last {each[-1]:.4f}, "
                      f"largest {max(each):.4f} -> {'stable' if ok else 'UNSTABLE'}", flush=True)
            except (RuntimeError, ValueError):
                print(f"  {dim}D maxCo {co}: run failed -> UNSTABLE", flush=True)


def main():
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    if "--study" in sys.argv:
        study(work)
        return 0
    maxCo = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
    print("[1] Elastic waves in material streaming at 0.25 c_p (diagonal)")
    check_waves(work)
    print("[2] Pre-stressed disc in rigid rotation, one turn")
    check_disc(work)
    print("[3] von Mises step loading along x on a 3D mesh (WP 1.2)")
    check_plastic_1d(work)
    print("[4] 45 degree elastic-plastic impact (WP 1.3b)")
    check_impact(work)
    print(f"[5] Stability at the default maxCo {maxCo}")
    check_stability(work, maxCo)
    if failures:
        print(f"\n{len(failures)} check(s) failed: " + "; ".join(failures))
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
