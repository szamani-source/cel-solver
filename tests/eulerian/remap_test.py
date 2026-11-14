#!/usr/bin/env python3
"""WP 2.4 remap checks (run by tests/test_remap.sh; needs OpenFOAM loaded and
celRemapCheck built).

Each case: a blockMesh box written here, initial fields written here as cell
averages (Gauss quadrature, or exact for aligned squares), a prescribed
motion in system/remapCheckDict, then celRemapCheck, which calls
eulerianSolid::remap alone. The checks are the "done when" list of
docs/wp/WP-2.4.md; the numbers are printed and asserted below.

Usage: remap_test.py [--quick] [--study]
  --study   also scan the material Courant number for stability (report only)
"""

import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "prototype"))
import lr1d                     # noqa: E402
import lr1d_plastic             # noqa: E402

DEFAULT_CO = 0.25               # material Courant number per direction
FAILS = []


def check(what, ok):
    print(("PASS  " if ok else "FAIL  ") + what, flush=True)
    if not ok:
        FAILS.append(what)


# ----------------------------------------------------------------------------
# Case writing
# ----------------------------------------------------------------------------

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def make_case(work, name, n, L=1.0, periodic=True, material=""):
    """Box of cubic cells of size h = L/n[0]: [0,L] x [0,n[1] h] x [0,n[2] h]
    (n[2] == 1: a 2D case with empty front and back). periodic: cyclic pairs,
    else walls."""
    case = os.path.join(work, name)
    if os.path.isdir(case):
        shutil.rmtree(case)
    nx, ny, nz = n
    two_d = nz == 1
    h = L / nx
    Ly, Lz = ny * h, nz * h
    v = [(0, 0, 0), (L, 0, 0), (L, Ly, 0), (0, Ly, 0),
         (0, 0, Lz), (L, 0, Lz), (L, Ly, Lz), (0, Ly, Lz)]
    faces = {
        "x0": "(0 4 7 3)", "x1": "(1 2 6 5)",
        "y0": "(0 1 5 4)", "y1": "(3 7 6 2)",
        "z0": "(0 3 2 1)", "z1": "(4 5 6 7)",
    }
    patches = []
    for a in "xyz":
        if a == "z" and two_d:
            patches.append(("frontBack", "empty",
                            [faces["z0"], faces["z1"]], None))
            continue
        if periodic:
            patches.append((a + "0", "cyclic", [faces[a + "0"]], a + "1"))
            patches.append((a + "1", "cyclic", [faces[a + "1"]], a + "0"))
        else:
            patches.append((a + "walls", "wall",
                            [faces[a + "0"], faces[a + "1"]], None))
    bnd = ""
    for pname, ptype, fl, nbr in patches:
        bnd += f"    {pname}\n    {{\n        type {ptype};\n"
        if nbr:
            bnd += f"        neighbourPatch {nbr};\n"
        bnd += "        faces (" + " ".join(fl) + ");\n    }\n"
    write(os.path.join(case, "system/blockMeshDict"),
          HEADER.format(cls="dictionary", obj="blockMeshDict")
          + "scale 1;\nvertices (" + " ".join("(%r %r %r)" % p for p in v) + ");\n"
          + f"blocks ( hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1) );\n"
          + "boundary\n(\n" + bnd + ");\n")
    write(os.path.join(case, "system/controlDict"),
          HEADER.format(cls="dictionary", obj="controlDict")
          + "application celRemapCheck;\nstartFrom startTime;\nstartTime 0;\n"
          "stopAt endTime;\nendTime 1;\ndeltaT 1;\nwriteControl timeStep;\n"
          "writeInterval 1;\nwriteFormat ascii;\nwritePrecision 17;\n"
          "timeFormat general;\ntimePrecision 12;\nrunTimeModifiable false;\n")
    write(os.path.join(case, "system/fvSchemes"),
          HEADER.format(cls="dictionary", obj="fvSchemes")
          + "ddtSchemes{default Euler;}\ngradSchemes{default Gauss linear;}\n"
          "divSchemes{default none;}\nlaplacianSchemes{default none;}\n"
          "interpolationSchemes{default linear;}\nsnGradSchemes{default corrected;}\n")
    write(os.path.join(case, "system/fvSolution"),
          HEADER.format(cls="dictionary", obj="fvSolution"))
    write(os.path.join(case, "constant/materialProperties"),
          HEADER.format(cls="dictionary", obj="materialProperties")
          + "type linearElastic;\nrho 2000;\nE 2e7;\nnu 0.3;\n" + material)
    run(case, "blockMesh")
    # cell centres (blockMesh order: x fastest, then y, then z)
    ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                             indexing="ij")
    idx = (ix.ravel(order="F"), iy.ravel(order="F"), iz.ravel(order="F"))
    xc = np.stack([(idx[0] + 0.5) * h, (idx[1] + 0.5) * h,
                   (idx[2] + 0.5) * h], axis=1)
    return dict(path=case, n=n, L=L, h=h, two_d=two_d, xc=xc, idx=idx,
                patches=patches)


def run(case, *cmd, np_=1, log=None, quiet_fail=False):
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1",
               OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    args = list(cmd) + ["-case", case]
    if np_ > 1:
        args = ["mpirun", "-np", str(np_), "--oversubscribe"] + args + ["-parallel"]
    p = subprocess.run(args, capture_output=True, text=True, env=env)
    if log:
        write(os.path.join(case, log), p.stdout + p.stderr)
    if p.returncode != 0:
        if not quiet_fail:
            print(p.stdout[-3000:], p.stderr[-3000:])
        raise RuntimeError(f"{' '.join(args)} failed")
    return p.stdout


CLASSES = {1: ("volScalarField", "scalar"), 3: ("volVectorField", "vector"),
           6: ("volSymmTensorField", "symmTensor")}
DIMS = {"rho": "[1 -3 0 0 0 0 0]", "U": "[0 1 -1 0 0 0 0]",
        "sigma": "[1 -1 -2 0 0 0 0]", "epsP": "[0 0 0 0 0 0 0]",
        "epsPEq": "[0 0 0 0 0 0 0]"}
NCMPT = {"rho": 1, "U": 3, "sigma": 6, "epsP": 6, "epsPEq": 1}


def write_fields(case, fields, time="0"):
    for name, vals in fields.items():
        vals = np.asarray(vals, dtype=float).reshape(len(case["xc"]), -1)
        cls, _ = CLASSES[vals.shape[1]]
        if vals.shape[1] == 1:
            body = "\n".join("%.17g" % x for x in vals[:, 0])
        else:
            body = "\n".join("(" + " ".join("%.17g" % x for x in row) + ")"
                             for row in vals)
        bf = ""
        for pname, ptype, _, _ in case["patches"]:
            t = ptype if ptype in ("cyclic", "empty") else "zeroGradient"
            bf += f"    {pname} {{ type {t}; }}\n"
        write(os.path.join(case["path"], time, name),
              HEADER.format(cls=cls, obj=name)
              + f"dimensions {DIMS[name]};\n"
              + f"internalField nonuniform List<{CLASSES[vals.shape[1]][1]}> "
              + f"{len(vals)}\n(\n{body}\n);\n"
              + "boundaryField\n{\n" + bf + "}\n")


def read_field(path, ncell, ncmpt):
    txt = open(path).read()
    m = re.search(r"internalField\s+uniform\s+(\([^)]*\)|\S+);", txt)
    if m:
        v = np.array([float(x) for x in m.group(1).strip("()").split()])
        return np.tile(v, (ncell, 1))
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", txt)
    start = m.end()
    end = txt.index("\n)", start)
    nums = np.array(txt[start:end].replace("(", " ").replace(")", " ").split(),
                    dtype=float)
    return nums.reshape(ncell, ncmpt)


def last_time(path):
    times = [d for d in os.listdir(path)
             if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"]
    return max(times, key=float)


def read_all(case, path=None):
    path = path or case["path"]
    t = last_time(path)
    return {f: read_field(os.path.join(path, t, f), len(case["xc"]), NCMPT[f])
            for f in NCMPT}


def remap_dict(case, **kw):
    body = "".join(f"{k} {v};\n" for k, v in kw.items())
    write(os.path.join(case["path"], "system/remapCheckDict"),
          HEADER.format(cls="dictionary", obj="remapCheckDict") + body)


def parse(out):
    res = {"ranges": {}}
    for line in out.splitlines():
        w = line.split()
        if not w:
            continue
        if w[0] == "steps":
            res["steps"], res["dt"] = int(w[1]), float(w[3])
        elif w[0] in ("maxCourantSum", "conservation", "cellSize"):
            res[w[0]] = float(w[1])
        elif w[0] == "range":
            res["ranges"][(w[1], int(w[2]))] = tuple(map(float, w[3:7]))
    return res


def vec(v):
    return "(" + " ".join(repr(float(x)) for x in v) + ")"


# ----------------------------------------------------------------------------
# Initial data
# ----------------------------------------------------------------------------

GAUSS3 = (np.array([-math.sqrt(0.6), 0.0, math.sqrt(0.6)]) / 2,
          np.array([5.0, 8.0, 5.0]) / 18)


def cell_average(case, f):
    """Cell averages of f(x, y, z) (vectorised) by 3-point Gauss per
    direction (z only in 3D)."""
    h = case["h"]
    xc = case["xc"]
    pts, wts = GAUSS3
    zs = [(0.0, 1.0)] if case["two_d"] else list(zip(pts, wts))
    acc = 0.0
    for px, wx in zip(pts, wts):
        for py, wy in zip(pts, wts):
            for pz, wz in zs:
                acc = acc + wx * wy * wz * f(xc[:, 0] + px * h,
                                             xc[:, 1] + py * h,
                                             xc[:, 2] + pz * h)
    return acc


def fill(s_rho, s_q, ncell):
    """Fields from shapes s_rho (density) and s_q (all other fields), each an
    array of cell values in [0, 1] (or [-1, 1]); constant offsets keep every
    component different."""
    s = s_q
    rho = 2000.0 * (1.0 + 0.25 * s_rho)
    U = np.stack([1.0 + 0.5 * s, -2.0 + s, 0.3 * s], axis=1)
    sig = np.stack([-5e4 + 1e4 * s, 1e3 + 2e3 * s, -3e3 * s,
                    -4e4 + 5e3 * s, 4e3 * s, -4e4 - 6e3 * s], axis=1)
    epsP = np.stack([0.01 * s, -0.004 * s, 0.002 * s, -0.003 * s,
                     0.001 * s, -0.003 * s], axis=1)
    epsPEq = 0.05 + 0.2 * s
    return {"rho": rho, "U": U, "sigma": sig, "epsP": epsP, "epsPEq": epsPEq}


def ones(case):
    return np.ones(len(case["xc"]))


# ----------------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------------

def check_conservation(tag, res):
    check(f"{tag}: mass, momentum, sum(m q) conserved to 1e-12 "
          f"(largest relative change {res['conservation']:.2e})",
          res["conservation"] < 1e-12)


def overshoot(res, field):
    """Largest new extremum over the run, relative to the initial range,
    over the field's components."""
    worst = 0.0
    for (f, c), (lo, hi, lo_r, hi_r) in res["ranges"].items():
        if f != field:
            continue
        span = hi - lo
        if span <= 0:
            continue
        worst = max(worst, (hi_r - hi) / span, (lo - lo_r) / span)
    return worst


def reference_1d(n, d, nsteps, rho, q_list, lims):
    """lr1d_plastic remap of a periodic row: rho by volume, each q mass-weighted
    (mass coordinates). d: face displacement per step."""
    dx = 1.0 / n
    h = np.full(n, dx)
    dd = np.full(n, d)
    qs = [q.copy() for q in q_list]
    for _ in range(nsteps):
        m = rho * h
        mass_flux = dd * lr1d.Bar1D.swept_average(rho, h, dd)
        m_new = m + mass_flux - np.roll(mass_flux, -1)
        for k, q in enumerate(qs):
            f = lr1d_plastic.mass_weighted_flux(q, m, mass_flux, lims[k])
            qs[k] = (m * q + f - np.roll(f, -1)) / m_new
        rho = m_new / dx
    return rho, qs


def case_1d(work):
    """A. Translation along x on a 3D mesh against the 1D prototype remap."""
    n = 64
    co = 0.4
    T = 0.3
    for variant in ("q", "rho"):
        case = make_case(work, "trans1d_" + variant, (n, 3, 3))
        x = case["xc"][:, 0]

        def gauss(x):
            return sum(np.exp(-((x - 0.3 + k) / 0.08) ** 2) for k in (-1, 0, 1))
        sq = ((case["idx"][0] >= 36) & (case["idx"][0] < 52)).astype(float)
        g = cell_average(case, lambda x, y, z: gauss(x))
        if variant == "q":
            f = fill(0 * g, g, len(x))
            # stress and history also get a square, to exercise the limiters
            f["sigma"][:, 0] += 2e4 * sq
            f["epsPEq"] = f["epsPEq"] + 0.3 * sq
            f["epsP"][:, 1] += 0.02 * sq
        else:
            f = fill(g + 2 * sq, 0 * g, len(x))
        write_fields(case, f)
        remap_dict(case, motion="translate", U="(1 0 0)", endTime=T, courant=co)
        res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
        out = read_all(case)
        nsteps, dt = res["steps"], res["dt"]

        # 1D reference along the row j = k = 0 (all rows are identical)
        row = (case["idx"][1] == 0) & (case["idx"][2] == 0)
        order = np.argsort(case["idx"][0][row])
        r1 = f["rho"][row][order]
        comps, lims, got, names = [], [], [], []
        for fname, lim in (("U", lr1d.limiter), ("sigma", lr1d.limiter),
                           ("epsP", lr1d_plastic.superbee),
                           ("epsPEq", lr1d_plastic.superbee)):
            a = np.asarray(f[fname]).reshape(len(x), -1)
            b = out[fname]
            for c in range(a.shape[1]):
                comps.append(a[row][order][:, c])
                got.append(b[row][order][:, c])
                lims.append(lim)
                names.append(f"{fname}[{c}]")
        rho_ref, q_ref = reference_1d(n, dt * 1.0, nsteps, r1, comps, lims)
        err = np.max(np.abs(out["rho"][row][order][:, 0] - rho_ref)) / np.max(np.abs(rho_ref))
        worst = ("rho", err)
        for nm, gq, rq in zip(names, got, q_ref):
            scale = max(np.max(np.abs(rq)), 1e-300)
            e = np.max(np.abs(gq - rq)) / scale
            if e > worst[1]:
                worst = (nm, e)
        # all rows equal
        spread = 0.0
        for fname in NCMPT:
            a = out[fname]
            ref = a[(case["idx"][1] == 0) & (case["idx"][2] == 0)][order]
            for j in range(3):
                for k in range(3):
                    rr = (case["idx"][1] == j) & (case["idx"][2] == k)
                    o = np.argsort(case["idx"][0][rr])
                    s = max(np.max(np.abs(ref)), 1e-300)
                    spread = max(spread, np.max(np.abs(a[rr][o] - ref)) / s)
        print(f"   A ({variant} profile): {nsteps} steps at Courant "
              f"{dt * n:.4f}; largest difference to lr1d/lr1d_plastic "
              f"{worst[1]:.2e} ({worst[0]}); rows differ by {spread:.2e}")
        check(f"A. x translation on 64x3x3 ({variant} profile) matches the 1D "
              f"prototype remap to 1e-10", worst[1] < 1e-10 and spread < 1e-10)
        check_conservation(f"A ({variant})", res)


def l1(a, b):
    return np.sum(np.abs(a - b)) / np.sum(np.abs(b))


def l1_dev(a, b):
    """L1 error relative to the L1 size of the deviation from the mean."""
    return np.sum(np.abs(a - b)) / np.sum(np.abs(b - np.mean(b)))


def smooth_order(work, dim, ns, co=DEFAULT_CO, tag="B"):
    """Smooth periodic profile translated along the diagonal for one box
    crossing; L1 errors of the MC fields against the initial data."""
    errs = {"rho": [], "U": [], "sigma": []}
    for n in ns:
        case = make_case(work, f"order{dim}d_{n}", (n, n, 1 if dim == 2 else n))

        def s(x, y, z):
            r = np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)
            if dim == 3:
                r = r * np.sin(2 * np.pi * z)
            return 0.5 + 0.5 * r

        sa = cell_average(case, s)

        def s2(x, y, z):
            r = np.cos(2 * np.pi * (x + 0.1)) * np.sin(2 * np.pi * (y - 0.2))
            if dim == 3:
                r = r * np.cos(2 * np.pi * z)
            return 0.5 + 0.5 * r
        # q is mass-weighted: initial cell values are int(rho q)/int(rho)
        rho_f = lambda x, y, z: 2000.0 * (1 + 0.25 * s2(x, y, z))  # noqa: E731
        rq = cell_average(case, lambda x, y, z: rho_f(x, y, z) * s(x, y, z))
        rr = cell_average(case, rho_f)
        f = fill(cell_average(case, s2), rq / rr, len(sa))
        write_fields(case, f)
        U = "(1 1 0)" if dim == 2 else "(1 1 1)"
        remap_dict(case, motion="translate", U=U, endTime=1, courant=co)
        res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
        out = read_all(case)
        errs["rho"].append(l1_dev(out["rho"][:, 0], f["rho"]))
        errs["U"].append(l1_dev(out["U"][:, 0], f["U"][:, 0]))
        errs["sigma"].append(l1_dev(out["sigma"][:, 0], f["sigma"][:, 0]))
        check_conservation(f"{tag} {dim}D n={n}", res)
    orders = {}
    for k, e in errs.items():
        p = [math.log(e[i] / e[i + 1]) / math.log(ns[i + 1] / ns[i])
             for i in range(len(ns) - 1)]
        fit = -np.polyfit(np.log(ns), np.log(e), 1)[0]
        orders[k] = (p, fit)
        print(f"   {tag} {dim}D {k}: L1 errors " + " ".join(f"{x:.3e}" for x in e)
              + "; orders " + " ".join(f"{x:.2f}" for x in p) + f"; fit {fit:.2f}")
    return orders


def case_order(work, quick):
    o2 = smooth_order(work, 2, (32, 64, 128))
    for k, (p, fit) in o2.items():
        check(f"B. 2D diagonal translation, {k}: order over 32/64/128 is "
              f"{fit:.2f} (>= 1.7)", fit >= 1.7)
    ns3 = (12, 24, 48)
    o3 = smooth_order(work, 3, ns3)
    for k, (p, fit) in o3.items():
        check(f"B. 3D body-diagonal translation, {k}: order over "
              f"{'/'.join(map(str, ns3))} is {fit:.2f} (>= 1.7)", fit >= 1.7)


def edge_width(profile, lo, hi):
    """Number of cells with value between 10 and 90 percent of the jump."""
    t = (profile - lo) / (hi - lo)
    return int(np.sum((t > 0.1) & (t < 0.9)))


def case_square(work, n=100, nprocs=(2, 4)):
    """C. Square profile, diagonal translation by 50 cells in x and y."""
    case = make_case(work, "square", (n, n, 1))
    i, j = case["idx"][0], case["idx"][1]
    sq = ((i >= 20) & (i < 60) & (j >= 20) & (j < 60)).astype(float)
    sq2 = ((i >= 30) & (i < 50) & (j >= 25) & (j < 55)).astype(float)
    f = fill(sq2, sq, len(sq))
    write_fields(case, f)
    shift = 50
    remap_dict(case, motion="translate", U="(1 1 0)", endTime=shift * case["h"],
               courant=DEFAULT_CO)
    res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
    out = read_all(case)
    for field in ("rho", "U", "sigma"):
        ov = overshoot(res, field)
        check(f"C. square, diagonal: no new extrema in {field} above 1e-10 "
              f"of the range (largest {ov:.1e})", ov <= 1e-10)
    ov = overshoot(res, "epsPEq")
    print(f"   C: history (superbee) new extrema {ov:.1e} of the range")
    # history edge width on the row through the middle of the moved square
    e = out["epsPEq"][:, 0]
    jm = (40 + shift) % n
    row = j == jm
    prof = e[row][np.argsort(i[row])]
    lo, hi = 0.05, 0.25
    # left edge near x = 20 + shift, right edge near x = 60 + shift
    w_left = edge_width(np.roll(prof, -(20 + shift - 10))[:20], lo, hi)
    w_right = edge_width(np.roll(prof, -(60 + shift - 10))[:20], lo, hi)
    colp = e[i == (40 + shift) % n][np.argsort(j[i == (40 + shift) % n])]
    w_bottom = edge_width(np.roll(colp, -(20 + shift - 10))[:20], lo, hi)
    w = max(w_left, w_right, w_bottom)
    print(f"   C: history edge widths (10-90 %) after {shift} cells in x and y "
          f"({shift * math.sqrt(2):.0f} cells of travel): left {w_left}, "
          f"right {w_right}, bottom {w_bottom}")
    check(f"C. history edge width at most 4 cells (largest {w})", w <= 4)
    check_conservation("C", res)

    # F. Parallel: same case on 2 and 4 processors
    for npr in nprocs:
        pc = case["path"] + f"_np{npr}"
        if os.path.isdir(pc):
            shutil.rmtree(pc)
        shutil.copytree(case["path"], pc,
                        ignore=shutil.ignore_patterns("[1-9]*", "0.*"))
        coeffs = "(2 1 1)" if npr == 2 else "(2 2 1)"
        write(os.path.join(pc, "system/decomposeParDict"),
              HEADER.format(cls="dictionary", obj="decomposeParDict")
              + f"numberOfSubdomains {npr};\nmethod simple;\n"
              f"coeffs {{ n {coeffs}; }}\n")
        run(pc, "decomposePar", log="log.decomposePar")
        resp = parse(run(pc, "celRemapCheck", np_=npr, log="log.remapCheck"))
        run(pc, "reconstructPar", "-latestTime", log="log.reconstructPar")
        outp = read_all(case, pc)
        worst = 0.0
        for fname in NCMPT:
            a, b = outp[fname], out[fname]
            scale = np.maximum(np.max(np.abs(b), axis=0), 1e-300)
            worst = max(worst, float(np.max(np.abs(a - b) / scale)))
        check(f"F. square diagonal on {npr} processors matches serial to 1e-12 "
              f"(largest difference {worst:.1e})", worst <= 1e-12)
        check_conservation(f"F ({npr} processors)", resp)


def case_disc(work, n=100):
    """D. Slotted disc of history, solid-body rotation, one turn."""
    case = make_case(work, "disc", (n, n, 1), periodic=False)
    h = case["h"]

    def disc(x, y):
        inside = (x - 0.5) ** 2 + (y - 0.75) ** 2 <= 0.15 ** 2
        slot = (np.abs(x - 0.5) <= 0.025) & (y <= 0.85)
        return (inside & ~slot).astype(float)
    xc = case["xc"]
    sub = (np.arange(4) + 0.5) / 4 - 0.5
    d = sum(disc(xc[:, 0] + a * h, xc[:, 1] + b * h) for a in sub for b in sub) / 16
    f = fill(0 * d, 0 * d, len(d))
    f["epsPEq"] = d
    f["epsP"][:, 0] = 0.01 * d
    write_fields(case, f)
    omega = 2 * np.pi
    remap_dict(case, motion="rotate", omega=omega, centre="(0.5 0.5 0)",
               rIn=0.42, rOut=0.49, endTime=1.0, courant=DEFAULT_CO)
    res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
    out = read_all(case)
    e = out["epsPEq"][:, 0]
    lo, hi = e.min(), e.max()
    ov = max(hi - 1.0, -lo)
    shape = np.sum(np.abs(e - d)) / np.sum(d)
    print(f"   D: {res['steps']} steps; history in [{lo:.4f}, {hi:.4f}]; "
          f"L1 shape error {shape:.3f} of the disc; total "
          f"{np.sum(e):.15g} vs {np.sum(d):.15g}")
    check_conservation("D. slotted disc", res)
    check(f"D. slotted disc: no new extrema above 1 percent (largest "
          f"{100 * ov:.2f} %)", ov <= 0.01)
    # cells beyond the moving region are untouched, bit for bit
    r = np.hypot(xc[:, 0] - 0.5, xc[:, 1] - 0.5)
    far = r > 0.49 + 2 * h
    same = all(np.array_equal(out[k].reshape(len(d), -1)[far],
                              np.asarray(f[k]).reshape(len(d), -1)[far])
               for k in NCMPT)
    check(f"D. {int(far.sum())} cells with no flux keep every value bit for bit",
          same)
    return shape


def case_uniform(work):
    """E. A uniform state stays uniform."""
    n = 40
    # also reads a non-default remapLimiters subdictionary
    case = make_case(work, "uniform_swirl", (n, n, 1),
                     material="remapLimiters { rho minmod; U upwind; epsP MC; }\n")
    f = fill(ones(case), ones(case), n * n)
    write_fields(case, f)
    remap_dict(case, motion="swirl", A=0.1, L=1, endTime=1.0, courant=0.4)
    res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
    out = read_all(case)
    worst = max(float(np.max(np.abs(out[k] - np.asarray(f[k]).reshape(n * n, -1)))
                      / max(np.max(np.abs(f[k])), 1e-300)) for k in NCMPT)
    print(f"   E: swirl {res['steps']} steps (max Courant sum "
          f"{res['maxCourantSum']:.2f}): largest change {worst:.1e}")
    check(f"E. uniform state under a non-uniform divergence-free flux stays "
          f"uniform to 1e-12 (largest change {worst:.1e})", worst < 1e-12)
    check_conservation("E (swirl)", res)

    case = make_case(work, "uniform_compress", (n, n, 1))
    f = fill(ones(case), ones(case), n * n)
    write_fields(case, f)
    remap_dict(case, motion="compress", A=1.0, L=1, endTime=0.01, nSteps=2,
               reverse="true")
    res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
    out = read_all(case)
    worst_q = max(float(np.max(np.abs(out[k] - np.asarray(f[k]).reshape(n * n, -1)))
                        / max(np.max(np.abs(f[k])), 1e-300))
                  for k in NCMPT if k != "rho")
    drho = float(np.max(np.abs(out["rho"][:, 0] - f["rho"])) / f["rho"][0])
    print(f"   E: compress then reverse (Courant sum {res['maxCourantSum']:.3f}): "
          f"velocity, stress, history change {worst_q:.1e}; density {drho:.1e}")
    check(f"E. uniform velocity, stress and history stay uniform to 1e-12 under a "
          f"compressing flux and its reverse (largest change {worst_q:.1e})",
          worst_q < 1e-12)
    check_conservation("E (compress)", res)
    return drho


def case_walls(work):
    """G. Compression along x in a box with walls: every row along x stays the
    same (round-off fluxes on faces parallel to the motion must not matter)."""
    case = make_case(work, "walls_rows", (40, 10, 1), periodic=False)
    i = case["idx"][0]
    sq = ((i >= 10) & (i < 20)).astype(float)
    write_fields(case, fill(sq, sq, len(sq)))
    remap_dict(case, motion="compress", A=0.5, L=1, endTime=0.1, nSteps=16)
    res = parse(run(case["path"], "celRemapCheck", log="log.remapCheck"))
    out = read_all(case)
    worst = 0.0
    for k, v in out.items():
        a = v.reshape(10, 40, -1)
        worst = max(worst, float(np.max(np.abs(a - a[0])))
                    / max(float(np.max(np.abs(a))), 1e-300))
    check(f"G. compression along x between walls: rows agree to 1e-12 "
          f"(largest difference {worst:.1e})", worst <= 1e-12)
    check_conservation("G", res)


def study(work):
    """Largest stable material Courant number (report only)."""
    for dim, n in ((2, 64), (3, 24)):
        for co in (0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6):
            case = make_case(work, f"study{dim}_{co}", (n, n, 1 if dim == 2 else n))
            i, j, k = case["idx"]
            sq = ((i >= n // 4) & (i < 3 * n // 4) & (j >= n // 4) & (j < 3 * n // 4)
                  & ((k >= n // 4) & (k < 3 * n // 4) if dim == 3 else True)).astype(float)
            write_fields(case, fill(sq, sq, len(sq)))
            U = "(1 1 0)" if dim == 2 else "(1 1 1)"
            remap_dict(case, motion="translate", U=U, endTime=1, courant=co)
            try:
                res = parse(run(case["path"], "celRemapCheck", log="log", quiet_fail=True))
                ov = max(overshoot(res, fld) for fld in ("rho", "U", "sigma"))
                print(f"   study {dim}D Courant {co} per direction (sum "
                      f"{res['maxCourantSum']:.2f}): new extrema {ov:.1e}")
            except RuntimeError:
                print(f"   study {dim}D Courant {co}: failed (non-positive density)")


def main():
    args = sys.argv[1:]
    work = tempfile.mkdtemp(prefix="remap_")
    keep = "--keep" in args
    try:
        case_1d(work)
        case_uniform(work)
        case_walls(work)
        case_square(work)
        case_disc(work)
        case_order(work, "--quick" in args)
        if "--study" in args:
            study(work)
    finally:
        if not keep:
            shutil.rmtree(work, ignore_errors=True)
        else:
            print("kept", work)
    if FAILS:
        print(f"{len(FAILS)} check(s) failed")
        sys.exit(1)
    print("all remap checks passed")


if __name__ == "__main__":
    main()
