"""WP N.1 done-when checks: nonlocal granular fluidity (decision 0111).

[1] Kernel (tests/nonlocal/kernelCheck.C, g++): homogeneous simple shear at
    fixed pressure, Mohr-Coulomb and Drucker-Prager with nonlocalA, g = g_loc,
    dt 1e-6 and 1e-4 s, shear rates 1, 10, 50 1/s: mu at the end equals the
    local law mu_s + b I to a relative 1e-6.
[2] Fluidity field (celFluidityCheck): a layer of 40 cells (H = 40 d) with
    P and mu linear in y, mu from below mu_s to above it, so the flowing
    part spreads fluidity into the part below mu_s; zero gradient at both
    ends. g within 1 % (of its largest value) of an independent fine-grid
    solution of k g - g'' = k g_loc (prototype/constitutive/ngf.py,
    4000 cells). Again with the top quarter void: the material part against
    the same solution over the material only.

Glass-bead constants of Henann and Kamrin (2013): mu_s 0.3819, mu_2 0.6435,
I0 0.278, A 0.48; d 1 mm, rho_s 2450.

    python3 nonlocal_checks.py WORKDIR
"""

import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "prototype", "constitutive"))
sys.path.insert(0, os.path.join(ROOT, "tests", "lagrange"))
import ngf  # noqa: E402
import cases as cs  # noqa: E402

LAW = dict(muS=0.3819, mu2=0.6435, I0=0.278, A=0.48, d=1.0e-3, rhoS=2450.0)
fails = []


def report(name, ok, detail):
    print(f"  {name}: {'PASS' if ok else 'FAIL'} ({detail})", flush=True)
    if not ok:
        fails.append(name)


def check_kernel(work):
    print("[1] Kernel: homogeneous simple shear", flush=True)
    k = os.path.join(ROOT, "src", "celConstitutive", "kernel")
    exe = os.path.join(work, "kernelCheck")
    srcs = [os.path.join(k, f) for f in sorted(os.listdir(k)) if f.endswith(".C")]
    subprocess.run(["g++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror", f"-I{k}",
                    os.path.join(HERE, "kernelCheck.C"), *srcs, "-o", exe], check=True)
    out = subprocess.run([exe], check=True, capture_output=True, text=True).stdout.split("\n")
    worst = 0.0
    for line in filter(None, out):
        name, dt, rate, mu, ref = line.split()
        worst = max(worst, abs(float(mu) / float(ref) - 1))
    report("mu = mu_s + b I", worst <= 1e-6 and len(out) > 10,
           f"{len(list(filter(None, out)))} runs, largest relative error {worst:.1e}")


def layer_case(p, n, H, P0, P1, mu0, mu1, empty=None):
    os.makedirs(f"{p}/system", exist_ok=True)
    os.makedirs(f"{p}/constant", exist_ok=True)
    w = H / n
    v = [(0, 0, 0), (w, 0, 0), (w, H, 0), (0, H, 0)]
    v += [(x, y, w) for x, y, _ in v]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    cs.write(f"{p}/system/blockMeshDict", cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n"
             f"    hex (0 1 2 3 4 5 6 7) (1 {n} 1) simpleGrading (1 1 1)\n);\n"
             "boundary\n(\n"
             "    bottom { type wall; faces ((0 1 5 4)); }\n"
             "    top { type patch; faces ((3 7 6 2)); }\n"
             "    sides { type empty; faces ((0 4 7 3) (1 2 6 5) (0 3 2 1) (4 5 6 7)); }\n);\n")
    cs.write(f"{p}/system/controlDict", cs.HEADER.format(cls="dictionary", obj="controlDict")
             + "application celFluidityCheck;\nstartFrom startTime; startTime 0; stopAt endTime;"
             " endTime 1; deltaT 1;\nwriteControl timeStep; writeInterval 1;\n")
    cs.write(f"{p}/system/fvSchemes", cs.HEADER.format(cls="dictionary", obj="fvSchemes")
             + "ddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
             "divSchemes { default none; }\nlaplacianSchemes { default none; }\n"
             "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    cs.write(f"{p}/system/fvSolution", cs.HEADER.format(cls="dictionary", obj="fvSolution"))
    law = dict(phi=math.degrees(math.atan(LAW["muS"])), mu2=LAW["mu2"], I0=LAW["I0"],
               nonlocalA=LAW["A"], d=LAW["d"], rhoS=LAW["rhoS"],
               P0=P0, P1=P1, mu0=mu0, mu1=mu1)
    if empty is not None:
        law["emptyAbove"] = empty
    cs.write(f"{p}/system/fluidityCheckDict",
             cs.HEADER.format(cls="dictionary", obj="fluidityCheckDict")
             + "".join(f"{k} {x!r};\n" for k, x in law.items()))
    for cmd in (["blockMesh"], ["celFluidityCheck"]):
        with open(f"{p}/log.{cmd[0]}", "w") as f:
            subprocess.run(cmd, cwd=p, stdout=f, stderr=subprocess.STDOUT, check=True)
    return np.loadtxt(f"{p}/fluidity.dat")


def check_field(work):
    print("[2] Fluidity across a layer with mu through mu_s", flush=True)
    law = ngf.NGF(**LAW)
    H, n = 40 * LAW["d"], 40
    P0, P1 = 2000.0, -20000.0            # 2000 Pa at the bottom, 1200 at the top
    mu0, mu1 = 0.30, 0.20 / H            # 0.30 at the bottom, 0.50 at the top
    for name, top in (("full layer", H), ("top quarter void", 0.75 * H)):
        r = layer_case(f"{work}/{name.replace(' ', '_')}", n, H, P0, P1, mu0, mu1,
                       None if top == H else top)
        y, g = r[:, 0], r[:, 1]
        nf = 4000
        yf = (np.arange(nf) + 0.5) * top / nf
        gf = ngf.solve_layer(law, yf, mu0 + mu1 * yf, P0 + P1 * yf)
        m = y < top
        ref = np.interp(y[m], yf, gf)
        err = np.abs(g[m] - ref).max() / ref.max()
        below = (mu0 + mu1 * y[m]) < LAW["muS"]
        creep = g[m][below].max() / ref.max()
        void = np.abs(g[~m]).max() if (~m).any() else 0.0
        report(name, err <= 0.01 and creep > 0.01 and void <= 1e-12 * ref.max(),
               f"largest error {err * 100:.2f} % of g_max {ref.max():.1f} 1/s; "
               f"fluidity below mu_s up to {creep * 100:.0f} % of g_max (nonlocal); "
               f"void {void:.1e}")


def main():
    work = os.path.abspath(sys.argv[1])
    os.makedirs(work, exist_ok=True)
    check_kernel(work)
    check_field(work)
    print(f"{len(fails)} failures")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
