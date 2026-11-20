"""WP N.2: a granular layer sheared by a moving lid, under gravity, on three
meshes (decision 0111). Simple case for the mesh-dependence law (Sia,
2026-11-20: not the silo).

Layer H = 40 d (d = 1 mm), one column periodic in x (cyclic), base fixed,
lid moving at V along x (both fixedValue U, so no slip), gravity -y, started
at rest under gravity with P0 under the lid (K0 = nu/(1-nu)); a uniform
start lost P0 (the top went into tension and was cut). The friction mu = tau/P falls with
depth (tau is uniform in steady shear, P grows), so the flow gathers under
the lid: with the local law the band is as thin as the cells; with the
nonlocal law its width is set by xi ~ A d / sqrt|mu - mu_s|.

    python3 shear_layer.py WORKDIR [local|nonlocal] [n ...]
"""

import math
import os
import re
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
RHO, E, NU, G = 1500.0, 1.0e6, 0.3, 9.81
H = 40 * LAW["d"]
P0, V = 200.0, 0.02
T_END = float(os.environ.get("CEL_T_END", "1.0"))


def field_file(cls, obj, dims, internal, bc):
    return (cs.HEADER.format(cls=cls, obj=obj) + f"dimensions {dims};\n"
            f"internalField uniform {internal};\nboundaryField\n{{\n{bc}}}\n")


def make(p, n, nonlocal_on):
    w = H / n
    os.makedirs(f"{p}/system", exist_ok=True)
    os.makedirs(f"{p}/constant", exist_ok=True)
    os.makedirs(f"{p}/0", exist_ok=True)
    v = [(0, 0, 0), (w, 0, 0), (w, H, 0), (0, H, 0)]
    v += [(x, y, w) for x, y, _ in v]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    cs.write(f"{p}/system/blockMeshDict", cs.HEADER.format(cls="dictionary", obj="blockMeshDict")
             + f"scale 1;\nvertices\n(\n{vs}\n);\nblocks\n(\n"
             f"    hex (0 1 2 3 4 5 6 7) (1 {n} 1) simpleGrading (1 1 1)\n);\n"
             "boundary\n(\n"
             "    left { type cyclic; neighbourPatch right; faces ((0 4 7 3)); }\n"
             "    right { type cyclic; neighbourPatch left; faces ((1 2 6 5)); }\n"
             "    bottom { type patch; faces ((0 1 5 4)); }\n"
             "    top { type patch; faces ((3 7 6 2)); }\n"
             "    frontAndBack { type empty; faces ((0 3 2 1) (4 5 6 7)); }\n);\n")
    cs.write(f"{p}/system/controlDict", cs.HEADER.format(cls="dictionary", obj="controlDict")
             + f"application celFoam;\nstartFrom startTime; startTime 0; stopAt endTime;\n"
             f"endTime {T_END!r}; deltaT 1e-6;\nwriteControl adjustableRunTime;\n"
             f"writeInterval {T_END / 2!r}; writeFormat ascii; writePrecision 12;\n"
             "maxCo 0.4;\n")
    cs.write(f"{p}/system/fvSchemes", cs.HEADER.format(cls="dictionary", obj="fvSchemes")
             + "ddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
             "divSchemes { default none; }\nlaplacianSchemes { default none; }\n"
             "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    cs.write(f"{p}/system/fvSolution", cs.HEADER.format(cls="dictionary", obj="fvSolution"))
    phi = math.degrees(math.atan(LAW["muS"]))
    mat = (f"rho {RHO};\nfreeSurface no;\ntype mohrCoulomb;\nE {E!r};\nnu {NU!r};\n"
           f"phi {phi!r};\nc 0;\npsi 0;\nmu2 {LAW['mu2']!r};\nI0 {LAW['I0']!r};\n"
           f"d {LAW['d']!r};\nrhoS {LAW['rhoS']!r};\n")
    if nonlocal_on:
        mat += f"nonlocalA {LAW['A']!r};\n"
    mat += os.environ.get("CEL_EXTRA_MATERIAL", "")
    mat += os.environ.get("CEL_EXTRA_MATERIAL", "")
    cs.write(f"{p}/constant/materialProperties",
             cs.HEADER.format(cls="dictionary", obj="materialProperties") + mat)
    cs.write(f"{p}/constant/g", cs.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G!r} 0);\n")
    cyc = "    left { type cyclic; }\n    right { type cyclic; }\n" \
          "    frontAndBack { type empty; }\n"
    cs.write(f"{p}/0/U", field_file(
        "volVectorField", "U", "[0 1 -1 0 0 0 0]", "(0 0 0)",
        "    bottom { type fixedValue; value uniform (0 0 0); }\n"
        f"    top {{ type fixedValue; value uniform ({V!r} 0 0); }}\n" + cyc))
    # at rest under gravity: s_yy = -(P0 + rho g (H - y)), s_xx = s_zz = K0 s_yy
    y = (np.arange(n) + 0.5) * w
    syy = -(P0 + RHO * G * (H - y))
    k0 = NU / (1 - NU)
    rows = "\n".join(f"({float(k0 * a)!r} 0 0 {float(a)!r} 0 {float(k0 * a)!r})" for a in syy)
    cs.write(f"{p}/0/sigma", field_file(
        "volSymmTensorField", "sigma", "[1 -1 -2 0 0 0 0]", "(0 0 0 0 0 0)",
        "    bottom { type zeroGradient; }\n    top { type zeroGradient; }\n" + cyc)
        .replace("internalField uniform (0 0 0 0 0 0);",
                 f"internalField nonuniform List<symmTensor> {n}\n(\n{rows}\n);"))
    return p


def run(p):
    for cmd in (["blockMesh"], ["celFoam"]):
        with open(f"{p}/log.{cmd[0]}", "w") as f:
            r = subprocess.run(cmd, cwd=p, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode:
            raise RuntimeError(open(f"{p}/log.{cmd[0]}").read()[-3000:])


def read_internal(fn, ncomp):
    """Internal field of an ascii vector or symmTensor file, one row per cell."""
    s = open(fn).read()
    body = s[s.index("internalField"):]
    body = body[:body.index(";")]
    if "nonuniform" not in body:
        raise ValueError(f"{fn}: uniform field")
    body = body[body.index("(", body.index(">")) + 1:]
    rows = re.findall(r"\(([^()]*)\)", body)
    return np.array([[float(x) for x in r.split()] for r in rows])


def profile(p):
    t = f"{T_END:g}"
    U = read_internal(f"{p}/{t}/U", 3)
    S = read_internal(f"{p}/{t}/sigma", 6)
    n = len(U)
    y = (np.arange(n) + 0.5) * H / n
    return y, U[:, 0], S[:, 1], -S[:, 3]        # y, u, tau (xy), P (-yy)


def reference(y, tau, P):
    """1D steady nonlocal solution for the run's own tau and P(y): g from
    k g - g'' = k g_loc, gamma_dot = g mu, u(0) = 0."""
    law = ngf.NGF(**LAW)
    nf = 4000
    yf = (np.arange(nf) + 0.5) * H / nf
    Pf = np.interp(yf, y, P)
    muf = abs(tau) / Pf
    gf = law.solve_layer(yf, muf, Pf) if hasattr(law, "solve_layer") else \
        ngf.solve_layer(law, yf, muf, Pf)
    rate = gf * muf
    uf = np.concatenate([[0], np.cumsum(rate) * H / nf])[:-1] + 0.5 * rate * H / nf
    return yf, uf


def band(y, u):
    """Height over which u goes from 0.1 to 0.9 of the lid speed [m]."""
    yy = np.concatenate([[0], y, [H]])
    uu = np.concatenate([[0], u, [V]])
    return float(np.interp(0.9 * V, uu, yy) - np.interp(0.1 * V, uu, yy))


def main():
    work = os.path.abspath(sys.argv[1])
    mode = sys.argv[2] if len(sys.argv) > 2 else "nonlocal"
    ns = [int(x) for x in sys.argv[3:]] or [10, 20, 40]
    for n in ns:
        p = make(f"{work}/{mode}_{n}", n, mode == "nonlocal")
        run(p)
        y, u, tau, P = profile(p)
        yf, uf = reference(y, tau.mean(), P)
        err = np.abs(u - np.interp(y, yf, uf)).max() / V
        print(f"{mode} {n} cells ({H / n / LAW['d']:g} d): band {band(y, u) / LAW['d']:.1f} d, "
              f"tau {tau.mean():.1f} Pa (spread {np.ptp(tau):.2g}), P {P.min():.0f}-{P.max():.0f} Pa, "
              f"mu {abs(tau.mean()) / P.max():.3f}-{abs(tau.mean()) / P.min():.3f}; "
              f"profile vs 1D nonlocal {err * 100:.1f} % of V", flush=True)


if __name__ == "__main__":
    main()
