"""WP 2.3: build and run box cases for celLagrangeStepCheck.

A case is a box of n cells per side (side L), 2D (one cell thick, thickness
L/n, empty front and back) or 3D. Each side is one of: cyclic (periodic
pair), empty, fixed (U fixedValue), free (U zeroGradient: traction-free) or
symmetry. Used by run_checks.py.
"""

import os
import subprocess

import numpy as np

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""

# hex (0 1 2 3 4 5 6 7) faces
FACES = {
    "xMin": "(0 4 7 3)", "xMax": "(1 2 6 5)",
    "yMin": "(0 1 5 4)", "yMax": "(3 7 6 2)",
    "zMin": "(0 3 2 1)", "zMax": "(4 5 6 7)",
}

# Material of decision 0013 (WP 1.3b)
RHO = 2000.0
E, NU = 20.0e6, 0.3
MATERIALS = {
    "linearElastic": "type linearElastic;\nE {E};\nnu {NU};\n",
    "vonMises": "type vonMises;\nE {E};\nnu {NU};\nsigma_y 10e3;\nH 0;\n",
    "druckerPrager": ("type druckerPrager;\nE {E};\nnu {NU};\nphi 30;\nc 5e3;\n"
                      "psi 0;\nfit compression;\n"),
}


def moduli():
    K = E / (3 * (1 - 2 * NU))
    G = E / (2 * (1 + NU))
    return K, G, K + 4 * G / 3


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def foam_value(v):
    if isinstance(v, (list, tuple, np.ndarray)):
        return "(" + " ".join(foam_value(x) for x in v) + ")"
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, np.integer):
        return str(int(v))
    return str(v)


def make_case(path, n, sides, check, dim=2, L=1.0, material="linearElastic",
              maxCo=None, solveMomentum=True, secondOrder=True, nProcs=1,
              split=None, ny=None):
    """sides: dict side name -> kind; missing sides default to cyclic
    (x, y and, in 3D, z) or empty (z in 2D). ny: cells in y (default n),
    same cell size."""
    h = L / n
    nz, Lz = (n, L) if dim == 3 else (1, h)
    ny = n if ny is None else ny
    Ly = ny * h
    kinds = {}
    for s in FACES:
        default = "empty" if (dim == 2 and s.startswith("z")) else "cyclic"
        kinds[s] = sides.get(s, default)

    patches = []
    for s, k in kinds.items():
        if k == "cyclic":
            other = s[0] + ("Max" if s.endswith("Min") else "Min")
            body = f"type cyclic;\n        neighbourPatch {other};"
        elif k == "empty":
            body = "type empty;"
        elif k == "symmetry":
            body = "type symmetry;"
        else:
            body = "type patch;"
        patches.append(f"    {s}\n    {{\n        {body}\n        faces ({FACES[s]});\n    }}")

    write(f"{path}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices
(
    (0 0 0) ({L} 0 0) ({L} {Ly} 0) (0 {Ly} 0)
    (0 0 {Lz}) ({L} 0 {Lz}) ({L} {Ly} {Lz}) (0 {Ly} {Lz})
);
blocks ( hex (0 1 2 3 4 5 6 7) ({n} {ny} {nz}) simpleGrading (1 1 1) );
boundary
(
{chr(10).join(patches)}
);
""")
    extra = ""
    if maxCo is not None:
        extra += f"maxCo {maxCo!r};\n"
    extra += f"solveMomentum {'yes' if solveMomentum else 'no'};\n"
    extra += f"secondOrder {'yes' if secondOrder else 'no'};\n"
    write(f"{path}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application celLagrangeStepCheck;
startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 17;
{extra}""")
    write(f"{path}/system/fvSchemes", HEADER.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default none; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write(f"{path}/system/fvSolution", HEADER.format(cls="dictionary", obj="fvSolution"))
    if split is None:
        split = (nProcs, 1, 1)
    write(f"{path}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict") + f"""
numberOfSubdomains {nProcs};
method simple;
coeffs {{ n ({split[0]} {split[1]} {split[2]}); }}
""")
    write(f"{path}/constant/materialProperties",
          HEADER.format(cls="dictionary", obj="materialProperties")
          + f"rho {RHO};\n" + MATERIALS[material].format(E=E, NU=NU))

    bcs = []
    for s, k in kinds.items():
        t = {"cyclic": "cyclic", "empty": "empty", "symmetry": "symmetry",
             "fixed": "fixedValue;\n        value uniform (0 0 0)",
             "free": "zeroGradient"}[k]
        bcs.append(f"    {s}\n    {{\n        type {t};\n    }}")
    write(f"{path}/0/U", HEADER.format(cls="volVectorField", obj="U") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{{
{chr(10).join(bcs)}
}}
""")
    body = "\n".join(f"{k} {foam_value(v)};" for k, v in check.items())
    write(f"{path}/system/lagrangeCheckDict",
          HEADER.format(cls="dictionary", obj="lagrangeCheckDict") + body + "\n")
    return path


def run(path, nProcs=1):
    """blockMesh, then the driver (decomposed on nProcs if > 1)."""
    def sh(cmd, log):
        with open(f"{path}/log.{log}", "w") as f:
            r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            tail = open(f"{path}/log.{log}").read()[-3000:]
            raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n{tail}")

    sh(["blockMesh"], "blockMesh")
    if nProcs > 1:
        sh(["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1",
                   OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        with open(f"{path}/log.driver", "w") as f:
            r = subprocess.run(["mpirun", "-np", str(nProcs), "--oversubscribe",
                                "celLagrangeStepCheck", "-parallel"], cwd=path,
                               stdout=f, stderr=subprocess.STDOUT, env=env)
        if r.returncode != 0:
            raise RuntimeError(open(f"{path}/log.driver").read()[-3000:])
    else:
        sh(["celLagrangeStepCheck"], "driver")
    return path


def history(path):
    return np.loadtxt(f"{path}/history.dat")


def cells(path, nProcs=1):
    if nProcs == 1:
        return np.loadtxt(f"{path}/cells.dat")
    return np.vstack([np.loadtxt(f"{path}/cells.dat.proc{p}") for p in range(nProcs)])


def log_value(path, key):
    for line in open(f"{path}/log.driver"):
        if line.startswith(key):
            return float(line.split("=")[1].split()[0])
    return None
