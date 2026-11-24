"""Cases for the pore-pressure tests (WP F.2, decision 0127).

Shared by tests/gas/pore_pressure_checks.py (done-when test) and the gas
validation pack test 1, tests/validation/gas/ergunBed. The solid is frozen:
no gravity, the bed at rest at a uniform packing, so only the gas moves.
All runs use the ergunWenYu drag law (Sia, 2026-11-24).
"""

import os
import re
import subprocess
import time

import numpy as np

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""

RHOP = 2650.0       # quartz grains [kg/m^3]
P0 = 101325.0
T = 293.15


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def make_case(path, *, size, cells, e, d, end, write_interval, p_patches,
              inflow=None, p_set=None, skew=0.0, decomp=(1, 1, 1), E=1.0e6,
              extra_gas="", coupled=False):
    """2D box size = (Lx, Ly), cells = (nx, ny), packing gas fraction e,
    grain diameter d. p_patches maps left/right/bottom/top to a 0/p entry
    body. inflow: dict patch -> superficial velocity at p0. p_set: setFields
    box for p, ((x0, y0), (x1, y1), value). skew shifts the top edge in x
    (non-orthogonal cells). coupled: gas acts on the solid (WP F.3); the
    F.2 cases keep the solid frozen (coupled no)."""
    Lx, Ly = size
    nx, ny = cells
    rho = (1.0 - e) * RHOP
    v = [(0, 0, 0), (Lx, 0, 0), (Lx + skew, Ly, 0), (skew, Ly, 0)]
    dz = min(Lx / nx, Ly / ny)
    v = v + [(x, y, dz) for x, y, _ in v]
    write(f"{path}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices ({' '.join('(%r %r %r)' % p for p in v)});
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1));
boundary
(
    left {{ type patch; faces ((0 4 7 3)); }}
    right {{ type patch; faces ((1 2 6 5)); }}
    bottom {{ type patch; faces ((0 1 5 4)); }}
    top {{ type patch; faces ((3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
""")
    write(f"{path}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
application     celFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end!r};
deltaT          1e-5;
writeControl    adjustableRunTime;
writeInterval   {write_interval!r};
purgeWrite      0;
writeFormat     ascii;
writePrecision  17;
timeFormat      general;
timePrecision   8;
runTimeModifiable false;
maxCo           0.4;
""")
    write(f"{path}/system/fvSchemes", HEADER.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default none; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    write(f"{path}/system/fvSolution", HEADER.format(cls="dictionary", obj="fvSolution"))
    n = decomp[0] * decomp[1] * decomp[2]
    write(f"{path}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict")
          + f"numberOfSubdomains {n};\nmethod simple;\ncoeffs {{ n ({decomp[0]} {decomp[1]} {decomp[2]}); }}\n")
    for f, dim, val in (("alpha", "[0 0 0 0 0 0 0]", 1.0), ("rho", "[1 -3 0 0 0 0 0]", rho)):
        write(f"{path}/0/{f}", HEADER.format(cls="volScalarField", obj=f) + f"""
dimensions {dim};
internalField uniform {val!r};
boundaryField {{ ".*" {{ type zeroGradient; }} frontAndBack {{ type empty; }} }}
""")
    bf = "\n".join(f"    {name} {{ {body} }}" for name, body in p_patches.items())
    write(f"{path}/0/p", HEADER.format(cls="volScalarField", obj="p") + f"""
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{{
    ".*" {{ type zeroGradient; }}
{bf}
    frontAndBack {{ type empty; }}
}}
""")
    if p_set:
        (x0, y0), (x1, y1), val = p_set
        write(f"{path}/system/setFieldsDict", HEADER.format(cls="dictionary", obj="setFieldsDict") + f"""
defaultFieldValues ();
regions ( boxToCell {{ box ({x0!r} {y0!r} -1) ({x1!r} {y1!r} 1);
    fieldValues (volScalarFieldValue p {val!r}); }} );
""")
    write(f"{path}/constant/g", HEADER.format(cls="uniformDimensionedVectorField", obj="g")
          + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    write(f"{path}/constant/materialProperties",
          HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {rho!r};
type linearElastic;
E {E!r};
nu 0.3;
""")
    inl = ""
    if inflow:
        inl = "inflow\n{\n" + "".join(f"    {k} {v!r};\n" for k, v in inflow.items()) + "}\n"
    write(f"{path}/constant/gasProperties", HEADER.format(cls="dictionary", obj="gasProperties") + f"""
model   poreGas;
rhoP    {RHOP!r};
d       {d!r};
T       {T!r};
p0      {P0!r};
drag    ergunWenYu;
coupled {'yes' if coupled else 'no'};
{inl}{extra_gas}""")


def sh(path, cmd, log, env=None):
    env = dict(env or os.environ, PWD=os.path.abspath(path))
    with open(f"{path}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n"
                           + open(f"{path}/log.{log}").read()[-3000:])


def run(path, nProcs=1, set_fields=False):
    """Mesh, set up and run; returns the celFoam wall time [s]."""
    sh(path, ["blockMesh"], "blockMesh")
    if set_fields:
        sh(path, ["setFields"], "setFields")
    t0 = time.perf_counter()
    if nProcs > 1:
        sh(path, ["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        t0 = time.perf_counter()
        sh(path, ["mpirun", "-np", str(nProcs), "--oversubscribe", "celFoam", "-parallel"],
           "celFoam", env)
        secs = time.perf_counter() - t0
        sh(path, ["reconstructPar"], "reconstructPar")
        return secs
    sh(path, ["celFoam"], "celFoam")
    return time.perf_counter() - t0


def times(path):
    return sorted((d for d in os.listdir(path) if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
                  key=float)


def field(path, t, name, n=None):
    """Internal field of an ascii volScalarField or volVectorField."""
    txt = open(f"{path}/{t}/{name}").read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\((.*?)\n\)\s*;", txt, re.S)
    if m:
        body = m.group(2).replace("(", " ").replace(")", " ")
        a = np.array(body.split(), dtype=float)
        k = int(m.group(1))
        return a.reshape(k, -1) if a.size != k else a
    m = re.search(r"internalField\s+uniform\s+(\S+);", txt)
    return np.full(n, float(m.group(1)))


def patch_value(path, t, name, patch):
    """Face values of a patch of a volScalarField (value entry)."""
    txt = open(f"{path}/{t}/{name}").read()
    blk = re.search(patch + r"\s*\{(.*?)\n    \}", txt, re.S).group(1)
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;", blk, re.S)
    if m:
        return np.array(m.group(2).split(), dtype=float)
    return np.array([float(re.search(r"value\s+uniform\s+(\S+);", blk).group(1))])


def patch_gradient(path, t, name, patch):
    """Gradient entry of a fixedGradient patch (uniform or first value)."""
    txt = open(f"{path}/{t}/{name}").read()
    blk = re.search(patch + r"\s*\{(.*?)\n    \}", txt, re.S).group(1)
    m = re.search(r"gradient\s+uniform\s+(\S+);", blk)
    if m:
        return float(m.group(1))
    m = re.search(r"gradient\s+nonuniform\s+List<scalar>\s*\d+\s*\((.*?)\)", blk, re.S)
    return np.array(m.group(1).split(), dtype=float)


def mass_log(path):
    """Gas mass and net inflow printed each step: arrays (mass, netIn)."""
    rows = re.findall(r"Pore gas: (\d+) solves \(\d+ Newton, \d+ iterations\), change (\S+) Pa, mass (\S+) kg, net in (\S+) kg",
                      open(f"{path}/log.celFoam").read())
    a = np.array([[float(r[2]), float(r[3]), float(r[0]), float(r[1])] for r in rows])
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def cell_centres(path):
    sh(path, ["postProcess", "-func", "writeCellCentres", "-time", "0"], "centres")
    return field(path, "0", "C")
