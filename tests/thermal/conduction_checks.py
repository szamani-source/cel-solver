"""WP T.1 done-when: conduction in the grains of a fixed bed (decision 0136).

The solid is an elastic bed at rest (no gravity, E 0.1 MPa so the solid step is
long), rho 1500 kg/m^3, with thermalProperties model equilibrium, c 1000
J/(kg K) and k 1500 W/(m K): diffusivity a = 1e-3 m^2/s, so heat moves
centimetres in a second (the closed forms do not care about the values).

Cases:
A  Bar 0.2 m, 100 cells of 2 mm, the left face held 10 K above the start
   (fixedValue), other faces adiabatic; t = 1 s. Against the erfc step of
   tools/thermal_reference.py: largest error under 1 % of the step; heat in
   within 1 % of 2 k dT sqrt(t/(pi a)) per area. Energy: thermal energy -
   start - heat in under 1e-12 of the thermal energy (round-off).
A2 The same bar heated by a constant flux q = 5000 W/m^2 (fixedGradient
   q/k). Against the constant-flux closed form: under 1 % of the face rise.
B  Box 0.1 m on cells leaning 27 deg (20 x 20, x shifted by y tan 27 deg),
   T = 300 + 50 x + 30 y everywhere at the start and on every face
   (fixedValue): a linear temperature is steady, and the scheme with the
   non-orthogonal correction keeps it to 1e-10 of its range over 200 steps.
C  The same leaning box, all faces adiabatic, a hot square in the middle, run
   0.1 s: energy constant to 1e-12; on 2 (2 x 1) and 4 (2 x 2) processors T
   matches serial to 1e-10 of the hot-cold difference.
D  Temperature off: case A run without thermalProperties and with model none
   gives rho, U and sigma byte for byte equal to the run with temperature;
   neither prints anything about temperature.

    python3 conduction_checks.py WORKDIR
"""

import math
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import thermal_reference as tr  # noqa: E402

RHO, E, NU = 1500.0, 1.0e5, 0.3
C_HEAT, K = 1000.0, 1500.0
A = K / (RHO * C_HEAT)
T0, DT_STEP, Q = 293.15, 10.0, 5000.0
SHEAR = math.tan(math.radians(27.0))

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""

failed = []


def report(name, good, detail):
    print(("ok  " if good else "BAD ") + name + ": " + detail, flush=True)
    if not good:
        failed.append(name)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def sh(path, cmd, log, env=None):
    env = dict(env or os.environ, PWD=os.path.abspath(path))
    with open(f"{path}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n"
                           + open(f"{path}/log.{log}").read()[-3000:])


def run(path, nProcs=1):
    if nProcs > 1:
        sh(path, ["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        sh(path, ["mpirun", "-np", str(nProcs), "--oversubscribe", "celFoam", "-parallel"],
           "celFoam", env)
        sh(path, ["reconstructPar"], "reconstructPar")
    else:
        sh(path, ["celFoam"], "celFoam")


def common(path, *, end, write_interval, decomp=(1, 1, 1), thermal="equilibrium"):
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
    write(f"{path}/constant/g", HEADER.format(cls="uniformDimensionedVectorField", obj="g")
          + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    write(f"{path}/constant/materialProperties",
          HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO!r};
type linearElastic;
E {E!r};
nu {NU!r};
walls {{ left {{ friction 0; }} right {{ friction 0; }} bottom {{ friction 0; }} top {{ friction 0; }} }}
""")
    if thermal is not None:
        write(f"{path}/constant/thermalProperties",
              HEADER.format(cls="dictionary", obj="thermalProperties") + f"""
model        {thermal};
c            {C_HEAT!r};
conductivity constant;
k            {K!r};
""")


def blockmesh(path, W, H, nx, ny, shear=0.0):
    dz = W / nx
    v = [(0, 0), (W, 0), (W + H * shear, H), (H * shear, H)]
    pts = [(x, y, 0.0) for x, y in v] + [(x, y, dz) for x, y in v]
    write(f"{path}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices ({' '.join('(%r %r %r)' % p for p in pts)});
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
    sh(path, ["blockMesh"], "blockMesh")


def write_T(path, internal, patches):
    """internal: list or scalar; patches: {name: (type, extra)}."""
    if np.isscalar(internal):
        field = f"uniform {internal!r}"
    else:
        field = "nonuniform List<scalar> %d (%s)" % (
            len(internal), " ".join(repr(float(x)) for x in internal))
    bf = []
    for name, (typ, extra) in patches.items():
        bf.append(f"    {name} {{ type {typ}; {extra} }}")
    bf.append("    frontAndBack { type empty; }")
    write(f"{path}/0/T", HEADER.format(cls="volScalarField", obj="T") + f"""
dimensions [0 0 0 1 0 0 0];
internalField {field};
boundaryField
{{
{chr(10).join(bf)}
}}
""")


def read_internal(path):
    s = open(path).read()
    s = s[s.index("internalField"):]
    m = re.match(r"internalField\s+uniform\s+(\S+);", s)
    if m:
        return None
    n = int(re.search(r"List<scalar>\s*(\d+)", s).group(1))
    body = s[s.index("(") + 1:]
    return np.array(body.split()[:n], dtype=float)


def read_centres(path):
    """Cell centres and the boundary face centres of each patch from 0/C."""
    s = open(f"{path}/0/C").read()
    def vec_list(txt):
        n = int(re.search(r"List<vector>\s*(\d+)", txt).group(1))
        rows = re.findall(r"\(([^()]*)\)", txt[txt.index("(", txt.index("List")):])[:n]
        return np.array([list(map(float, r.split())) for r in rows])
    inner = s[s.index("internalField"):s.index("boundaryField")]
    cells = vec_list(inner)
    patches = {}
    bf = s[s.index("boundaryField"):]
    for name in ("left", "right", "bottom", "top"):
        i = bf.index(name)
        j = bf.index("}", bf.index("value", i))
        patches[name] = vec_list(bf[i:j])
    return cells, patches


def log_value(path, key):
    vals = re.findall(rf"^{re.escape(key)}\s*=\s*(\S+)", open(f"{path}/log.celFoam").read(), re.M)
    return [float(v) for v in vals]


def steps(path):
    return len(re.findall(r"^Time = ", open(f"{path}/log.celFoam").read(), re.M))


# ---------------------------------------------------------------- case A, A2

def bar(path, kind, thermal="equilibrium"):
    common(path, end=1.0, write_interval=0.5, thermal=thermal)
    blockmesh(path, 0.2, 0.002, 100, 1)
    left = ("fixedValue", f"value uniform {T0 + DT_STEP!r};") if kind == "step" \
        else ("fixedGradient", f"gradient uniform {Q / K!r};")
    write_T(path, T0, {"left": left, "right": ("zeroGradient", ""),
                       "bottom": ("zeroGradient", ""), "top": ("zeroGradient", "")})
    t0 = time.perf_counter()
    run(path)
    return time.perf_counter() - t0


def check_bar(work):
    print("[A] Bar with the left face held hot (erfc)", flush=True)
    p = os.path.join(work, "barStep")
    wall = bar(p, "step")
    T = read_internal(f"{p}/1/T")
    x = (np.arange(100) + 0.5) * 0.002
    exact = np.array([tr.erfc_step(xi, 1.0, A, T0, DT_STEP) for xi in x])
    err = np.abs(T - exact).max() / DT_STEP
    report("erfc step at t = 1 s", err < 0.01,
           f"largest error {err * 100:.3f} % of the step (limit 1 %), {steps(p)} steps, {wall:.0f} s")
    E, Qin, bal = log_value(p, "Thermal energy"), log_value(p, "Heat in"), log_value(p, "Energy balance")
    rel = max(abs(b) for b in bal) / E[0]
    report("energy balance with a fixed-temperature face", rel < 1e-12,
           f"energy - start - heat in = {bal[-1]:.3e} J (heat in {Qin[-1]:.6e} J), "
           f"{rel:.1e} of the energy {E[0]:.6e} J")
    # heat in against the closed form: k dT sqrt(t/(pi a)) * 2 * area
    area = 0.002 * 0.002
    q_exact = 2 * K * DT_STEP * math.sqrt(1.0 / (math.pi * A)) * area
    report("heat in against 2 k dT sqrt(t/(pi a)) per area", abs(Qin[-1] / q_exact - 1) < 0.01,
           f"{Qin[-1]:.6e} J against {q_exact:.6e} J ({(Qin[-1] / q_exact - 1) * 100:+.3f} %)")

    print("[A2] Bar heated by a constant flux", flush=True)
    p2 = os.path.join(work, "barFlux")
    bar(p2, "flux")
    T = read_internal(f"{p2}/1/T")
    exact = np.array([tr.flux_step(xi, 1.0, A, K, Q, T0) for xi in x])
    rise = tr.flux_step(0.0, 1.0, A, K, Q) - 0.0
    err = np.abs(T - exact).max() / rise
    report("constant-flux face at t = 1 s", err < 0.01,
           f"largest error {err * 100:.3f} % of the face rise {rise:.4f} K (limit 1 %)")
    Qin = log_value(p2, "Heat in")[-1]
    want = Q * 0.002 * 0.002 * 1.0
    report("heat in equals q A t", abs(Qin / want - 1) < 1e-10,
           f"{Qin:.12e} J against {want:.12e} J")
    return p


# ------------------------------------------------------------------- case B, C

def leaning_box(path, end, write_interval, decomp=(1, 1, 1)):
    common(path, end=end, write_interval=write_interval, decomp=decomp)
    blockmesh(path, 0.1, 0.1, 20, 20, SHEAR)
    write_T(path, T0, {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})
    sh(path, ["postProcess", "-func", "writeCellCentres", "-time", "0"], "centres")
    return read_centres(path)


def check_linear(work):
    print("[B] Linear temperature on cells leaning 27 deg", flush=True)
    p = os.path.join(work, "linear")
    # the solid step on 5 mm cells is about 2.1e-4 s: 200 steps
    cells, faces = leaning_box(p, 0.042, 0.042)
    lin = lambda c: 300.0 + 50.0 * c[:, 0] + 30.0 * c[:, 1]
    patches = {n: ("fixedValue", "value nonuniform List<scalar> %d (%s);"
                   % (len(c), " ".join(repr(float(v)) for v in lin(c))))
               for n, c in faces.items()}
    write_T(p, lin(cells), patches)
    run(p)
    T = read_internal(f"{p}/0.042/T")
    rng = lin(cells).max() - lin(cells).min()
    err = np.abs(T - lin(cells)).max() / rng
    report("linear temperature stays steady on leaning cells", err < 1e-10,
           f"largest change {err:.2e} of the range over {steps(p)} steps (limit 1e-10)")


def hot_square(path, decomp=(1, 1, 1)):
    cells, _ = leaning_box(path, 0.1, 0.05, decomp)
    x = cells[:, 0] - cells[:, 1] * SHEAR
    y = cells[:, 1]
    hot = (np.abs(x - 0.05) < 0.02) & (np.abs(y - 0.05) < 0.02)
    write_T(path, np.where(hot, T0 + DT_STEP, T0),
            {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})


def check_adiabatic(work):
    print("[C] Adiabatic leaning box, serial and 2/4 processors", flush=True)
    res = {}
    for n, dec in ((1, (1, 1, 1)), (2, (2, 1, 1)), (4, (2, 2, 1))):
        p = os.path.join(work, f"adiabatic{n}")
        hot_square(p, dec)
        run(p, n)
        res[n] = p
    p = res[1]
    E = log_value(p, "Thermal energy")
    bal = log_value(p, "Energy balance")
    rel = max(abs(b) for b in bal) / E[0]
    report("energy constant with adiabatic walls", rel < 1e-12,
           f"largest |energy - start| {max(abs(b) for b in bal):.3e} J of {E[0]:.6e} J ({rel:.1e})")
    Ts = read_internal(f"{p}/0.1/T")
    spread = Ts.max() - Ts.min()
    report("heat has spread", spread < 0.9 * DT_STEP and spread > 0.01,
           f"hot-cold difference {spread:.4f} K at 0.1 s (was {DT_STEP} K)")
    for n in (2, 4):
        Tp = read_internal(f"{res[n]}/0.1/T")
        d = np.abs(Tp - Ts).max() / DT_STEP
        report(f"{n} processors match serial", d < 1e-10,
               f"largest difference {d:.2e} of the step (limit 1e-10)")
        bp = log_value(res[n], "Energy balance")
        report(f"{n} processors conserve energy", max(abs(b) for b in bp) / E[0] < 1e-12,
               f"largest |balance| {max(abs(b) for b in bp):.3e} J")


# ---------------------------------------------------------------------- case D

def check_off(work, pA):
    print("[D] Temperature off: solid unchanged", flush=True)
    same = True
    detail = []
    for label, thermal in (("no file", None), ("model none", "none")):
        p = os.path.join(work, "bar_" + label.replace(" ", ""))
        bar(p, "step", thermal=thermal)
        log = open(f"{p}/log.celFoam").read()
        for t in ("0.5", "1"):
            for f in ("rho", "U", "sigma"):
                a = os.path.join(pA, t, f)
                b = os.path.join(p, t, f)
                if os.path.exists(a) != os.path.exists(b) or \
                        (os.path.exists(a) and open(a, "rb").read() != open(b, "rb").read()):
                    same = False
                    detail.append(f"{label} {t}/{f} differs")
        quiet = "Thermal energy" not in log and "Temperature (WP" not in log
        if not quiet:
            same = False
            detail.append(f"{label} printed temperature output")
        if os.path.exists(f"{p}/1/T"):
            same = False
            detail.append(f"{label} wrote 1/T")
    report("temperature off leaves the solid byte for byte", same,
           "; ".join(detail) or "rho, U, sigma identical at 0.5 and 1 s, no temperature output")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    pA = check_bar(work)
    check_linear(work)
    check_adiabatic(work)
    check_off(work, pA)
    print(f"conduction checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_conduction: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
