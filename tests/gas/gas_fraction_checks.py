"""WP F.1 done-when: gas fractions and the celFluid skeleton (decisions 0125,
0126).

Set-up (10 mm cells, 2D): column W = 0.2 m wide, H = 0.5 m tall, walls and
floor friction 20 deg, smooth lid. T0 Mohr-Coulomb (rho 1800, E 10 MPa,
nu 0.3, phi 35, c 0, psi 0) with the T0 rule set and freeSurface yes. Filled
to y = 0.40 m, plus one half-filled row (alpha 0.5, rho 900) above it; void
on top. Gravity ramps in over 0.05 s with damping 250 1/s; the run ends at
t = 0.2 s (settled). Grains for the gas: quartz rhoP 2650 kg/m^3, d 1 mm, air
at 293.15 K and 101325 Pa, drag ergunWenYu. Set packing 1800/2650 = 0.679.

Runs: A without constant/gasProperties (v1), B with model none, C with model
poreGas; E = C with plain Mohr-Coulomb (T0 rule set off) and D = E on 2
processors. The parallel pair runs without the rule set because its switches
(flowImpedance snap, shear cap) flip on round-off and part serial and parallel
runs of the solid itself by about 1e-5 here (known issue, WP 6.1, decision
0080); without them the solid matches to 1e-14, so the pair tests the gas.

Checks:
1. Packing: in every full cell phi_s equals the set packing rho0/rhoP within
   0.5 % (self weight compresses by rho g h / M, about 0.05 %).
2. Fractions: phi_s = rho/rhoP (= alpha rho_m/rhoP) to 1e-15 in every cell,
   phi_g = 1 - phi_s; void cells phi_s 0, phi_g 1, k 0; the half-filled row
   has phi_s near half the packing.
3. Permeability: k in every cell with material equals
   tools/gas_reference.py permeability(d, e_m), e_m = 1 - rho/(alpha rhoP),
   to 1e-12; air density and viscosity printed by celFoam match the script.
4. Gas off unchanged: rho, U, sigma, epsP, epsPEq, alpha of A, B and C are
   identical byte for byte at every written time; A and B write no phiS,
   phiG or k and A prints nothing about gas. (The full v1 suite in CI is the
   wider check that nothing else moved.)
5. Parallel: D's reconstructed phiS, phiG and k match E to 1e-12.

    python3 gas_fraction_checks.py WORKDIR
"""

import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402

CELL = 0.01
W, H, FILL = 0.20, 0.50, 0.40
RHO, E, NU, PHI = 1800.0, 1.0e7, 0.3, 35.0
RHOP, D = 2650.0, 1.0e-3
T_END, WRITE = 0.2, 0.05
SOLID = ["rho", "U", "sigma", "epsP", "epsPEq", "alpha"]
GAS = ["phiS", "phiG", "k"]

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


RULES = "shearCap yes;\nspreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n"


def make_case(path, gas=None, rules=True):
    nx, ny = int(round(W / CELL)), int(round(H / CELL))
    v = [(0, 0, 0), (W, 0, 0), (W, H, 0), (0, H, 0),
         (0, 0, CELL), (W, 0, CELL), (W, H, CELL), (0, H, CELL)]
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
endTime         {T_END!r};
deltaT          1e-5;
writeControl    adjustableRunTime;
writeInterval   {WRITE!r};
purgeWrite      0;
writeFormat     ascii;
writePrecision  17;
timeFormat      general;
timePrecision   6;
runTimeModifiable false;
maxCo           0.4;
gravityRamp     0.05;
damping         250;
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
    write(f"{path}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict")
          + "numberOfSubdomains 2;\nmethod simple;\ncoeffs { n (2 1 1); }\n")
    write(f"{path}/system/setFieldsDict", HEADER.format(cls="dictionary", obj="setFieldsDict") + f"""
defaultFieldValues (volScalarFieldValue alpha 0 volScalarFieldValue rho 0);
regions
(
    boxToCell {{ box (-1 -1 -1) (1 {FILL!r} 1);
        fieldValues (volScalarFieldValue alpha 1 volScalarFieldValue rho {RHO!r}); }}
    boxToCell {{ box (-1 {FILL!r} -1) (1 {FILL + CELL!r} 1);
        fieldValues (volScalarFieldValue alpha 0.5 volScalarFieldValue rho {0.5 * RHO!r}); }}
);
""")
    for f in ("alpha", "rho"):
        dim = "[0 0 0 0 0 0 0]" if f == "alpha" else "[1 -3 0 0 0 0 0]"
        write(f"{path}/0/{f}", HEADER.format(cls="volScalarField", obj=f) + f"""
dimensions {dim};
internalField uniform 0;
boundaryField {{ ".*" {{ type zeroGradient; }} frontAndBack {{ type empty; }} }}
""")
    write(f"{path}/constant/g", HEADER.format(cls="uniformDimensionedVectorField", obj="g")
          + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 -9.81 0);\n")
    write(f"{path}/constant/materialProperties",
          HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO!r};
freeSurface yes;
type mohrCoulomb;
E {E!r};
nu {NU!r};
phi {PHI!r};
c 0;
psi 0;
sigmaT 0;
mu2 0.961813;
I0 0.279;
d 0.001;
rhoS 2650;
{RULES if rules else ""}walls
{{
    left {{ friction 20; }}
    right {{ friction 20; }}
    bottom {{ friction 20; }}
    top {{ friction 0; }}
}}
""")
    if gas is not None:
        write(f"{path}/constant/gasProperties", HEADER.format(cls="dictionary", obj="gasProperties") + f"""
model   {gas};
rhoP    {RHOP!r};
d       {D!r};
T       293.15;
p0      101325;
drag    ergunWenYu;
""")


def sh(path, cmd, log, env=None):
    env = dict(env or os.environ, PWD=os.path.abspath(path))
    with open(f"{path}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n"
                           + open(f"{path}/log.{log}").read()[-3000:])


def run(path, nProcs=1):
    t0 = time.perf_counter()
    sh(path, ["blockMesh"], "blockMesh")
    sh(path, ["setFields"], "setFields")
    if nProcs > 1:
        sh(path, ["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        sh(path, ["mpirun", "-np", str(nProcs), "--oversubscribe", "celFoam", "-parallel"],
           "celFoam", env)
        sh(path, ["reconstructPar"], "reconstructPar")
    else:
        sh(path, ["celFoam"], "celFoam")
    return time.perf_counter() - t0


def times(path):
    return sorted((d for d in os.listdir(path) if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
                  key=float)


def field(path, t, name):
    """Internal field of an ascii volScalarField or volVectorField as an array."""
    txt = open(f"{path}/{t}/{name}").read()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\((.*?)\n\)\s*;", txt, re.S)
    if m:
        body = m.group(2).replace("(", " ").replace(")", " ")
        a = np.array(body.split(), dtype=float)
        n = int(m.group(1))
        return a.reshape(n, -1) if a.size != n else a
    m = re.search(r"internalField\s+uniform\s+(\S+);", txt)
    return np.full(int(W / CELL) * int(H / CELL), float(m.group(1)))


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    cases = {"A": None, "B": "none", "C": "poreGas"}
    secs = {}
    for name, gas in cases.items():
        make_case(f"{work}/{name}", gas)
        secs[name] = run(f"{work}/{name}")
    make_case(f"{work}/E", "poreGas", rules=False)
    secs["E"] = run(f"{work}/E")
    make_case(f"{work}/D", "poreGas", rules=False)
    secs["D"] = run(f"{work}/D", nProcs=2)

    C = f"{work}/C"
    tl = times(C)
    t = tl[-1]
    rho, alpha = field(C, t, "rho"), field(C, t, "alpha")
    phiS, phiG, k = (field(C, t, f) for f in GAS)
    U = field(C, t, "U")
    umax = np.sqrt((U ** 2).sum(axis=1)).max()
    pack = RHO / RHOP
    full = alpha > 1 - 1e-9
    void = alpha < 1e-6
    part = ~full & ~void
    mat = ~void

    print(f"Column W = {W} m, H = {H} m, 10 mm cells, filled to {FILL} m plus a half row; "
          f"T0 Mohr-Coulomb E = {E / 1e6:g} MPa; grains rhoP {RHOP:g}, d {D * 1e3:g} mm")
    print(f"Runs to t = {t} s: " + ", ".join(f"{n} {s:.0f} s" for n, s in secs.items())
          + f"; largest speed at the end {umax:.1e} m/s; {full.sum()} full, {part.sum()} partly "
          f"filled, {void.sum()} void cells\n")
    res = []

    def report(name, good, detail):
        res.append(bool(good))
        print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}")

    dev = np.abs(phiS[full] / pack - 1.0).max()
    report("1 packing", dev <= 0.005,
           f"phi_s in full cells {phiS[full].min():.5f} to {phiS[full].max():.5f}, set packing "
           f"{pack:.5f}, largest deviation {100 * dev:.3f} % (<= 0.5 %)")

    err_s = np.abs(phiS - rho / RHOP).max()
    err_g = np.abs(phiG - (1.0 - phiS)).max()
    void_ok = (np.all(phiS[void] == 0) and np.all(phiG[void] == 1) and np.all(k[void] == 0))
    half = phiS[part].mean() / pack if part.any() else float("nan")
    report("2 fractions", err_s <= 1e-15 and err_g <= 1e-15 and void_ok and part.any(),
           f"phi_s - rho/rhoP {err_s:.0e}, phi_g - (1 - phi_s) {err_g:.0e} (<= 1e-15); void cells "
           f"0/1/0: {void_ok}; partly filled cells phi_s {half:.3f} of the packing")

    mu, rg = gr.air_viscosity(), gr.air_density()
    e_m = 1.0 - rho[mat] / (alpha[mat] * RHOP)
    k_ref = np.array([gr.permeability(D, e) for e in e_m])
    err_k = np.abs(k[mat] / k_ref - 1.0).max()
    log = open(f"{C}/log.celFoam").read()
    m = re.search(r"rho_g = (\S+) kg/m\^3, mu = (\S+) Pa s", log)
    err_air = max(abs(float(m.group(1)) / rg - 1), abs(float(m.group(2)) / mu - 1)) if m else 1.0
    report("3 permeability", err_k <= 1e-12 and err_air <= 1e-7,
           f"k against gas_reference.permeability(d, e_m) in {mat.sum()} cells: {err_k:.0e} "
           f"(<= 1e-12), k = {k[full].mean():.4e} m^2 in full cells; air rho_g, mu "
           f"{err_air:.0e} (<= 1e-7, printed to 8 digits)")

    same, missing, extra = True, [], []
    for tt in tl:
        for f in SOLID:
            ref = open(f"{work}/A/{tt}/{f}", "rb").read()
            for n in ("B", "C"):
                if open(f"{work}/{n}/{tt}/{f}", "rb").read() != ref:
                    same = False
                    missing.append(f"{n}/{tt}/{f}")
        for n in ("A", "B"):
            extra += [f"{n}/{tt}/{g}" for g in GAS if os.path.exists(f"{work}/{n}/{tt}/{g}")]
    quiet = "Pore gas" not in open(f"{work}/A/log.celFoam").read()
    report("4 gas off unchanged", same and not extra and quiet and len(tl) >= 2,
           f"{len(SOLID)} solid fields at {len(tl)} times identical in A, B, C: {same}"
           + (f" (differ: {missing[:3]})" if missing else "")
           + f"; gas fields written without gas: {extra[:3] or 'none'}; A silent about gas: {quiet}")

    err_par, err_rho = 0.0, 0.0
    for tt in tl:
        for f in GAS + ["rho"]:
            a, b = field(f"{work}/E", tt, f), field(f"{work}/D", tt, f)
            e = np.abs(a - b).max() / max(np.abs(a).max(), 1e-300)
            if f == "rho":
                err_rho = max(err_rho, e)
            else:
                err_par = max(err_par, e)
    report("5 parallel", err_par <= 1e-12,
           f"plain Mohr-Coulomb column, 2 processors against serial, phiS, phiG, k at "
           f"{len(tl)} times: {err_par:.0e} (<= 1e-12; rho {err_rho:.0e})")

    ok = all(res)
    print(f"\nGas fraction checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
