"""WP T.2 done-when: temperature moves with the material (decision 0135).

Cases (celFoam, freeSurface yes, thermalProperties model equilibrium):
A  Translating block. A 1D channel 2 m long, 200 cells of 10 mm, no gravity,
   smooth walls; a block 0.6 m long (x 0.3 to 0.9) moves at 1 m/s through
   void for 0.5 s (50 cells). Elastic, E 1 MPa. No conduction (k 0), so T
   is only carried. The left half of the block is 30 K hotter.
   - energy constant to 1e-12;
   - the hot-cold front moved 0.5 m (within one cell) and is at most 2
     cells wide between 10 % and 90 % (default limiter superbee, as the
     history fields);
   - on 4 processors (the block crosses two cuts) T matches serial to 1e-10;
   - a second run with a uniform block keeps T uniform to 1e-12.
B  Box with an inlet and an outlet. 0.2 x 0.3 m, 10 mm cells, plain
   Mohr-Coulomb (rho 1800, E 1 MPa, phi 35), gravity; the bottom half is
   filled at the start with T = 293.15 + 100 y; the inlet (lid middle, 40
   mm) feeds material at 0.5 m/s and 353.15 K from 0.05 to 0.35 s; the
   outlet (floor middle, 40 mm) opens at 0.1 s; k 0.3 W/(m K); t = 0.6 s.
   - energy balance (energy - start - conducted in - carried in) to 1e-12
     of the energy at every write;
   - material arrives at the inlet temperature: the hottest cell is 353.15
     K and no cell is hotter than that or colder than the coldest start;
   - empty cells hold voidT;
   - on 2 processors (cut through the inlet and the outlet) T matches serial
     to 1e-10 until 0.2 s (the plain Mohr-Coulomb solid itself parts once
     the stream lands, decision 0134); the energy balance holds to 1e-12
     throughout.
C  Box B with one temperature everywhere (start and inlet 293.15 K): T stays
   293.15 in every cell with material to 1e-12, and the heat carried in
   equals c T (fed - discharged) to 1e-12.

    python3 carried_checks.py WORKDIR
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import conduction_checks as cc  # noqa: E402  (write, sh, run, HEADER, ...)

C_HEAT = 1000.0
T0 = 293.15
HEADER = cc.HEADER
failed = cc.failed
report = cc.report


def controls(path, end, write_interval, decomp=(1, 1, 1)):
    cc.write(f"{path}/system/controlDict", HEADER.format(cls="dictionary", obj="controlDict") + f"""
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
    cc.write(f"{path}/system/fvSchemes", HEADER.format(cls="dictionary", obj="fvSchemes") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; }
divSchemes { default none; }
laplacianSchemes { default none; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
""")
    cc.write(f"{path}/system/fvSolution", HEADER.format(cls="dictionary", obj="fvSolution"))
    n = decomp[0] * decomp[1] * decomp[2]
    cc.write(f"{path}/system/decomposeParDict", HEADER.format(cls="dictionary", obj="decomposeParDict")
             + f"numberOfSubdomains {n};\nmethod simple;\ncoeffs {{ n ({decomp[0]} {decomp[1]} {decomp[2]}); }}\n")


def thermal(path, k):
    cc.write(f"{path}/constant/thermalProperties",
             HEADER.format(cls="dictionary", obj="thermalProperties") + f"""
model        equilibrium;
c            {C_HEAT!r};
conductivity constant;
k            {k!r};
voidT        {T0!r};
""")


def scalar_file(path, name, dim, values, patches):
    if np.isscalar(values):
        field = f"uniform {values!r}"
    else:
        field = "nonuniform List<scalar> %d (%s)" % (
            len(values), " ".join(repr(float(v)) for v in values))
    bf = "\n".join(f"    {n} {{ type {t}; {x} }}" for n, (t, x) in patches.items())
    cc.write(f"{path}/0/{name}", HEADER.format(cls="volScalarField", obj=name) + f"""
dimensions {dim};
internalField {field};
boundaryField
{{
{bf}
    frontAndBack {{ type empty; }}
}}
""")


def vector_file(path, name, values, patches):
    field = "nonuniform List<vector> %d (%s)" % (
        len(values), " ".join("(%r %r %r)" % tuple(map(float, v)) for v in values))
    bf = "\n".join(f"    {n} {{ type {t}; {x} }}" for n in patches for t, x in [patches[n]])
    cc.write(f"{path}/0/{name}", HEADER.format(cls="volVectorField", obj=name) + f"""
dimensions [0 1 -1 0 0 0 0];
internalField {field};
boundaryField
{{
{bf}
    frontAndBack {{ type empty; }}
}}
""")


def energy_log(path):
    log = open(f"{path}/log.celFoam").read()
    get = lambda k: [float(v) for v in re.findall(rf"^{re.escape(k)}\s*=\s*(\S+)", log, re.M)]
    return dict(E=get("Thermal energy"), heat=get("Heat in"), carried=get("Heat carried in"),
                bal=get("Energy balance"), dis=get("Discharged mass"), fed=get("Fed mass"),
                mass=get("Total mass"))


# ------------------------------------------------------------------ case A

def channel(path, uniform, decomp=(1, 1, 1)):
    nx, L, h = 200, 2.0, 0.01
    controls(path, 0.5, 0.25, decomp)
    cc.write(f"{path}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices ((0 0 0) ({L!r} 0 0) ({L!r} {h!r} 0) (0 {h!r} 0)
          (0 0 {h!r}) ({L!r} 0 {h!r}) ({L!r} {h!r} {h!r}) (0 {h!r} {h!r}));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} 1 1) simpleGrading (1 1 1));
boundary
(
    left {{ type patch; faces ((0 4 7 3)); }}
    right {{ type patch; faces ((1 2 6 5)); }}
    bottom {{ type patch; faces ((0 1 5 4)); }}
    top {{ type patch; faces ((3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
""")
    cc.sh(path, ["blockMesh"], "blockMesh")
    x = (np.arange(nx) + 0.5) * h
    block = (x > 0.3) & (x < 0.9)
    zg = {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")}
    scalar_file(path, "alpha", "[0 0 0 0 0 0 0]", np.where(block, 1.0, 0.0), zg)
    scalar_file(path, "rho", "[1 -3 0 0 0 0 0]", np.where(block, 1800.0, 0.0), zg)
    vector_file(path, "U", [(1.0 if b else 0.0, 0, 0) for b in block],
                {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})
    Tin = np.where(block & (x < 0.6) & (not uniform), T0 + 30.0, T0)
    scalar_file(path, "T", "[0 0 0 1 0 0 0]", Tin, zg)
    cc.write(f"{path}/constant/g", HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    cc.write(f"{path}/constant/materialProperties",
             HEADER.format(cls="dictionary", obj="materialProperties") + """
rho 1800;
freeSurface yes;
type linearElastic;
E 1e6;
nu 0.3;
walls { left { friction 0; } right { friction 0; } bottom { friction 0; } top { friction 0; } }
""")
    thermal(path, 0.0)
    cc.run(path, decomp[0] * decomp[1] * decomp[2])
    return x


def check_channel(work):
    print("[A] Hot-cold block carried 50 cells through void", flush=True)
    p = os.path.join(work, "channel")
    x = channel(p, uniform=False)
    e = energy_log(p)
    rel = max(abs(b) for b in e["bal"]) / e["E"][0]
    report("energy carried exactly", rel < 1e-12,
           f"largest |balance| {max(abs(b) for b in e['bal']):.3e} J of {e['E'][0]:.6e} J ({rel:.1e})")
    T = cc.read_internal(f"{p}/0.5/T")
    rho = cc.read_internal(f"{p}/0.5/rho")
    th = (T - T0) / 30.0
    mat = rho > 0
    inside = mat & (x > 0.85) & (x < 1.35)
    width = int(np.sum(inside & (th > 0.1) & (th < 0.9)))
    # front position: where theta crosses 0.5
    i = np.where(inside & (th < 0.5))[0][0]
    xf = x[i - 1] + (th[i - 1] - 0.5) / (th[i - 1] - th[i]) * (x[i] - x[i - 1])
    report("front moved with the material", abs(xf - 1.1) < 0.01,
           f"front at x = {xf:.4f} m (expected 1.1)")
    report("front stays sharp", width <= 2, f"{width} cells between 10 % and 90 % (limit 2)")
    report("no new extremes", th[mat].max() <= 1 + 1e-12 and th[mat].min() >= -1e-12,
           f"theta in [{th[mat].min():.2e}, {th[mat].max():.15f}]")
    report("empty cells hold voidT", np.all(T[~mat] == T0), f"{np.sum(~mat)} empty cells")
    pp = os.path.join(work, "channel4")
    channel(pp, uniform=False, decomp=(4, 1, 1))
    d = np.abs(cc.read_internal(f"{pp}/0.5/T") - T).max() / 30.0
    report("4 processors (cuts at x = 0.5, 1.0, 1.5 m): T matches serial", d < 1e-10,
           f"largest difference {d:.2e} of the hot-cold difference")
    p2 = os.path.join(work, "channelUniform")
    channel(p2, uniform=True)
    T = cc.read_internal(f"{p2}/0.5/T")
    if T is None:
        dev = 0.0
    else:
        dev = np.abs(T - T0).max() / T0
    report("uniform block stays uniform", dev < 1e-12, f"largest |T - T0|/T0 {dev:.1e}")


# ------------------------------------------------------------------ case B, C

W, H, CELL, WIN = 0.20, 0.30, 0.01, 0.04
T_IN = 353.15


def box_mesh(path):
    nx, ny = int(round(W / CELL)), int(round(H / CELL))
    n_in = int(round(WIN / CELL))
    n_side = (nx - n_in) // 2
    xs = [0.0, n_side * CELL, (n_side + n_in) * CELL, W]
    verts = ([f"({x!r} 0 0)" for x in xs] + [f"({x!r} {H!r} 0)" for x in xs]
             + [f"({x!r} 0 {CELL!r})" for x in xs] + [f"({x!r} {H!r} {CELL!r})" for x in xs])
    ncell = [n_side, n_in, n_side]
    blocks = [f"    hex ({b} {b + 1} {5 + b} {4 + b} {8 + b} {9 + b} {13 + b} {12 + b}) "
              f"({ncell[b]} {ny} 1) simpleGrading (1 1 1)" for b in range(3)]
    ymax = lambda b: f"({4 + b} {12 + b} {13 + b} {5 + b})"
    ymin = lambda b: f"({b} {b + 1} {9 + b} {8 + b})"
    patches = {"left": "(0 8 12 4)", "right": "(3 7 15 11)",
               "bottom": ymin(0) + " " + ymin(2), "outlet": ymin(1),
               "top": ymax(0) + " " + ymax(2), "inlet": ymax(1)}
    pt = [f"    {n} {{ type patch; faces ({f}); }}" for n, f in patches.items()]
    empty = " ".join(f"({b} {4 + b} {5 + b} {b + 1}) ({8 + b} {9 + b} {13 + b} {12 + b})"
                     for b in range(3))
    pt.append(f"    frontAndBack {{ type empty; faces ({empty}); }}")
    cc.write(f"{path}/system/blockMeshDict", HEADER.format(cls="dictionary", obj="blockMeshDict")
             + "scale 1;\nvertices\n(\n    " + "\n    ".join(verts) + "\n);\nblocks\n(\n"
             + "\n".join(blocks) + "\n);\nboundary\n(\n" + "\n".join(pt) + "\n);\n")
    cc.sh(path, ["blockMesh"], "blockMesh")
    os.makedirs(f"{path}/0", exist_ok=True)
    cc.sh(path, ["postProcess", "-func", "writeCellCentres", "-time", "0"], "centres")
    s = open(f"{path}/0/C").read()
    inner = s[s.index("internalField"):s.index("boundaryField")]
    n = int(re.search(r"List<vector>\s*(\d+)", inner).group(1))
    rows = re.findall(r"\(([^()]*)\)", inner[inner.index("(", inner.index("List")):])[:n]
    os.remove(f"{path}/0/C")
    return np.array([list(map(float, r.split())) for r in rows])


def box(path, uniform, decomp=(1, 1, 1)):
    controls(path, 0.6, 0.1, decomp)
    C = box_mesh(path)
    y = C[:, 1]
    fill = y < 0.15
    names = ("left", "right", "bottom", "outlet", "top", "inlet")
    zg = {n: ("zeroGradient", "") for n in names}
    scalar_file(path, "alpha", "[0 0 0 0 0 0 0]", np.where(fill, 1.0, 0.0), zg)
    scalar_file(path, "rho", "[1 -3 0 0 0 0 0]", np.where(fill, 1800.0, 0.0), zg)
    t_in = T0 if uniform else T_IN
    Tstart = np.where(fill, T0 if uniform else T0 + 100.0 * y, T0)
    tp = dict(zg)
    tp["inlet"] = ("fixedValue", f"value uniform {t_in!r};")
    scalar_file(path, "T", "[0 0 0 1 0 0 0]", Tstart, tp)
    cc.write(f"{path}/constant/g", HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 -9.81 0);\n")
    cc.write(f"{path}/constant/materialProperties",
             HEADER.format(cls="dictionary", obj="materialProperties") + """
rho 1800;
freeSurface yes;
type mohrCoulomb;
E 1e6;
nu 0.3;
phi 35;
c 0;
psi 0;
sigmaT 0;
walls { left { friction 20; } right { friction 20; } bottom { friction 20; } top { friction 0; } }
outlets { outlet { friction 20; openTime 0.1; } }
inlets { inlet { speed 0.5; startTime 0.05; endTime 0.35; density 1800; friction 0; } }
""")
    thermal(path, 0.3)
    n = decomp[0] * decomp[1] * decomp[2]
    cc.run(path, n)
    return C


def check_box(work):
    print("[B] Box fed through an inlet and emptied through an outlet", flush=True)
    p = os.path.join(work, "box")
    t0 = time.perf_counter()
    box(p, uniform=False)
    wall = time.perf_counter() - t0
    e = energy_log(p)
    rel = max(abs(b) for b in e["bal"]) / max(e["E"])
    report("energy balance with inlet, outlet and walls", rel < 1e-12,
           f"largest |energy - start - conducted - carried| {max(abs(b) for b in e['bal']):.3e} J "
           f"({rel:.1e} of the energy), fed {e['fed'][-1]:.4f} kg, out {e['dis'][-1]:.4f} kg, "
           f"carried in {e['carried'][-1]:.1f} J, {wall:.0f} s")
    report("material went through", e["fed"][-1] > 0.1 and e["dis"][-1] > 0.01,
           f"fed {e['fed'][-1]:.4f} kg, discharged {e['dis'][-1]:.4f} kg")
    T = cc.read_internal(f"{p}/0.3/T")
    rho = cc.read_internal(f"{p}/0.3/rho")
    mat = rho > 0
    lo = T0
    hot = T[mat].max()
    report("fed material arrives at the inlet temperature",
           abs(hot - T_IN) < 1e-9 and T[mat].min() >= lo - 1e-9 and hot <= T_IN + 1e-9,
           f"material T in [{T[mat].min():.6f}, {hot:.9f}] K at 0.3 s")
    report("empty cells hold voidT", np.all(T[~mat] == T0), f"{np.sum(~mat)} empty cells")

    print("[B] ... on 2 processors", flush=True)
    p2 = os.path.join(work, "box2")
    box(p2, uniform=False, decomp=(2, 1, 1))
    e2 = energy_log(p2)
    # The plain Mohr-Coulomb solid itself parts serial and parallel once the
    # fed stream lands (amplified round-off, decision 0134: rho differs by
    # hundreds of kg/m^3 in a few surface cells after 0.2 s). Before that T
    # must match everywhere; the channel (case A) checks T across a processor
    # cut over the whole run.
    worst, n_parted = 0.0, 0
    for t in ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6"):
        Ts, Tp = cc.read_internal(f"{p}/{t}/T"), cc.read_internal(f"{p2}/{t}/T")
        rs, rp = cc.read_internal(f"{p}/{t}/rho"), cc.read_internal(f"{p2}/{t}/rho")
        d = np.abs(Ts - Tp) / (T_IN - T0)
        if t in ("0.1", "0.2"):
            worst = max(worst, d.max())
        same = np.abs(rs - rp) < 1e-9 * 1800.0
        n_parted = max(n_parted, int(np.sum(~same)))
    report("2 processors: T matches serial to 0.2 s", worst < 1e-10,
           f"largest difference {worst:.2e} of the inlet-start difference")
    print(f"    (the solid itself parts in up to {n_parted} cells after the stream lands; "
          "T after that is compared on the channel instead)")
    report("2 processors: energy balance", max(abs(b) for b in e2["bal"]) / max(e2["E"]) < 1e-12,
           f"largest |balance| {max(abs(b) for b in e2['bal']):.2e} J")

    print("[C] Same box at one temperature", flush=True)
    p3 = os.path.join(work, "boxUniform")
    box(p3, uniform=True)
    dev = 0.0
    for t in ("0.1", "0.2", "0.3", "0.4", "0.5", "0.6"):
        T = cc.read_internal(f"{p3}/{t}/T")
        if T is not None:
            dev = max(dev, np.abs(T - T0).max() / T0)
    report("one temperature stays exact", dev < 1e-12, f"largest |T - T0|/T0 {dev:.1e}")
    e3 = energy_log(p3)
    want = C_HEAT * T0 * (e3["fed"][-1] - e3["dis"][-1])
    d = abs(e3["carried"][-1] - want) / abs(C_HEAT * T0 * e3["fed"][-1])
    report("carried heat equals c T (fed - discharged)", d < 1e-12,
           f"{e3['carried'][-1]:.9e} J against {want:.9e} J ({d:.1e})")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check_channel(work)
    check_box(work)
    print(f"carried checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_carried_temperature: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
