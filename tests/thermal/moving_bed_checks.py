"""WP T.6 done-when: a moving bed past heated walls (decision 0135).

The case of Adapa et al. (2024, arXiv 2403.19892): ceramic particles (CARBO
HSP 40/70, d 404 microns) flow in plug flow down a 5 mm gap between plates
heated over 500 mm. Their measured flowing-bed inputs at 480 C and 10 mm/s
(Table 5): k_eff 0.28 W/(m K), near-wall gap D_air 35 microns; flowing bulk
density 2090 kg/m^3 (their frozen/tapped bed). The specific heat is not in
the paper (it cites Siegel et al. 2015, not read): 1100 J/(kg K) is used,
and the answer hardly depends on it (900 or 1300 move h by 1 %).

Model: half the gap (b 2.5 mm, 10 cells across; the mid-plane is a
smooth adiabatic wall), 500 mm long (200 cells), material fed through the
inlet and leaving through an open outlet, elastic, no gravity, no friction:
plug flow. The wall takes a uniform flux q 5000 W/m^2 (thermalProperties
walls { wall { q; } }). To keep the run short the bed moves at 1 m/s with
c scaled by 0.01/1, so the Peclet number U D_h rho c/k (821) and every
dimensionless result are those of the experiment; t = 0.7 s (1.4 transits).
The near-wall gap is the series resistance R_nw = D_air/k_air, as in their
model (eq. 2): the wall is T_bed,wall + q R_nw.

The remap carries T with the default limiter MC (superbee turns the smooth
profile into steps: local Nu off by 20 %).

A  Local Nu_Dh at x = 0.1, 0.25 and 0.45 m and the mean h over the heated
   length (q over the length-averaged wall minus mean bed temperature) within
   5 % of the slug-flow solution with the gap (Adapa et al. eqs. 4 and 9,
   tools/thermal_reference.py).
B  The mean h lies in the 225 to 350 W/(m^2 K) Adapa et al. give for a 5 mm,
   500 mm plate heat exchanger (their model, fitted to their measurements).
C  Energy balance (thermal energy - start - heat in - heat carried in) under
   1e-12 of the thermal energy; heat in = q A t to 1e-4 (the inlet patch
   conducts too: about 5e-5 of it).

    python3 moving_bed_checks.py WORKDIR
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import conduction_checks as cc  # noqa: E402
import carried_checks as ca  # noqa: E402
import thermal_reference as tr  # noqa: E402

report = cc.report
failed = cc.failed

B, LEN = 2.5e-3, 0.5
NX, NY = 200, 10
DY = B / NY
K_EFF, D_AIR, T_BED = 0.28, 35e-6, 753.15       # Adapa et al. Table 5, 480 C, 10 mm/s
RHO, CP = 2090.0, 1100.0
U_REAL, U_SIM = 0.010, 1.0
C_SIM = CP * U_REAL / U_SIM
Q = 5000.0
END = 0.7
K_AIR = tr.air_conductivity(T_BED)
R_NW = D_AIR / K_AIR
DH = 4 * B
PE = U_REAL * DH * RHO * CP / K_EFF
GAP = R_NW * K_EFF / (4 * B)


def make(path, decomp=(1, 1, 1)):
    ca.controls(path, END, 0.1, decomp)
    dz = DY
    cc.write(f"{path}/system/blockMeshDict", cc.HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices ((0 0 0) ({LEN!r} 0 0) ({LEN!r} {B!r} 0) (0 {B!r} 0)
          (0 0 {dz!r}) ({LEN!r} 0 {dz!r}) ({LEN!r} {B!r} {dz!r}) (0 {B!r} {dz!r}));
blocks (hex (0 1 2 3 4 5 6 7) ({NX} {NY} 1) simpleGrading (1 1 1));
boundary
(
    inlet {{ type patch; faces ((0 4 7 3)); }}
    outlet {{ type patch; faces ((1 2 6 5)); }}
    wall {{ type patch; faces ((0 1 5 4)); }}
    mid {{ type patch; faces ((3 7 6 2)); }}
    frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }}
);
""")
    cc.sh(path, ["blockMesh"], "blockMesh")
    n = NX * NY
    names = ("inlet", "outlet", "wall", "mid")
    zg = {k: ("zeroGradient", "") for k in names}
    ca.scalar_file(path, "alpha", "[0 0 0 0 0 0 0]", 1.0, zg)
    ca.scalar_file(path, "rho", "[1 -3 0 0 0 0 0]", RHO, zg)
    ca.vector_file(path, "U", np.tile([U_SIM, 0.0, 0.0], (n, 1)), zg)
    tp = dict(zg, inlet=("fixedValue", f"value uniform {T_BED!r};"))
    ca.scalar_file(path, "T", "[0 0 0 1 0 0 0]", T_BED, tp)
    cc.write(f"{path}/constant/g", cc.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n")
    cc.write(f"{path}/constant/materialProperties",
             cc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO!r};
freeSurface yes;
type linearElastic;
E 1e4;
nu 0.3;
walls {{ wall {{ friction 0; }} mid {{ friction 0; }} }}
outlets {{ outlet {{ friction 0; openTime 0; }} }}
inlets {{ inlet {{ speed {U_SIM!r}; startTime 0; endTime 100; density {RHO!r}; friction 0; }} }}
""")
    cc.write(f"{path}/constant/thermalProperties",
             cc.HEADER.format(cls="dictionary", obj="thermalProperties") + f"""
model        equilibrium;
c            {C_SIM!r};
conductivity constant;
k            {K_EFF!r};
voidT        {T_BED!r};
walls {{ wall {{ q {Q!r}; }} }}
""")


def metrics(p, verbose=False):
    """Largest local Nu error at x = 0.1, 0.25, 0.45 m, the mean h and the
    mean h of eq. 9."""
    T = cc.read_internal(f"{p}/{END:g}/T").reshape(NY, NX)
    x = (np.arange(NX) + 0.5) * LEN / NX
    Tm = T.mean(axis=0)
    Tw = T[0] + Q * (DY / 2) / K_EFF + Q * R_NW
    dT = Tw - Tm
    worst = 0.0
    for xs in (0.1, 0.25, 0.45):
        i = int(np.argmin(np.abs(x - xs)))
        nu = Q / dT[i] * DH / K_EFF
        want = tr.slug_plates_nu_local(DH / x[i] * PE, GAP)
        worst = max(worst, abs(nu / want - 1))
        if verbose:
            print(f"    x {x[i]:.4f} m: Nu {nu:.3f} against {want:.3f}", flush=True)
    h = Q / dT.mean()
    hw = tr.slug_plates_h_mean(K_EFF, D_AIR, K_AIR, RHO, CP, U_REAL, B, LEN)
    return worst, h, hw


def check(work):
    print("[A] Moving bed past a heated wall against the slug-flow solution", flush=True)
    print(f"    Pe {PE:.1f}, k_air {K_AIR:.4f} W/(m K), R_nw k/(4b) {GAP:.4f}, "
          f"c scaled to {C_SIM:g} J/(kg K)", flush=True)
    p = os.path.join(work, "bed")
    make(p)
    t0 = time.perf_counter()
    cc.run(p)
    secs = time.perf_counter() - t0
    worst, h, hw = metrics(p, verbose=True)
    rho = cc.read_internal(f"{p}/{END:g}/rho")
    rho = np.full(NX * NY, RHO) if rho is None else rho
    report("local Nu against the slug-flow solution", worst < 0.05,
           f"largest difference {worst * 100:.2f} % (limit 5 %), {secs:.0f} s")
    report("mean h against eq. 9", abs(h / hw - 1) < 0.05,
           f"{h:.1f} against {hw:.1f} W/(m^2 K) ({(h / hw - 1) * 100:+.2f} %, limit 5 %); "
           f"density {rho.min():.2f} to {rho.max():.2f} kg/m^3")
    print("[B] Against Adapa et al. (2024)", flush=True)
    report("mean h within 225 to 350 W/(m^2 K)", 225 <= h <= 350, f"{h:.1f} W/(m^2 K)")
    print("[C] Energy", flush=True)
    log = ca.energy_log(p)
    rel = max(abs(b) for b in log["bal"]) / log["E"][0]
    report("energy balance", rel < 1e-12, f"largest |balance| {rel:.1e} of the thermal energy")
    want = Q * LEN * DY * END
    # The inlet patch (fixedValue T) also conducts: the first cell warms by
    # about 0.2 K and loses about 5e-5 of q A t back through it
    report("heat in equals q A t", abs(log["heat"][-1] / want - 1) < 1e-4,
           f"{log['heat'][-1]:.15e} J against {want:.15e} J")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check(work)
    print(f"moving bed checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_moving_bed_heat: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
