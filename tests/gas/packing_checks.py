"""WP F.8f: packing-dependent frictional pressure (materialProperties
princetonFr; decision 0133; Princeton / Srivastava-Sundaresan law as in
Breard, Dufek and Roche 2019):

    Pc = Fr (es - esMin)^r / (esMax - es)^s   (0 below esMin)

with es = rho/(alpha rhoS) the solid fraction of the filled part of a cell;
the effective mean stress P = -tr(sigma)/3 is kept at or below Pc.

Case: a dry column 0.3 m high in a box 4 x 40 cells of 10 mm, beads at
es 0.58 (Pc 520 Pa) with geostatic stress (K0 = nu/(1 - nu)), Mohr-Coulomb
phi 24, princetonFr 0.05, r 2, s 5, esMin 0.49, esMax 0.64, no damping.

L  At every write and in every material cell P <= Pc(es) within 1e-3 of Pc
   (the bound is applied after the Lagrange step; the remap then moves es
   by one step's advection).
C  The column compacts: the base cell's es is above 0.58 and at or above the
   packing whose Pc equals its final P (es only grows as it compacts: the
   law is a limit, unloading is elastic), and below esMax.
O  Without princetonFr the column keeps es 0.58 (option off changes nothing).
P  2 processors (cut across the column at y = 0.1 m) equal serial to 1e-10.

    python3 packing_checks.py WORKDIR
"""

import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fluid_cases as fl  # noqa: E402
import pore_cases as pc  # noqa: E402

NX, NY, CELL, HB = 4, 40, 0.01, 0.3
N = NX * NY
RHO_B, RHOS = 1450.0, 2500.0
FR, R, S, EMIN, EMAX = 0.05, 2.0, 5.0, 0.49, 0.64
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def pc_of(es):
    es = np.asarray(es, dtype=float)
    return np.where(es > EMIN, FR * np.maximum(es - EMIN, 0) ** R / np.maximum(EMAX - es, 1e-6) ** S, 0.0)


def make(path, princeton=True, decomp=(1, 1, 1)):
    pc.make_case(path, size=(NX * CELL, NY * CELL), cells=(NX, NY), e=1 - RHO_B / RHOS, d=7.5e-5,
                 end=0.5, write_interval=0.05, p_patches={"top": "type fixedValue; value uniform 0;"},
                 decomp=decomp, coupled=True)
    mat = "type mohrCoulomb;\nE 1e6;\nnu 0.3;\nphi 24;\nc 0;\npsi 0;\nrhoS 2500;\n"
    if princeton:
        mat += f"princetonFr {FR!r};\nprincetonR {R!r};\nprincetonS {S!r};\nprincetonMin {EMIN!r};\nprincetonMax {EMAX!r};\n"
    fl._common(path, rho=RHO_B, E=1e6, material=mat, extra_gas="",
               walls="    left { friction 0; }\n    right { friction 0; }\n    bottom { friction 0; }")
    os.remove(f"{path}/constant/gasProperties")
    alpha, rho, sig = [], [], []
    for j in range(NY):
        w = RHO_B * 9.81 * max(HB - (j + 0.5) * CELL, 0.0)
        for _ in range(NX):
            full = (j + 0.5) * CELL < HB
            alpha.append(1.0 if full else 0.0)
            rho.append(RHO_B if full else 0.0)
            sig.append((-0.3 / 0.7 * w, -w))
    fl._write_fields(path, alpha, rho, sig, N)
    if decomp != (1, 1, 1):
        pc.write(f"{path}/system/decomposeParDict", pc.HEADER.format(cls="dictionary", obj="decomposeParDict")
                 + "numberOfSubdomains 2;\nmethod simple;\ncoeffs { n (1 2 1); }\n")


def mean_stress(path, t):
    txt = open(f"{path}/{t}/sigma").read().split("internalField")[1].split("boundaryField")[0]
    if txt.split()[0] == "uniform":
        return np.zeros(N)
    s = np.array([[float(x) for x in v.split()] for v in re.findall(r"\(([^()]+)\)", txt)][:N])
    return -(s[:, 0] + s[:, 3] + s[:, 5]) / 3


def state(path, t):
    a = pc.field(path, t, "alpha", N)
    es = np.where(a > 1e-6, pc.field(path, t, "rho", N) / np.maximum(a, 1e-12) / RHOS, 0.0)
    return a, es, mean_stress(path, t)


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    print(f"Column {HB} m, {NX} x {NY} cells; es 0.58, Pc(0.58) = {float(pc_of(0.58)):.0f} Pa; "
          f"base geostatic P = {RHO_B * 9.81 * (HB - CELL / 2) * (1 + 2 * 0.3 / 0.7) / 3:.0f} Pa")
    a1 = f"{work}/on"
    make(a1)
    pc.run(a1, 1)
    worst = 0.0
    for t in pc.times(a1):
        a, es, P = state(a1, t)
        m = a > 1e-3
        worst = max(worst, float(np.max((P[m] - pc_of(es[m])) / np.maximum(pc_of(es[m]), 1.0))))
    report("L P <= Pc(es) at every write", worst <= 1e-3, f"largest (P - Pc)/Pc {worst:.1e}")
    a, es, P = state(a1, pc.times(a1)[-1])
    base = es[:NX].mean()
    grid = np.linspace(0.58, EMAX - 1e-4, 20001)
    need = float(grid[np.searchsorted(pc_of(grid), P[:NX].mean())])
    report("C column compacts", 0.58 < base < EMAX and base >= need - 1e-4,
           f"base es {base:.4f} (packing carrying its final P {P[:NX].mean():.0f} Pa: {need:.4f}); "
           f"top material cell es {es[(int(HB / CELL) - 1) * NX]:.4f}")
    a2 = f"{work}/off"
    make(a2, princeton=False)
    pc.run(a2, 1)
    _, es2, _ = state(a2, pc.times(a2)[-1])
    dev = float(np.abs(es2[es2 > 0] - 0.58).max())
    report("O option off keeps the packing", dev <= 2e-3, f"largest |es - 0.58| {dev:.1e} (elastic, <= 2e-3)")
    a3 = f"{work}/on2"
    make(a3, decomp=(1, 2, 1))
    pc.run(a3, 2)
    d = max(float(np.abs(pc.field(a3, t, f, N) - pc.field(a1, t, f, N)).max()
                  / max(np.abs(pc.field(a1, t, f, N)).max(), 1e-30))
            for t in pc.times(a1) for f in ("rho", "alpha"))
    report("P 2 processors", d <= 1e-10, f"rho, alpha within {d:.1e} of serial")
    print(f"\ntest_packing_pressure: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)}), "
          f"{time.perf_counter() - t0:.0f} s")
    return 0 if all(res) else 1


if __name__ == "__main__":
    sys.exit(main())
