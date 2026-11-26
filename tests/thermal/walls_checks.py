"""WP T.3 done-when: convective and heat-flux walls, Zehner-Schluender bed
conductivity (decision 0135).

Cases (celFoam, an elastic bed at rest without gravity, rho 1500, c 1000):
A  Slab 2L = 0.1 m (50 cells of 2 mm), k 1500 W/(m K) (a = 1e-3 m^2/s),
   both faces convective (thermalProperties walls { h; Tinf; }) at
   Bi = h L/k = 1 and 10, from 100 K above Tinf. Against the Carslaw-Jaeger
   series of tools/thermal_reference.py at Fo = 0.2, 0.4, 0.6, 0.8, 1.0:
   largest error under 1 % of the start difference. Energy balance to 1e-12.
   Bi 10 also on 2 processors (cut at the middle): matches serial to 1e-10.
B  The same slab heated through one face by walls { left { q 5000; } }:
   heat in = q A t to 1e-12.
C  Column 0.2 m (20 cells) with freeSurface: packing rho 1500 to 1800 kg/m^3
   in rows, a half-filled row (alpha 0.5) and void above; conductivity
   zehnerSchlunder (kParticle 3, kGas 0.02569, rhoP 2650). The written kBed
   equals tools/thermal_reference.py zs_bed at the filled part's gas fraction
   1 - rho/(alpha rhoP) to 1e-12 in every cell, 0 in void.

    python3 walls_checks.py WORKDIR
"""

import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import conduction_checks as cc  # noqa: E402
import carried_checks as ca  # noqa: E402
import thermal_reference as tr  # noqa: E402

HEADER = cc.HEADER
report = cc.report
failed = cc.failed
RHO, C_HEAT, K = 1500.0, 1000.0, 1500.0
A = K / (RHO * C_HEAT)
L, N = 0.05, 50
TINF, DT0 = 293.15, 100.0


def slab(path, walls, T_start, decomp=(1, 1, 1), end=2.5, write=0.5):
    cc.common(path, end=end, write_interval=write, decomp=decomp)
    cc.blockmesh(path, 2 * L, 2 * L / N, N, 1)
    cc.write_T(path, T_start, {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})
    tp = f"{path}/constant/thermalProperties"
    open(tp, "a").write("walls\n{\n" + walls + "}\n")
    cc.run(path, decomp[0] * decomp[1] * decomp[2])


def check_slab(work):
    print("[A] Slab cooled by convection on both faces", flush=True)
    x = (np.arange(N) + 0.5) * (2 * L / N) - L
    res = {}
    for Bi in (1.0, 10.0):
        h = Bi * K / L
        p = os.path.join(work, f"slab{Bi:g}")
        w = f"    left {{ h {h!r}; Tinf {TINF!r}; }}\n    right {{ h {h!r}; Tinf {TINF!r}; }}\n"
        t0 = time.perf_counter()
        slab(p, w, TINF + DT0)
        worst = 0.0
        for t in (0.5, 1.0, 1.5, 2.0, 2.5):
            Fo = A * t / L ** 2
            T = cc.read_internal(f"{p}/{t:g}/T")
            th = (T - TINF) / DT0
            ex = np.array([tr.slab_theta(abs(xi) / L, Fo, Bi) for xi in x])
            worst = max(worst, np.abs(th - ex).max())
        report(f"Bi {Bi:g} against Carslaw-Jaeger, Fo 0.2 to 1", worst < 0.01,
               f"largest error {worst * 100:.3f} % of the start difference (limit 1 %), "
               f"{time.perf_counter() - t0:.0f} s")
        log = ca.energy_log(p)
        rel = max(abs(b) for b in log["bal"]) / log["E"][0]
        report(f"Bi {Bi:g} energy balance", rel < 1e-12,
               f"largest |balance| {max(abs(b) for b in log['bal']):.2e} J ({rel:.1e}), "
               f"heat in {log['heat'][-1]:.4e} J")
        res[Bi] = (p, w)
    p, w = res[10.0]
    p2 = os.path.join(work, "slab10par")
    slab(p2, w, TINF + DT0, decomp=(2, 1, 1))
    d = max(np.abs(cc.read_internal(f"{p2}/{t:g}/T") - cc.read_internal(f"{p}/{t:g}/T")).max()
            for t in (0.5, 2.5)) / DT0
    report("Bi 10 on 2 processors matches serial", d < 1e-10, f"largest difference {d:.2e}")

    print("[B] Heat flux through one face", flush=True)
    p3 = os.path.join(work, "flux")
    slab(p3, "    left { q 5000; }\n", TINF, end=1.0, write=1.0)
    q = ca.energy_log(p3)["heat"][-1]
    want = 5000.0 * (2 * L / N) ** 2 * 1.0
    report("heat in equals q A t", abs(q / want - 1) < 1e-12, f"{q:.15e} J against {want:.15e} J")


def check_zs(work):
    print("[C] Zehner-Schluender bed conductivity", flush=True)
    p = os.path.join(work, "zs")
    ny, H, w = 20, 0.2, 0.01
    cc.common(p, end=0.002, write_interval=0.001)
    cc.blockmesh(p, w, H, 1, ny)
    y = (np.arange(ny) + 0.5) * (H / ny)
    rho_m = np.linspace(1500.0, 1800.0, ny)
    alpha = np.where(y < 0.15, 1.0, np.where(y < 0.16, 0.5, 0.0))
    rho = alpha * rho_m
    zg = {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")}
    ca.scalar_file(p, "alpha", "[0 0 0 0 0 0 0]", alpha, zg)
    ca.scalar_file(p, "rho", "[1 -3 0 0 0 0 0]", rho, zg)
    cc.write_T(p, TINF, zg)
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read()
    assert "rho 1500.0;" in txt
    open(mp, "w").write(txt.replace("rho 1500.0;", "rho 1800.0;\nfreeSurface yes;"))
    cc.write(f"{p}/constant/thermalProperties", HEADER.format(cls="dictionary", obj="thermalProperties") + """
model        equilibrium;
c            1000;
conductivity zehnerSchlunder;
kParticle    3;
kGas         0.02569;
rhoP         2650;
""")
    cc.run(p)
    kb = cc.read_internal(f"{p}/0.002/kBed")
    rho2 = cc.read_internal(f"{p}/0.002/rho")
    al2 = cc.read_internal(f"{p}/0.002/alpha")
    worst = 0.0
    for i in range(ny):
        if rho2[i] <= 0:
            worst = max(worst, abs(kb[i]))
            continue
        e = 1 - rho2[i] / (al2[i] * 2650.0)
        worst = max(worst, abs(kb[i] / tr.zs_bed(e, 0.02569, 3.0) - 1))
    report("kBed equals zs_bed of the filled part", worst < 1e-12,
           f"largest relative difference {worst:.1e}; k from {kb[rho2 > 0].min():.4f} to "
           f"{kb.max():.4f} W/(m K), void {kb[rho2 <= 0].max() if np.any(rho2 <= 0) else 0}")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check_slab(work)
    check_zs(work)
    print(f"walls checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_thermal_walls: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
