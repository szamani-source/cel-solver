"""WP F.2 done-when: pore pressure in a rigid bed (decisions 0125, 0127).

The solid is frozen (no gravity, a uniform bed at rest), so only the gas
moves. Every run uses the ergunWenYu drag law. Cases (2D, one cell deep):

A  Permeameter: column 0.2 m tall, 40 cells, quartz grains d 1 mm, gas
   fraction 0.40; air blown in at the bottom (inflow) at 0.05 and 0.5 m/s
   (superficial, at p0), top open (p = 0 gauge), sides sealed; run to
   0.03 s (steady after about 0.003 s). The pressure at the bottom face is
   checked against Ergun integrated for a compressible isothermal gas,
   p_in^2 = p_top^2 + 2 R T (a m + b m^2) L with m the mass flux.
W  The same column loosened to gas fraction 0.85 at 0.5 m/s: the Wen-Yu
   branch of the law and its root finder; reference: the Wen-Yu pressure
   gradient integrated down the column (tools/gas_reference.py law).
B  Step at a face: column 0.5 m, 500 cells, d 0.1 mm, gas fraction 0.4;
   the bottom face is raised to 0.001 p0 at t = 0, the top is sealed.
   Profiles at 0.0025, 0.005 and 0.01 s against
   dp erfc(x/(2 sqrt(c t))), c = k p0/(mu e).
C  Sealed box: 0.1 x 0.1 m, 20 x 20 cells, sheared (top edge shifted by
   half the width, so the faces are 27 deg non-orthogonal), d 1 mm; a
   20 kPa blob in one corner, every patch sealed; run to 0.01 s. Gas mass
   conserved to 1e-12; the pressure ends uniform at the mass-conserving
   value. Large gradients, so Ergun's inertial term is strong.

A, B and C also run on 2 and 4 processors (cuts through the column, the
front and the blob): fields match serial to 1e-10 of the pressure scale,
gas mass is conserved to 1e-12 on each. The cost of gas on against gas
off is measured on a 100 x 100 box (report only; F.9 sets the limit).

    python3 pore_pressure_checks.py WORKDIR
"""

import math
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

RT = gr.R_AIR * pc.T
MU, RG = gr.air_viscosity(pc.T), gr.air_density(pc.P0, pc.T)

res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def ergun_inlet(U, e, d, L):
    """Gauge inlet pressure, Ergun with compressible isothermal gas."""
    m = U * RG
    a = 150.0 * MU * (1 - e) ** 2 / (e ** 3 * d * d)
    b = 1.75 * (1 - e) / (e ** 3 * d)
    return math.sqrt(pc.P0 ** 2 + 2 * RT * (a * m + b * m * m) * L) - pc.P0


def wen_yu_inlet(U, e, d, L, n=4000):
    """Gauge inlet pressure from the Ergun-Wen-Yu law, RK4 down the column."""
    m = U * RG

    def dpdx(p):            # -dp/dx at absolute pressure p
        rho = p / RT
        u = m / (rho * e)   # interstitial slip
        return gr.drag_ergun_wen_yu(e, u, d, MU, rho) * u / e

    p, h = pc.P0, L / n
    for _ in range(n):
        k1 = dpdx(p)
        k2 = dpdx(p + 0.5 * h * k1)
        k3 = dpdx(p + 0.5 * h * k2)
        k4 = dpdx(p + h * k3)
        p += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return p - pc.P0


def drift(path):
    M, N, ns, r = pc.mass_log(path)
    return np.abs(M - N - (M[0] - N[0])).max() / M[0], ns.mean(), r.max()


def column(path, U, e, d, nprocs):
    L, n = 0.2, 40
    pc.make_case(path, size=(L / n, L), cells=(1, n), e=e, d=d, end=0.03, write_interval=0.01,
                 p_patches={"top": "type fixedValue; value uniform 0;",
                            "bottom": "type fixedGradient; gradient uniform 0;"},
                 inflow={"bottom": U}, decomp=(1, nprocs, 1))
    secs = pc.run(path, nprocs)
    t = pc.times(path)[-1]
    p = pc.field(path, t, "p", n)
    pin = p[0] + pc.patch_gradient(path, t, "p", "bottom") * (L / n) / 2
    return secs, p, pin


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    print(f"Air at {pc.T} K, p0 {pc.P0} Pa: rho_g {RG:.6f} kg/m^3, mu {MU:.6e} Pa s; "
          "drag ergunWenYu; solid frozen\n", flush=True)

    # ---------------------------------------------------------------- A, W
    print("A  permeameter (column 0.2 m, 40 cells, d 1 mm, e 0.40)")
    par = {}
    for U in (0.05, 0.5):
        ref = ergun_inlet(U, 0.4, 1e-3, 0.2)
        inc = (150 * MU * 0.36 / (0.064 * 1e-6) * U + 1.75 * RG * 0.6 / (0.064 * 1e-3) * U * U) * 0.2
        out = {}
        for n in (1, 2, 4):
            path = f"{work}/A_{U}_{n}"
            secs, p, pin = column(path, U, 0.4, 1e-3, n)
            dr, ns, _ = drift(path)
            out[n] = (secs, p, pin, dr, ns)
        secs, p, pin, dr, ns = out[1]
        err = pin / ref - 1
        report(f"A U = {U} m/s", abs(err) <= 0.01,
               f"inlet {pin:.3f} Pa against Ergun {ref:.3f} Pa (compressible; incompressible "
               f"{inc:.1f}): {100 * err:+.4f} % (<= 1 %); {ns:.1f} solves per step")
        par[f"A {U}"] = out
    # Wen-Yu branch
    ref = wen_yu_inlet(0.5, 0.85, 1e-3, 0.2)
    secs, p, pin = column(f"{work}/W", 0.5, 0.85, 1e-3, 1)
    err = pin / ref - 1
    report("W Wen-Yu branch e = 0.85, U = 0.5 m/s", abs(err) <= 0.01,
           f"inlet {pin:.4f} Pa against the integrated law {ref:.4f} Pa: {100 * err:+.4f} % (<= 1 %)")

    # ------------------------------------------------------------------- B
    print("\nB  step at a face (column 0.5 m, 500 cells, d 0.1 mm, e 0.40, step 0.001 p0)")
    L, n, e, d = 0.5, 500, 0.4, 1e-4
    dp = 0.001 * pc.P0
    c = gr.diffusivity_rigid(d, e, MU, pc.P0)
    y = (np.arange(n) + 0.5) * L / n
    out = {}
    for k in (1, 2, 4):
        path = f"{work}/B_{k}"
        pc.make_case(path, size=(L / n, L), cells=(1, n), e=e, d=d, end=0.01, write_interval=0.0025,
                     p_patches={"bottom": f"type fixedValue; value uniform {dp!r};"},
                     decomp=(1, k, 1))
        secs = pc.run(path, k)
        dr, ns, _ = drift(path)
        out[k] = (secs, None, None, dr, ns)
    par["B"] = out
    worst = 0.0
    parts = []
    for t in ("0.0025", "0.005", "0.01"):
        p = pc.field(f"{work}/B_1", t, "p", n)
        ref = dp * np.array([math.erfc(yy / (2 * math.sqrt(c * float(t)))) for yy in y])
        err = np.abs(p - ref).max() / dp
        worst = max(worst, err)
        parts.append(f"t {t} s {100 * err:.3f} %")
    report("B erfc profiles", worst <= 0.02,
           f"c = {c:.4f} m^2/s; largest error / step: " + ", ".join(parts) + " (<= 2 %)")

    # ------------------------------------------------------------------- C
    print("\nC  sealed sheared box (0.1 m, 20 x 20 cells, 27 deg non-orthogonal, d 1 mm, 20 kPa blob)")
    out = {}
    for k, dec in ((1, (1, 1, 1)), (2, (2, 1, 1)), (4, (2, 2, 1))):
        path = f"{work}/C_{k}"
        pc.make_case(path, size=(0.1, 0.1), cells=(20, 20), e=0.4, d=1e-3, end=0.01,
                     write_interval=0.0005, skew=0.05, decomp=dec, p_patches={},
                     p_set=((-1, -1), (0.04, 0.04), 20000.0))
        secs = pc.run(path, k, set_fields=True)
        dr, ns, _ = drift(path)
        out[k] = (secs, None, None, dr, ns)
    par["C"] = out
    p0f = pc.field(f"{work}/C_1", "0", "p", 400)
    pend = pc.field(f"{work}/C_1", pc.times(f"{work}/C_1")[-1], "p", 400)
    peq = p0f.mean()        # equal cell volumes, uniform packing
    spread = np.abs(pend - peq).max() / 20000.0
    report("C sealed box mass", out[1][3] <= 1e-12,
           f"gas mass drift {out[1][3]:.1e} of the mass (<= 1e-12) over the run; "
           f"{out[1][4]:.1f} solves per step")
    report("C sealed box equilibrium", spread <= 1e-6,
           f"end pressure {pend.min():.4f} to {pend.max():.4f} Pa against the mass-conserving "
           f"{peq:.4f} Pa: {spread:.1e} of the blob (<= 1e-6)")

    # ------------------------------------------------------------ parallel
    print("\nParallel (2 and 4 processors against serial)")
    worst_f, worst_m, rows = 0.0, 0.0, []
    for name, out in par.items():
        base = name.replace(" ", "_")
        for k in (2, 4):
            ser, pa = f"{work}/{base}_1", f"{work}/{base}_{k}"
            if name.startswith("A"):
                ser, pa = f"{work}/A_{name.split()[1]}_1", f"{work}/A_{name.split()[1]}_{k}"
            scale, diff, qscale, qdiff = 0.0, 0.0, 0.0, 0.0
            for t in pc.times(ser):
                a, b = pc.field(ser, t, "p"), pc.field(pa, t, "p")
                scale = max(scale, np.abs(a).max())
                diff = max(diff, np.abs(a - b).max())
                a, b = pc.field(ser, t, "qGas"), pc.field(pa, t, "qGas")
                qscale = max(qscale, np.abs(a).max())
                qdiff = max(qdiff, np.abs(a - b).max())
            rel, qrel = diff / scale, qdiff / max(qscale, 1e-300)
            worst_f = max(worst_f, rel, qrel)
            worst_m = max(worst_m, out[k][3], out[1][3])
            rows.append(f"{name} np {k}: p {rel:.0e}, qGas {qrel:.0e}")
    report("parallel fields", worst_f <= 1e-10,
           "; ".join(rows) + f" (largest difference at any write over the largest value of the run; worst {worst_f:.0e}, <= 1e-10)")
    report("parallel gas mass", worst_m <= 1e-12,
           f"worst mass drift over all runs and processor counts {worst_m:.1e} (<= 1e-12)")
    print("  run times [s] (celFoam only): " + "; ".join(
        f"{name} " + "/".join(f"{out[k][0]:.1f}" for k in (1, 2, 4)) for name, out in par.items()))

    # ---------------------------------------------------------------- cost
    print("\nCost: 100 x 100 sealed box, 300 steps, gas off and on (report)")
    secs = {}
    for tag in ("off", "on"):
        path = f"{work}/cost_{tag}"
        pc.make_case(path, size=(0.2, 0.2), cells=(100, 100), e=0.4, d=1e-3, end=1.0,
                     write_interval=1.0, p_patches={}, p_set=((-1, -1), (0.08, 0.08), 20000.0))
        # 300 steps exactly: fixed end time from the stable step
        dt = 0.4 * 0.002 / math.sqrt(1e6 * 0.7 / (1.3 * 0.4) / (0.6 * pc.RHOP))
        txt = open(f"{path}/system/controlDict").read()
        txt = txt.replace("endTime         1.0;", f"endTime         {300 * dt * 0.999!r};")
        open(f"{path}/system/controlDict", "w").write(txt)
        if tag == "off":
            os.remove(f"{path}/constant/gasProperties")
        secs[tag] = pc.run(path, 1, set_fields=(tag == "on"))
    steps = {t: open(f"{work}/cost_{t}/log.celFoam").read().count("\nTime = ") for t in secs}
    per = {t: secs[t] / max(steps[t], 1) for t in secs}
    print(f"  gas off {1e3 * per['off']:.2f} ms per step, gas on {1e3 * per['on']:.2f} ms per step "
          f"({steps['on']} steps, 10,000 cells): ratio {per['on'] / per['off']:.2f}")

    ok = all(res)
    print(f"\nPore pressure checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
