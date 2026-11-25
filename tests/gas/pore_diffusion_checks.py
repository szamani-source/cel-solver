"""WP F.8e: pore pressure as an advected-diffused field (gasProperties
diffusivity; decision 0133, route B after Iverson & Denlinger 2001,
Gueugneau et al. 2017, Aravena et al. 2021).

The solid is frozen (no gravity, a uniform bed at rest), so only the
diffusion acts. Every run uses the ergunWenYu drag law.

D  Step at a face: column 0.5 m, 250 cells, d 75 um, gas fraction 0.42;
   the bottom face is raised to 0.001 p0 at t = 0, the top is sealed.
   With diffusivity kappa = 0.015 and 0.05 m^2/s the profiles at 0.25 and 0.5 s
   match dp erfc(x/(2 sqrt(kappa t))) within 2 % of the step (the
   diffusivity is the given one, not the rigid-bed k p0/(mu e)).
P  The kappa = 0.015 column on 2 and 4 processors (cuts through the front):
   p within 1e-10 of the step of serial at every write.

The advection (pressure carried by the grains, no compression or suction)
is exercised by the dam break (tests/manual/test_fluidised_collapse.sh).

    python3 pore_diffusion_checks.py WORKDIR
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def main():
    work = sys.argv[1]
    t_all = time.perf_counter()
    L, n, e, d = 0.5, 250, 0.42, 7.5e-5
    dp = 0.001 * pc.P0
    y = (np.arange(n) + 0.5) * L / n
    times = ("0.25", "0.5")
    c_rigid = gr.diffusivity_rigid(d, e, gr.air_viscosity(pc.T), pc.P0)
    print(f"D  step at a face (column {L} m, {n} cells, d {d * 1e6:g} um, e {e}, step 0.001 p0); "
          f"rigid-bed k p0/(mu e) = {c_rigid:.4f} m^2/s")
    for kap in (0.015, 0.05):
        out = {}
        for k in ((1, 2, 4) if kap == 0.015 else (1,)):
            path = f"{work}/D_{kap}_{k}"
            pc.make_case(path, size=(L / n, L), cells=(1, n), e=e, d=d, end=0.5, write_interval=0.25,
                         p_patches={"bottom": f"type fixedValue; value uniform {dp!r};"},
                         decomp=(1, k, 1), extra_gas=f"diffusivity {kap!r};\n")
            out[k] = pc.run(path, k)
        worst, parts = 0.0, []
        for t in times:
            p = pc.field(f"{work}/D_{kap}_1", t, "p", n)
            ref = dp * np.array([math.erfc(yy / (2 * math.sqrt(kap * float(t)))) for yy in y])
            err = np.abs(p - ref).max() / dp
            worst = max(worst, err)
            parts.append(f"t {t} s {100 * err:.3f} %")
        report(f"D kappa = {kap} m^2/s erfc profiles", worst <= 0.02,
               "largest error / step: " + ", ".join(parts) + f" (<= 2 %); serial {out[1]:.0f} s")
        if kap == 0.015:
            for k in (2, 4):
                diff = max(np.abs(pc.field(f"{work}/D_{kap}_{k}", t, "p", n)
                                  - pc.field(f"{work}/D_{kap}_1", t, "p", n)).max() / dp for t in times)
                report(f"P {k} processors", diff <= 1e-10, f"p within {diff:.1e} of the step of serial (<= 1e-10)")
    print(f"\ntest_pore_diffusion: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)}), "
          f"{time.perf_counter() - t_all:.0f} s")
    return 0 if all(res) else 1


if __name__ == "__main__":
    sys.exit(main())
