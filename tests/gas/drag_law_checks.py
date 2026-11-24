"""WP F.7 done-when part 1: the gas-grain drag laws (decisions 0125 item 7, 0132).

A  Each law of celFluid (celDragTable prints K) against its closed form in
   tools/gas_reference.py to 1e-10 relative, gas fraction 0.35 to 0.999 and
   particle Reynolds number rho_g u d/mu from 1e-3 to 1e3, plus u = 0.
B  Limits: ergunWenYu is ergun below 0.8 and wenYu above (1e-14); at e -> 1
   and u -> 0 every law but ergun gives one sphere's Stokes drag
   18 mu (1-e)/d^2 (Dallavalle's Cd, used by diFelice and syamlalOBrien,
   gives 23.04/24 of it); beetstra's packed bed at rest is Carman-Kozeny
   (180 (1-e)^2 mu/(e^3 d^2) per superficial velocity, to the 2nd term).
C  The solver uses the law it is given: a permeameter (frozen bed, e 0.40,
   1 mm grains, 0.2 m, 40 cells, air blown in at 0.3 m/s) with each law; its
   inlet pressure against the law integrated down the column (RK4,
   isothermal gas) within 1e-4.

Every other test uses ergunWenYu (Sia, 2026-11-24).
"""

import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

LAWS = ["ergunWenYu", "ergun", "wenYu", "diFelice", "syamlalOBrien", "beetstra"]
RT = gr.R_AIR * pc.T
MU, RG = gr.air_viscosity(pc.T), gr.air_density(pc.P0, pc.T)
D = 1.0e-3

res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def table(law, pairs, d=D, mu=MU, rho=RG):
    inp = "\n".join(f"{float(e)!r} {float(u)!r}" for e, u in pairs) + "\n"
    out = subprocess.run(["celDragTable", law, repr(d), repr(mu), repr(rho)], input=inp,
                         capture_output=True, text=True, check=True).stdout
    return np.array([list(map(float, line.split())) for line in out.split("\n") if line.strip()])


def inlet(law, U, e, d, L, n=4000):
    """Gauge inlet pressure of a column blown at U (superficial at p0)."""
    K = gr.drag_law(law)
    m = U * RG

    def dpdx(p):
        rho = p / RT
        u = m / (rho * e)
        return K(e, u, d, MU, rho) * u / e

    p, h = pc.P0, L / n
    for _ in range(n):
        k1 = dpdx(p)
        k2 = dpdx(p + 0.5 * h * k1)
        k3 = dpdx(p + 0.5 * h * k2)
        k4 = dpdx(p + h * k3)
        p += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return p - pc.P0


def main():
    work = sys.argv[1]
    print(f"Air at {pc.T} K: rho_g {RG:.6f} kg/m^3, mu {MU:.6e} Pa s; d {D} m\n")

    print("A  each law against its closed form (gas fraction 0.35-0.999, Re 1e-3-1e3, u = 0)")
    es = list(np.round(np.arange(0.35, 0.96, 0.05), 10)) + [0.79999999, 0.8, 0.85, 0.9, 0.99, 0.999]
    res_ = list(10.0 ** np.arange(-3, 3.01, 0.5))
    pairs = [(e, 0.0) for e in es] + [(e, r * MU / (RG * D)) for e in es for r in res_]
    for law in LAWS:
        t = table(law, pairs)
        ref = np.array([gr.drag_law(law)(e, u, D, MU, RG) for e, u in pairs])
        err = np.max(np.abs(t[:, 2] / ref - 1))
        report(f"A {law}", err <= 1e-10, f"{len(pairs)} points, largest relative difference {err:.2e} (<= 1e-10)")

    print("\nB  limits")
    t_ewy = table("ergunWenYu", pairs)
    t_e, t_w = table("ergun", pairs), table("wenYu", pairs)
    lo = np.array([e < 0.8 for e, _ in pairs])
    err = max(np.max(np.abs(t_ewy[lo, 2] / t_e[lo, 2] - 1)), np.max(np.abs(t_ewy[~lo, 2] / t_w[~lo, 2] - 1)))
    report("B ergunWenYu = ergun below 0.8, wenYu above", err <= 1e-14, f"{err:.1e}")
    e1 = 1 - 1e-7
    stokes = 18 * MU * (1 - e1) / D ** 2
    ph = 1 - e1
    f_b = e1 * (10 * ph / e1 ** 2 + e1 ** 2 * (1 + 1.5 * math.sqrt(ph)))
    for law, want in (("wenYu", e1 ** -2.65), ("ergunWenYu", e1 ** -2.65), ("beetstra", f_b),
                      ("diFelice", 23.04 / 24 * e1 ** -2.7),
                      ("syamlalOBrien", 23.04 / 24 * e1 ** (1 - 4.14))):
        k = table(law, [(e1, 0.0)])[0, 2]
        report(f"B {law} one sphere at e -> 1", abs(k / stokes / want - 1) <= 1e-12,
               f"K / (18 mu (1-e)/d^2) = {k / stokes:.7f} (closed form {want:.7f}; e = 1 - 1e-7)")
    e = 0.4
    k = table("beetstra", [(e, 0.0)])[0, 2]
    carman = 180 * (1 - e) ** 2 * MU / (e ** 3 * D * D) * e ** 2
    second = 18 * MU * e * (1 - e) * e * e * (1 + 1.5 * math.sqrt(1 - e)) / D ** 2
    report("B beetstra packed bed at rest", abs((k - second) / carman - 1) <= 1e-12,
           f"K = {k:.6g}, Carman-Kozeny {carman:.6g} plus the e^2 (1 + 1.5 sqrt(phi)) term {second:.6g}")

    print("\nC  permeameter with each law (e 0.40, d 1 mm, 0.2 m, 40 cells, U 0.3 m/s)")
    L, n, U = 0.2, 40, 0.3
    for law in LAWS:
        path = f"{work}/C_{law}"
        pc.make_case(path, size=(L / n, L), cells=(1, n), e=0.4, d=D, end=0.04, write_interval=0.04,
                     p_patches={"top": "type fixedValue; value uniform 0;",
                                "bottom": "type fixedGradient; gradient uniform 0;"},
                     inflow={"bottom": U})
        g = f"{path}/constant/gasProperties"
        txt = open(g).read().replace("drag    ergunWenYu;", f"drag    {law};")
        pc.write(g, txt)
        pc.run(path)
        t = pc.times(path)[-1]
        p = pc.field(path, t, "p", n)
        pin = p[0] + pc.patch_gradient(path, t, "p", "bottom") * (L / n) / 2
        ref = inlet(law, U, 0.4, D, L)
        err = pin / ref - 1
        report(f"C {law}", abs(err) <= 1e-4, f"inlet {pin:.4f} Pa against {ref:.4f} Pa: {err:+.2e} (<= 1e-4)")

    print(f"\ntest_drag_laws: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)})")
    return 0 if all(res) else 1


if __name__ == "__main__":
    sys.exit(main())
