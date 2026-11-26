"""WP T.4 done-when: gas temperature and grain-gas exchange (decision 0135).

A rigid packed column (poreGas, coupled no, no gravity) 0.6 m long, 300
cells of 2 mm, grains d 20 mm, gas fraction 0.4, air blown in through the
left face at U 1 m/s (superficial, at p0) and out at the right (p0). The
grains' specific heat is small (c 5 J/(kg K)) so the thermal front crosses
the column in seconds; nothing in the closed forms depends on the values.

A  twoTemperature: air enters 100 K hotter than the bed (0/Tg fixedValue at
   the inlet), no conduction (k 0; the gas conduction is negligible at a
   Peclet number of 1e4). Against Schumann (1929), tools/thermal_reference.py
   schumann_gas/solid, with h from Gunn (1978) at the column's Reynolds and
   Prandtl numbers: xi = h a_s x/(G c_g), tau = h a_s (t - x/u_i)/C_s,
   a_s = 6 (1 - e)/d. theta_g and theta_s within 1 % of the inlet step at
   t = 1 and 2 s (tau about 2.6 and 5.3); 1 % catches an error of about
   7 % in h (20 % in h moves theta by 3 %). On 2 processors (cut at the
   middle) T and Tg match serial to 1e-8 of the step.
B  equilibrium: T at the inlet 100 K hotter (0/T fixedValue). The front
   (theta 0.5) moves at e rho_g c_g u_i/(e rho_g c_g + (1 - e) rho_P c)
   (equilibrium_front_speed) within 1 %, measured between t = 1 and 2 s.
C  Both models with the inlet at the bed's temperature: the temperature
   stays uniform to 1e-12 under the flow. Energy balance (thermal energy -
   start - heat in - heat carried in) of A and B: changes by less than 1e-9
   of the thermal energy after the gas flow has settled (first write).

    python3 gas_heat_checks.py WORKDIR
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "gas"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import conduction_checks as cc  # noqa: E402
import carried_checks as ca  # noqa: E402
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402
import thermal_reference as tr  # noqa: E402

report = cc.report
failed = cc.failed

L, N = 0.6, 300
DX = L / N
E_GAS, D, U = 0.4, 0.02, 1.0
C_S, C_G, K_G = 5.0, 1006.0, 0.02569
T0, DT0 = pc.T, 100.0
END, WRITE = 2.0, 0.5

MU = gr.air_viscosity(pc.T)
RHO_G = pc.P0 / (287.05 * pc.T)
G = RHO_G * U                              # gas mass flux [kg/(m^2 s)]
U_I = U / E_GAS
RHO_BULK = (1 - E_GAS) * pc.RHOP
PR = MU * C_G / K_G
RE = G * D / MU
H = tr.gunn_nu(E_GAS, RE, PR) * K_G / D
A_S = 6 * (1 - E_GAS) / D
HA = H * A_S
CS = RHO_BULK * C_S


def column(path, model, T_in, decomp=(1, 1, 1), end=END, write=WRITE):
    pc.make_case(path, size=(L, DX), cells=(N, 1), e=E_GAS, d=D, end=end,
                 write_interval=write,
                 p_patches={"right": "type fixedValue; value uniform 0;",
                            "left": "type fixedGradient; gradient uniform 0;"},
                 inflow={"left": U}, decomp=decomp, E=1.0e4, coupled=False)
    cc.write(f"{path}/constant/thermalProperties",
             cc.HEADER.format(cls="dictionary", obj="thermalProperties") + f"""
model        {model};
c            {C_S!r};
conductivity constant;
k            0;
kGas         {K_G!r};
cGas         {C_G!r};
""")
    zg = {n: ("zeroGradient", "") for n in ("right", "bottom", "top")}
    two = model == "twoTemperature"
    inlet = ("fixedValue", f"value uniform {T_in!r};")
    cc.write_T(path, T0, dict(zg, left=("zeroGradient", "") if two else inlet))
    if two:
        cc.write_T(path, T0, dict(zg, left=inlet))
        os.rename(f"{path}/0/T", f"{path}/0/Tg")
        txt = open(f"{path}/0/Tg").read().replace("object      T;", "object      Tg;")
        open(f"{path}/0/Tg", "w").write(txt)
        cc.write_T(path, T0, dict(zg, left=("zeroGradient", "")))
    n = decomp[0] * decomp[1] * decomp[2]
    pc.run(path, n)


def field(path):
    """Internal field as an array (a uniform field gives its one value)."""
    v = cc.read_internal(path)
    if v is None:
        s = open(path).read()
        v = np.array([float(s.split("internalField")[1].split()[1].rstrip(";"))])
    return v


def balance_drift(path):
    log = ca.energy_log(path)
    bal = log["bal"]
    return max(abs(b - bal[0]) for b in bal) / log["E"][0], log


def check_schumann(work):
    print("[A] Hot air into a cold bed against Schumann", flush=True)
    print(f"    Re {RE:.1f}, Pr {PR:.4f}, Gunn h {H:.2f} W/(m^2 K), h a_s {HA:.1f} W/(m^3 K), "
          f"xi per cell {HA * DX / (G * C_G):.4f}", flush=True)
    p = os.path.join(work, "schumann")
    t0 = time.perf_counter()
    column(p, "twoTemperature", T0 + DT0)
    secs = time.perf_counter() - t0
    x = (np.arange(N) + 0.5) * DX
    sample = np.arange(4, N, 10)
    worst_g = worst_s = 0.0
    for t in (1.0, 2.0):
        Tg = cc.read_internal(f"{p}/{t:g}/Tg")
        Ts = cc.read_internal(f"{p}/{t:g}/T")
        for i in sample:
            xi = HA * x[i] / (G * C_G)
            tau = HA * (t - x[i] / U_I) / CS
            if tau <= 0:
                eg, es = 0.0, 0.0
            else:
                eg = tr.schumann_gas(xi, tau, 1000)
                es = tr.schumann_solid(xi, tau, 1000)
            worst_g = max(worst_g, abs((Tg[i] - T0) / DT0 - eg))
            worst_s = max(worst_s, abs((Ts[i] - T0) / DT0 - es))
        print(f"    t {t:g} s: tau at the inlet {HA * t / CS:.2f}", flush=True)
    report("gas temperature against Schumann", worst_g < 0.01,
           f"largest error {worst_g * 100:.2f} % of the step (limit 1 %), {secs:.0f} s")
    report("grain temperature against Schumann", worst_s < 0.01,
           f"largest error {worst_s * 100:.2f} % of the step (limit 1 %)")
    drift, log = balance_drift(p)
    report("twoTemperature energy balance", drift < 1e-9,
           f"balance changed by {drift:.1e} of the thermal energy after the first write; "
           f"carried in {log['carried'][-1]:.4e} J")

    p2 = os.path.join(work, "schumann_par")
    column(p2, "twoTemperature", T0 + DT0, decomp=(2, 1, 1), end=0.5, write=0.5)
    d = max(np.abs(cc.read_internal(f"{p2}/0.5/{f}") - cc.read_internal(f"{p}/0.5/{f}")).max()
            for f in ("T", "Tg")) / DT0
    report("2 processors match serial", d < 1e-8, f"largest difference {d:.2e} of the step")


def front(T):
    th = (T - T0) / DT0
    i = int(np.argmax(th < 0.5))
    x = (np.arange(N) + 0.5) * DX
    return x[i - 1] + (th[i - 1] - 0.5) / (th[i - 1] - th[i]) * DX


def check_equilibrium(work):
    print("[B] Equilibrium front speed", flush=True)
    p = os.path.join(work, "equilibrium")
    column(p, "equilibrium", T0 + DT0)
    x1 = front(cc.read_internal(f"{p}/1/T"))
    x2 = front(cc.read_internal(f"{p}/2/T"))
    want = tr.equilibrium_front_speed(E_GAS, U_I, RHO_G, C_G, pc.RHOP, C_S)
    got = x2 - x1
    report("front speed", abs(got / want - 1) < 0.01,
           f"{got:.5f} m/s against {want:.5f} m/s ({(got / want - 1) * 100:+.2f} %)")
    drift, _ = balance_drift(p)
    report("equilibrium energy balance", drift < 1e-9,
           f"balance changed by {drift:.1e} of the thermal energy after the first write")


def check_uniform(work):
    print("[C] Uniform temperature under flow", flush=True)
    for model in ("twoTemperature", "equilibrium"):
        p = os.path.join(work, "uniform_" + model)
        column(p, model, T0, end=0.2, write=0.2)
        fields = ("T", "Tg") if model == "twoTemperature" else ("T",)
        d = max(np.abs(field(f"{p}/0.2/{f}") / T0 - 1).max() for f in fields)
        report(f"{model} stays uniform", d < 1e-12, f"largest relative change {d:.1e}")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check_uniform(work)
    check_schumann(work)
    check_equilibrium(work)
    print(f"gas heat checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_gas_heat: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
