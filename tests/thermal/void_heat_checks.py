"""WP T.5 done-when: gas heat in the void (decision 0135).

A plane channel of void (gasProperties model gasFlow, no material), height
H 2 mm, length 60 mm, 200 x 16 cells, air at a mean speed of 0.53 m/s
(Peclet number u_m H/a about 50, Reynolds about 140, laminar). The flow
starts from the Poiseuille profile it keeps (0/Ug read, the pressure drop
12 mu u_m L/H^2 between the open ends). Air enters 10 K above the walls,
which are held at the starting temperature (0/Tg fixedValue); model
twoTemperature.

A  Thermally developed laminar flow between plates at a fixed wall
   temperature has Nu = h D_h/k = 7.5407 (D_h = 2H; Shah and London 1978;
   tools/thermal_reference.py graetz_plates_nu). The bulk (flow weighted)
   temperature difference to the walls decays as exp(-lambda x) with
   lambda = 2 h/(G c_g H), G the mass flux per area; Nu is taken from the
   slope of ln(theta_b) between x = 20 and 55 mm at the end, and must be
   within 3 % of 7.5407 (upwind transport adds about 1.1 %, axial
   conduction about 0.3 %). The profile is steady: theta_b at 0.2 and 0.3 s
   agree to 1e-3 of the inlet step.
B  Energy balance (thermal energy - start - heat in through the walls -
   heat carried in by the air) under 1e-10 of the thermal energy.

    python3 void_heat_checks.py WORKDIR
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
import flow_cases as fc  # noqa: E402
import pore_cases as pc  # noqa: E402
import thermal_reference as tr  # noqa: E402

report = cc.report
failed = cc.failed

H, LC = 2.0e-3, 60.0e-3
NX, NY = 200, 16
U_M = 0.53
K_G, C_G = 0.02569, 1006.0
MU = gr.air_viscosity(pc.T)
DP = 12 * MU * U_M * LC / H ** 2
T0, DT0 = pc.T, 10.0
END = 0.3


def make(path):
    pc.make_case(path, size=(LC, H), cells=(NX, NY), e=0.4, d=1e-3, end=END,
                 write_interval=0.1,
                 p_patches={"left": f"type fixedValue; value uniform {DP!r};",
                            "right": "type fixedValue; value uniform 0;"},
                 coupled=True)
    fc._void(path, "gasFlow", "pTolerance 1e-9;\n")
    # The solid (no material) only sets the time step: E 3300 gives about
    # 4e-5 s (gasFlow is not stable at 2e-4 s here)
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read()
    assert f"E {fc.E_SOLID!r};" in txt
    open(mp, "w").write(txt.replace(f"E {fc.E_SOLID!r};", "E 3300.0;"))
    xc = (np.arange(NX) + 0.5) * LC / NX
    yc = (np.arange(NY) + 0.5) * H / NY
    X, Y = np.meshgrid(xc, yc)                      # cell i = ix + NX iy
    u = (DP / (2 * MU * LC) * Y * (H - Y)).ravel()
    p = (DP * (1 - X / LC)).ravel()
    pc.write(f"{path}/0/Ug", pc.HEADER.format(cls="volVectorField", obj="Ug") + f"""
dimensions [0 1 -1 0 0 0 0];
internalField nonuniform List<vector> {u.size} ({' '.join(f'({float(v)!r} 0 0)' for v in u)});
boundaryField
{{
    left {{ type zeroGradient; }}
    right {{ type zeroGradient; }}
    bottom {{ type fixedValue; value uniform (0 0 0); }}
    top {{ type fixedValue; value uniform (0 0 0); }}
    frontAndBack {{ type empty; }}
}}
""")
    txt = open(f"{path}/0/p").read()
    txt = txt.replace("internalField uniform 0;", "internalField nonuniform List<scalar> "
                      f"{p.size} ({' '.join(repr(float(v)) for v in p)});")
    open(f"{path}/0/p", "w").write(txt)
    cc.write(f"{path}/constant/thermalProperties",
             cc.HEADER.format(cls="dictionary", obj="thermalProperties") + f"""
model        twoTemperature;
c            1000;
conductivity constant;
k            0;
kGas         {K_G!r};
cGas         {C_G!r};
""")
    walls = ("fixedValue", f"value uniform {T0!r};")
    cc.write_T(path, T0, {"left": ("fixedValue", f"value uniform {T0 + DT0!r};"),
                          "right": ("zeroGradient", ""), "bottom": walls, "top": walls})
    os.rename(f"{path}/0/T", f"{path}/0/Tg")
    txt = open(f"{path}/0/Tg").read().replace("object      T;", "object      Tg;")
    open(f"{path}/0/Tg", "w").write(txt)
    cc.write_T(path, T0, {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})


def bulk(path, t):
    Tg = cc.read_internal(f"{path}/{t:g}/Tg").reshape(NY, NX)
    U = pc.field(path, f"{t:g}", "Ug", NX * NY)[:, 0].reshape(NY, NX)
    return ((Tg - T0) * U).sum(axis=0) / U.sum(axis=0) / DT0, U


def check(work):
    print("[A] Developed laminar heat transfer between plates (Graetz)", flush=True)
    p = os.path.join(work, "graetz")
    make(p)
    t0 = time.perf_counter()
    pc.run(p)
    secs = time.perf_counter() - t0
    th, U = bulk(p, END)
    th4, _ = bulk(p, END - 0.1)
    x = (np.arange(NX) + 0.5) * LC / NX
    sel = (x > 0.020) & (x < 0.055)
    lam = -np.polyfit(x[sel], np.log(th[sel]), 1)[0]
    rho = (pc.P0 + 0.5 * DP) / (287.05 * pc.T)
    G = rho * U.mean(axis=0)[sel].mean()
    h = lam * G * C_G * H / 2
    nu = h * 2 * H / K_G
    want = tr.graetz_plates_nu()
    report("Nu against 7.5407", abs(nu / want - 1) < 0.03,
           f"Nu {nu:.4f} ({(nu / want - 1) * 100:+.2f} %, limit 3 %), mean speed "
           f"{U.mean():.4f} m/s, {secs:.0f} s")
    st = np.abs(th - th4).max()
    report("steady profile", st < 1e-3, f"theta_b changed by {st:.1e} from {END - 0.1:g} to {END:g} s")

    print("[B] Energy balance", flush=True)
    log = ca.energy_log(p)
    rel = max(abs(b) for b in log["bal"]) / log["E"][0]
    report("energy balance", rel < 1e-10,
           f"largest |balance| {rel:.1e} of the thermal energy; heat in {log['heat'][-1]:.4e} J, "
           f"carried in {log['carried'][-1]:.4e} J")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check(work)
    print(f"void heat checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_void_heat: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
