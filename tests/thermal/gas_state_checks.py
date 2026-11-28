"""WP T.7 done-when: gas properties from the gas temperature, and plastic
heating (decision 0136).

thermalProperties variableGas yes: the pore gas is an ideal gas at its own
temperature, rho_g = (p0 + p)/(R Tg), mu(Tg) and kGas(Tg) by Sutherland's
laws; after each heat step each cell's gas pressure follows its temperature
at constant gas mass. plasticHeating beta: beta sigma : d(epsP) per unit
material volume heats the material in each Lagrange step.

A  Sealed cell of gas (one cell of a closed bed, rigid, coupled no),
   twoTemperature, grains at 120 C and gas at 20 C: the gas heats up by Gunn
   exchange at constant volume, so (p0 + p)/Tg stays p0/Tg0 to 1e-10 at
   every write (ideal gas law); energy balance to 1e-12. (One cell: in a
   uniform box the round-off of the gas flux, through Re^0.2 in Gunn's
   Nusselt number, starts small flows between cells.)
B  Sealed box heated through one wall (equilibrium, conduction, h 1000 to
   120 C): the gas flows from the hot side to the cold; the gas mass sum V
   phi_g (p0 + p)/(R T), from the written fields, stays at its start to
   1e-12; energy balance 1e-12; 2 processors match serial to 1e-10.
C  Plastic heating, closed form: a closed box of perfectly plastic von
   Mises material starting with a uniform compressive stress outside the
   yield surface (q0 = 2 sigma_y). In the first step the return puts q on
   sigma_y at rest: epsPEq = (q0 - sigma_y)/(3 G) and, the plastic work of
   J2 plasticity being sigma_y d(epsPEq), T - T0 = beta sigma_y epsPEq/(rho c),
   both to 1e-10. Over 0.02 s the box moves: energy balance with the
   plastic heat 1e-12. With plasticHeating 0 T does not change.
D  Plastic heating in a moving case: tutorials/solidFluidThermal/hotSiloDischarge (E 0.1 MPa,
   to 0.8 s) with plasticHeating 0.9: energy balance (energy - start - heat
   in - carried in - plastic heat) to 1e-12, the plastic heat positive.
E  T.4 and T.5 with variableGas yes and a 1 K step (the properties then
   vary by 0.3 %): Schumann within 1 % and the developed Nusselt number
   within 3 %, their limits.

    python3 gas_state_checks.py WORKDIR
"""

import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "gas"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import conduction_checks as cc  # noqa: E402
import carried_checks as ca  # noqa: E402
import pore_cases as pc  # noqa: E402

report = cc.report
failed = cc.failed

R_AIR = 287.05
E_GAS, D = 0.4, 1.0e-3
T_GAS, T_HOT = pc.T, pc.T + 100.0


def field(path, n=None):
    """Internal field of a scalar file (a uniform field repeated n times)."""
    v = cc.read_internal(path)
    if v is None:
        s = open(path).read()
        val = float(s.split("internalField")[1].split()[1].rstrip(";"))
        v = np.full(n or 1, val)
    return v


def write_tp(path, body):
    cc.write(f"{path}/constant/thermalProperties",
             cc.HEADER.format(cls="dictionary", obj="thermalProperties") + body)


def rename_T(path, name):
    os.rename(f"{path}/0/T", f"{path}/0/{name}")
    txt = open(f"{path}/0/{name}").read().replace("object      T;", f"object      {name};")
    open(f"{path}/0/{name}", "w").write(txt)


# ------------------------------------------------------------------ A

def check_sealed_exchange(work):
    print("[A] Sealed cell: gas heated by the grains at constant volume", flush=True)
    p = os.path.join(work, "sealed")
    nx = ny = 1
    pc.make_case(p, size=(0.04, 0.04), cells=(nx, ny), e=E_GAS, d=D, end=0.01,
                 write_interval=0.002, p_patches={}, coupled=False)
    write_tp(p, f"""
model        twoTemperature;
c            800;
conductivity constant;
k            0;
variableGas  yes;
""")
    zg = {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")}
    cc.write_T(p, T_GAS, zg)
    rename_T(p, "Tg")
    cc.write_T(p, T_HOT, zg)
    pc.run(p)
    n = nx * ny
    worst, rise = 0.0, 0.0
    for t in pc.times(p):
        pg = field(f"{p}/{t}/p", n)
        Tg = field(f"{p}/{t}/Tg", n)
        worst = max(worst, np.abs((pc.P0 + pg) / Tg / (pc.P0 / T_GAS) - 1).max())
        rise = max(rise, Tg.max() - T_GAS)
    report("ideal gas at constant volume: (p0 + p)/Tg = p0/Tg0", worst < 1e-10 and rise > 50,
           f"largest relative difference {worst:.1e} (limit 1e-10); the gas rose by {rise:.1f} K")
    log = ca.energy_log(p)
    rel = max(abs(b) for b in log["bal"]) / log["E"][0]
    report("energy balance", rel < 1e-12, f"largest |balance| {rel:.1e} of the thermal energy")


# ------------------------------------------------------------------ B

NXB, NYB, LXB, LYB = 20, 10, 0.1, 0.05


def heated_box(path, decomp=(1, 1, 1)):
    pc.make_case(path, size=(LXB, LYB), cells=(NXB, NYB), e=E_GAS, d=D, end=0.3,
                 write_interval=0.1, p_patches={}, decomp=decomp, coupled=False)
    write_tp(path, f"""
model        equilibrium;
c            5;
conductivity constant;
k            2.5;
variableGas  yes;
walls {{ left {{ h 1000; Tinf {T_HOT!r}; }} }}
""")
    cc.write_T(path, T_GAS, {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})
    pc.run(path, decomp[0] * decomp[1] * decomp[2])


def gas_mass(path, t):
    n = NXB * NYB
    pg = field(f"{path}/{t}/p", n)
    T = field(f"{path}/{t}/T", n)
    rho = field(f"{path}/{t}/rho", n) if os.path.exists(f"{path}/{t}/rho") else \
        np.full(n, (1 - E_GAS) * pc.RHOP)
    dx, dy = LXB / NXB, LYB / NYB
    V = dx * dy * min(dx, dy)
    phiG = 1 - rho / pc.RHOP
    return (V * phiG * pc.P0 / (R_AIR * T)).sum() + (V * phiG * pg / (R_AIR * T)).sum(), pg, T


def check_heated_box(work):
    print("[B] Sealed box heated through one wall", flush=True)
    p = os.path.join(work, "heated")
    t0 = time.perf_counter()
    heated_box(p)
    secs = time.perf_counter() - t0
    ts = pc.times(p)
    m0, _, _ = gas_mass(p, "0")
    worst = 0.0
    for t in ts:
        m, pg, T = gas_mass(p, t)
        worst = max(worst, abs(m / m0 - 1))
    report("gas mass kept", worst < 1e-12,
           f"largest change {worst:.1e} (limit 1e-12); at {ts[-1]} s T {T.min() - 273.15:.2f} to "
           f"{T.max() - 273.15:.2f} C, p {pg.min():.1f} to {pg.max():.1f} Pa, {secs:.0f} s")
    log = ca.energy_log(p)
    rel = max(abs(b) for b in log["bal"]) / log["E"][0]
    report("energy balance", rel < 1e-12, f"largest |balance| {rel:.1e} of the thermal energy")

    p2 = os.path.join(work, "heated_par")
    heated_box(p2, decomp=(2, 1, 1))
    t = ts[-1]
    dT = np.abs(field(f"{p2}/{t}/T") - field(f"{p}/{t}/T")).max() / (T_HOT - T_GAS)
    dp = np.abs(field(f"{p2}/{t}/p") - field(f"{p}/{t}/p")).max() / max(np.abs(pg).max(), 1.0)
    report("2 processors match serial", max(dT, dp) < 1e-10,
           f"largest difference T {dT:.1e} of the step, p {dp:.1e} of the largest p")


# ------------------------------------------------------------------ C

RHO_S, E_S, NU_S = 1800.0, 1.0e6, 0.3
C_S = 1.0       # small, so the rise (0.04 K) is far above the round-off of T
SIG_Y, BETA = 10.0e3, 0.9


def stressed_box(path, material, sigma, beta, end=0.02):
    cc.common(path, end=end, write_interval=end)
    cc.write(f"{path}/constant/materialProperties",
             cc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO_S!r};
{material}
walls {{ left {{ friction 0; }} right {{ friction 0; }} bottom {{ friction 0; }} top {{ friction 0; }} }}
""")
    write_tp(path, f"""
model          equilibrium;
c              {C_S!r};
conductivity   constant;
k              0;
plasticHeating {beta!r};
""")
    cc.blockmesh(path, 0.1, 0.1, 4, 4)
    cc.write_T(path, T_GAS, {n: ("zeroGradient", "") for n in ("left", "right", "bottom", "top")})
    sx, sy, sz = sigma
    cc.write(f"{path}/0/sigma", cc.HEADER.format(cls="volSymmTensorField", obj="sigma") + f"""
dimensions [1 -1 -2 0 0 0 0];
internalField uniform ({sx!r} 0 0 {sy!r} 0 {sz!r});
boundaryField {{ ".*" {{ type zeroGradient; }} frontAndBack {{ type empty; }} }}
""")
    cc.run(path)


def check_plastic_closed_form(work):
    print("[C] Plastic heating against the closed form (von Mises)", flush=True)
    G = E_S / (2 * (1 + NU_S))
    q0 = 2 * SIG_Y
    a = float(q0 / np.sqrt(3))                # s = diag(-a, a, 0): q = sqrt(3) a
    vm_sigma = (-3 * a, -a, -2 * a)           # all compressive
    vm = f"type vonMises;\nE {E_S!r};\nnu {NU_S!r};\nsigma_y {SIG_Y!r};\nH 0;"
    n = 16

    # One step (shorter than the stable one): the return of the starting
    # stress, at rest
    p = os.path.join(work, "vm_step")
    stressed_box(p, vm, vm_sigma, BETA, end=3e-4)
    T = field(f"{p}/0.0003/T", n)
    eq = field(f"{p}/0.0003/epsPEq", n)
    eq_want = (q0 - SIG_Y) / (3 * G)
    dT_want = BETA * SIG_Y * eq_want / (RHO_S * C_S)
    yielded = eq > 0
    e1 = np.abs(eq[yielded] / eq_want - 1).max() if yielded.any() else 1.0
    e2 = np.abs((T - T_GAS) - dT_want * eq / eq_want).max() / dT_want
    report("one step: epsPEq = (q0 - sigma_y)/(3 G), T - T0 = beta sigma_y epsPEq/(rho c)",
           yielded.any() and max(e1, e2) < 1e-10,
           f"{yielded.sum()} of {n} cells returned in the first step (the others at the "
           f"next, with no strain increment the update is skipped); rise {dT_want:.6e} K, "
           f"largest error {e1:.1e} (epsPEq), {e2:.1e} (T)")

    # Over 0.02 s: the box moves, the energy balance holds
    p = os.path.join(work, "vm")
    stressed_box(p, vm, vm_sigma, BETA)
    T = field(f"{p}/0.02/T", n)
    eq = field(f"{p}/0.02/epsPEq", n)
    approx = np.abs((T - T_GAS) / (BETA * SIG_Y * eq / (RHO_S * C_S)) - 1).max()
    log = ca.energy_log(p)
    ph = [float(v) for v in re.findall(r"^Plastic heat\s*=\s*(\S+)",
                                         open(f"{p}/log.celFoam").read(), re.M)]
    rel = max(abs(b) for b in log["bal"]) / log["E"][0]
    report("moving: energy balance with the plastic heat", rel < 1e-12 and ph[-1] > 0,
           f"largest |balance| {rel:.1e} of the thermal energy, plastic heat {ph[-1]:.4e} J; "
           f"T - T0 within {approx:.1e} of beta sigma_y epsPEq/(rho c) at the starting "
           f"density (the density changes as it moves)")

    p = os.path.join(work, "vm_off")
    stressed_box(p, vm, vm_sigma, 0.0)
    T = field(f"{p}/0.02/T", n)
    report("plasticHeating 0: T unchanged", np.abs(T - T_GAS).max() == 0,
           f"largest change {np.abs(T - T_GAS).max():.1e} K")


# ------------------------------------------------------------------ D

def check_silo(work):
    print("[D] Plastic heating in the hot silo tutorial", flush=True)
    p = os.path.join(work, "hotSiloDischarge")
    shutil.copytree(os.path.join(ROOT, "tutorials", "solidFluidThermal", "hotSiloDischarge"), p)
    for k, v in (("endTime", "0.8"), ("writeInterval", "0.1")):
        subprocess.run(["foamDictionary", "-entry", k, "-set", v, "system/controlDict"],
                       cwd=p, check=True, stdout=subprocess.DEVNULL)
    tp = f"{p}/constant/thermalProperties"
    txt = open(tp).read()
    open(tp, "w").write(txt.replace("walls\n{", "plasticHeating 0.9;\n\nwalls\n{", 1))
    t0 = time.perf_counter()
    r = subprocess.run(["./Allrun", "0.1"], cwd=p, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    secs = time.perf_counter() - t0
    if r.returncode != 0:
        report("hot silo runs", False, "Allrun failed")
        return
    log = open(f"{p}/log.celFoam").read()
    e0 = float(re.search(r"^Thermal energy at start\s*=\s*(\S+)", log, re.M).group(1))
    bal = [abs(float(v)) for v in re.findall(r"^Energy balance\s*=\s*(\S+)", log, re.M)]
    ph = [float(v) for v in re.findall(r"^Plastic heat\s*=\s*(\S+)", log, re.M)]
    rel = max(bal) / e0
    report("energy balance with the plastic heat", rel < 1e-12 and len(ph) > 0 and ph[-1] > 0,
           f"largest |balance| {rel:.1e} of the thermal energy; plastic heat {ph[-1]:.4e} J "
           f"({ph[-1] / e0:.1e} of the thermal energy), {secs:.0f} s")


# ------------------------------------------------------------------ E

def with_variable_gas(fn):
    """Run fn with every thermalProperties written by cc.write getting
    variableGas yes."""
    orig = cc.write

    def write(path, text):
        if path.endswith("constant/thermalProperties"):
            text += "variableGas  yes;\n"
        orig(path, text)
    cc.write = write
    try:
        return fn()
    finally:
        cc.write = orig


def check_t4_t5(work):
    print("[E] T.4 and T.5 with variableGas and a 1 K step", flush=True)
    import gas_heat_checks as gh
    import void_heat_checks as vh
    import thermal_reference as tr
    gh.DT0 = 1.0
    p = os.path.join(work, "schumann_var")
    with_variable_gas(lambda: gh.column(p, "twoTemperature", gh.T0 + 1.0))
    wg, ws = gh.schumann_errors(p)
    report("T.4 Schumann with variableGas", max(wg, ws) < 0.01,
           f"largest error gas {wg * 100:.2f} %, grains {ws * 100:.2f} % of the step (limit 1 %)")

    vh.DT0 = 1.0
    p = os.path.join(work, "graetz_var")
    with_variable_gas(lambda: vh.make(p))
    pc.run(p)
    nu, _, _ = vh.nusselt(p)
    want = tr.graetz_plates_nu()
    report("T.5 Nusselt with variableGas", abs(nu / want - 1) < 0.03,
           f"Nu {nu:.4f} against {want:.4f} ({(nu / want - 1) * 100:+.2f} %, limit 3 %)")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check_sealed_exchange(work)
    check_heated_box(work)
    check_plastic_closed_form(work)
    check_silo(work)
    check_t4_t5(work)
    print(f"gas state checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_thermal_gas_state: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
