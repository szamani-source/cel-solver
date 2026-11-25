"""WP F.5 done-when: air-impeded hopper discharge (decisions 0125, 0130).
Cases in tests/gas/hopper_gas_cases.py: plane-flow hopper, 60 mm slot,
544 body-fitted cells, outlet opens at 0.15 s, run to 0.5 s; drag
ergunWenYu, two-way coupling, freeSurface solid.

A  Grain size sweep, v1's T0 rule set (the hopper set, mu(I) with the
   case's d): gas on and gas off for d = 4, 1, 0.5, 0.2 and 0.1 mm.
   Steady rate W (least squares of the discharged mass, 0.3 to 0.5 s).
   Checks: dry rate within 20 % of Beverloo (slot, C 1.03, k 1.5); dry
   rate the same for every d <= 1 mm to 2 % (there d acts only through
   the gas; at 4 mm mu(I) and Beverloo's k d lower it, reported);
   d >= 1 mm with gas within 5 % of dry; the rate falls with d; 0.1 mm at
   least 20 % slower with gas (published: standard correlations
   over-predict below about 0.5 mm, Hsiau et al. 2010; Crewdson et al.
   1977 could not be read, see docs/spec/gas-references.md); gas mass
   less the inflow constant to 1e-12; the outlet holds p0.
   Reported: Muite's N = U_t / U_out for each d, and Crewdson's model
   W/W_dry = sqrt(1 + G/(rho_b g)) with G the gas pressure gradient over
   the outlet row of the run.
P  2 and 4 processors against serial: plain Mohr-Coulomb (no T0 rules,
   whose switches part the solid serial vs parallel, decision 0080),
   d 0.1 mm, gas on, to 0.25 s (material leaving); cuts through the
   outlet, the bed and the free surface (2 x 1 and 2 x 2). p within 1e-10
   of the largest |p| (solve to 1e-12 Pa), U within 1e-10, gas mass 1e-12.
C  Cost per step, gas on against gas off (report; F.9 sets the 1.5x limit).

    python3 air_impeded_checks.py WORKDIR
"""

import math
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import hopper_gas_cases as hc  # noqa: E402
import pore_cases as pc  # noqa: E402

SIZES = (4e-3, 1e-3, 5e-4, 2e-4, 1e-4)
END, WRITE = 0.5, 0.025
WINDOW = (0.3, 0.5)
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def drift(path):
    m, ni, _, _ = pc.mass_log(path)
    return float(np.abs((m - ni) - (m[0] - ni[0])).max() / m[0])


def steps(path):
    return len(re.findall(r"^Time = ", open(f"{path}/log.celFoam").read(), re.M))


def discharge(path):
    """Write times and discharged mass per metre of slot [kg/m]."""
    txt = open(f"{path}/log.celFoam").read()
    m = np.array([float(x) for x in re.findall(r"Discharged mass = (\S+)", txt)])
    t = np.array([0.0] + [float(x) for x in pc.times(path)])
    n = min(len(t), len(m))
    return t[:n], m[:n] / hc.depth()


def rate(path):
    t, m = discharge(path)
    s = (t >= WINDOW[0] - 1e-9) & (t <= WINDOW[1] + 1e-9)
    return float(np.polyfit(t[s], m[s], 1)[0])


def outlet_gradient(path, t):
    """Mean gas pressure gradient over the outlet row [Pa/m] (negative:
    the bed is below p0, the air pushes the grains up)."""
    _, h0 = hc.outlet_row()
    p = pc.field(path, t, "p", hc.ncells())
    return float(p[: hc.NCOL].mean()) / (h0 / 2)


def muite_n(d, w):
    rho_g, mu = gr.air_density(), gr.air_viscosity()
    u_out = w / (hc.RHO * 2 * hc.RO)
    return gr.terminal_velocity(d, pc.RHOP, rho_g, mu) / u_out


def main(work):
    nc = hc.ncells()
    print(f"Hopper: {nc} cells, slot {2 * hc.RO * 1e3:g} mm, walls {hc.BETA:g} deg, "
          f"opens at {hc.OPEN_TIME} s, run to {END} s")

    # A
    print("A: grain size sweep, gas on and off (T0 rules)")
    jobs = []
    for d in SIZES:
        for gas in (False, True):
            p = f"{work}/A_{d * 1e3:g}mm_{'gas' if gas else 'dry'}"
            hc.make_hopper(p, d=d, end=END, write_interval=WRITE, gas=gas)
            jobs.append(p)
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda p: pc.run(p, set_fields=True), jobs))
    W = {}
    for d in SIZES:
        W[d] = (rate(f"{work}/A_{d * 1e3:g}mm_dry"), rate(f"{work}/A_{d * 1e3:g}mm_gas"))
    bev = hc.beverloo_slot(1e-3)
    wd = W[1e-3][0]
    report("dry rate against Beverloo (slot)", abs(wd - bev) <= 0.2 * bev,
           f"{wd:.2f} kg/(m s) against {bev:.2f} ({(wd - bev) / bev:+.1%})")
    # d enters the dry run only through mu(I) and the size of the grains at
    # the outlet edge (Beverloo's k d): small for d <= 1 mm
    spread = max(abs(W[d][0] - wd) for d in SIZES if d <= 1e-3) / wd
    report("dry rate independent of d (d <= 1 mm)", spread <= 0.02,
           f"largest change {spread:.2%}; 4 mm {(W[4e-3][0] - wd) / wd:+.2%} "
           f"(Beverloo {(hc.beverloo_slot(4e-3) - bev) / bev:+.2%})")
    print("  d [mm]   dry      gas      change   Muite N   Crewdson sqrt(1+G/rho g)")
    for d in SIZES:
        g = f"{work}/A_{d * 1e3:g}mm_gas"
        G = np.mean([outlet_gradient(g, t) for t in pc.times(g) if float(t) >= WINDOW[0]])
        cw = math.sqrt(max(0.0, 1 + G / (hc.RHO * hc.G_ACC)))
        print(f"  {d * 1e3:6.2f}  {W[d][0]:7.2f}  {W[d][1]:7.2f}  "
              f"{(W[d][1] - W[d][0]) / W[d][0]:+7.1%}  {muite_n(d, W[d][1]):7.2f}   {cw:.3f}")
    for d in SIZES:
        if d >= 1e-3:
            ch = (W[d][1] - W[d][0]) / W[d][0]
            report(f"{d * 1e3:g} mm with gas within 5 % of dry", abs(ch) <= 0.05, f"{ch:+.2%}")
    ratios = [W[d][1] / W[d][0] for d in SIZES]
    report("rate falls with grain size", all(a > b for a, b in zip(ratios, ratios[1:])),
           " > ".join(f"{r:.3f}" for r in ratios))
    r01 = ratios[-1]
    report("0.1 mm at least 20 % slower with gas", r01 <= 0.8, f"W/W_dry {r01:.3f}")
    dm = max(drift(f"{work}/A_{d * 1e3:g}mm_gas") for d in SIZES)
    report("gas mass balance (A)", dm <= 1e-12, f"{dm:.1e}")
    g = f"{work}/A_0.1mm_gas"
    pout = max(float(np.abs(pc.patch_value(g, t, "p", "outlet")).max()) for t in pc.times(g))
    report("outlet holds p0", pout == 0.0, f"largest |p| on the outlet {pout:.1e} Pa")

    # P
    print("P: 2 and 4 processors against serial (plain Mohr-Coulomb, 0.1 mm, gas on)")
    runs = {}
    for n, dc in ((1, (1, 1, 1)), (2, (2, 1, 1)), (4, (2, 2, 1))):
        p = f"{work}/P{n}"
        hc.make_hopper(p, d=1e-4, end=0.25, write_interval=0.05, decomp=dc, rules=False,
                       extra_gas="pTolerance 1e-12;\n")
        runs[n] = pc.run(p, n, set_fields=True)
    pmax = max(float(np.abs(pc.field(f"{work}/P1", t, "p", nc)).max()) for t in pc.times(f"{work}/P1"))
    for n in (2, 4):
        dp = du = 0.0
        for t in pc.times(f"{work}/P1"):
            a, b = f"{work}/P1", f"{work}/P{n}"
            dp = max(dp, float(np.abs(pc.field(a, t, "p", nc) - pc.field(b, t, "p", nc)).max()) / pmax)
            ua, ub = pc.field(a, t, "U", nc), pc.field(b, t, "U", nc)
            du = max(du, float(np.abs(ua - ub).max() / max(np.abs(ua).max(), 1e-300)))
        report(f"{n} processors (cuts through the outlet, bed and surface)",
               dp <= 1e-10 and du <= 1e-10 and drift(f"{work}/P{n}") <= 1e-12,
               f"p within {dp:.1e} of the largest |p|, U within {du:.1e}, "
               f"gas mass {drift(f'{work}/P{n}'):.1e}")
    print(f"  wall times: 1 proc {runs[1]:.1f} s, 2 procs {runs[2]:.1f} s, 4 procs {runs[4]:.1f} s "
          f"({nc} cells: too small to scale)")

    # C
    print("C: cost per step, gas on against gas off (report, run one at a time)")
    for rules in (True, False):
        secs = {}
        for gas in (False, True):
            p = f"{work}/C_{int(rules)}{int(gas)}"
            hc.make_hopper(p, d=1e-4, end=END, write_interval=END, gas=gas, rules=rules)
            secs[gas] = pc.run(p, set_fields=True) / steps(p)
        _, _, sv, _ = pc.mass_log(f"{work}/C_{int(rules)}1")
        print(f"  {'T0 rules' if rules else 'plain Mohr-Coulomb'}, 0.1 mm: "
              f"{secs[True] * 1e3:.3f} ms against {secs[False] * 1e3:.3f} ms per step: "
              f"{secs[True] / secs[False]:.2f}x, {sv.mean():.2f} solves per step")

    print(f"{sum(res)}/{len(res)} checks passed")
    return all(res)


if __name__ == "__main__":
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    for d in os.listdir(work):
        shutil.rmtree(os.path.join(work, d), ignore_errors=True)
    sys.exit(0 if main(work) else 1)
