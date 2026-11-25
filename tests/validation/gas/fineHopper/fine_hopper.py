"""Gas validation pack, test 4 (WP F.5, decisions 0125 item 9, 0130):
fine-powder hopper discharge with the pore gas on and off on the same mesh,
so the slowdown from the air shows.

Plane-flow hopper (tests/gas/hopper_gas_cases.py): slot 60 mm, walls 30 deg
from vertical, bin 0.2 m wide, 544 body-fitted cells (outlet cells 3.75 x
3.1 mm); fine quartz powder d 0.1 mm, grain density 2650 kg/m^3, gas
fraction 0.40 (bulk 1590 kg/m^3); Mohr-Coulomb phi 30, c 0, psi 1, E 1 MPa,
plain (without v1's T0 flow rules, whose switches part the solid between
processor counts, decision 0080; with them the dry rate is 74 instead of
33 kg/(m s), see tests/test_air_impeded.sh); fullBlend void; walls 20 deg;
air at 293.15 K, 101325 Pa; drag ergunWenYu. Outlet opens at 0.15 s; run
to 0.5 s.

Reference: Crewdson, Ormond and Nedderman (1977) measured that fine powders
discharge slower than coarse ones because the air must flow in to fill the
dilating bed; their model adds the air pressure gradient at the outlet to
gravity. The paper itself could not be read (docs/spec/gas-references.md),
so the check is the published direction and size of the effect: standard
correlations over-predict below about 0.5 mm (Hsiau et al. 2010), and at
0.1 mm the air slows the discharge by at least 20 %. Also reported:
Crewdson's sqrt(1 + G/(rho_b g)) with G the run's gradient over the outlet
row, and Muite's N = U_t/U_out (the gas matters below about 10).

Checks: gas on at least 20 % slower than gas off; 2 and 4 processors
(2 x 1 and 2 x 2: cuts through the outlet, the bed and the free surface)
match serial at every write to 1e-10 (U relative, p of the largest |p|,
with the pressure solve converged to 1e-12 Pa); gas mass less the boundary
flow constant to 1e-12. Writes times.csv and results.md.

    python3 fine_hopper.py [WORKDIR]      # default ./run
"""

import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tests", "gas"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import gas_reference as gr  # noqa: E402
import hopper_gas_cases as hc  # noqa: E402
import pore_cases as pc  # noqa: E402

D = 1e-4
END, WRITE = 0.5, 0.025
WINDOW = (0.3, 0.5)
DECOMP = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1)}


def steps(path):
    return len(re.findall(r"^Time = ", open(f"{path}/log.celFoam").read(), re.M))


def rate(path):
    txt = open(f"{path}/log.celFoam").read()
    m = np.array([float(x) for x in re.findall(r"Discharged mass = (\S+)", txt)]) / hc.depth()
    t = np.array([0.0] + [float(x) for x in pc.times(path)])
    n = min(len(t), len(m))
    t, m = t[:n], m[:n]
    s = (t >= WINDOW[0] - 1e-9) & (t <= WINDOW[1] + 1e-9)
    return float(np.polyfit(t[s], m[s], 1)[0])


def main():
    work = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "run")
    os.makedirs(work, exist_ok=True)
    nc = hc.ncells()
    rows, ok, W, secs, nst, notes = [], True, {}, {}, {}, []
    for gas in (False, True):
        serial = None
        for n, dec in DECOMP.items():
            path = f"{work}/{'gas' if gas else 'dry'}_np{n}"
            hc.make_hopper(path, d=D, end=END, write_interval=WRITE, gas=gas, decomp=dec,
                           rules=False, extra_gas="pTolerance 1e-12;\n")
            secs[gas, n] = pc.run(path, n, set_fields=True)
            nst[gas, n] = steps(path)
            W[gas, n] = rate(path)
            ts = pc.times(path)
            if serial is None:
                serial = path
                pmax = (max(float(np.abs(pc.field(path, t, "p", nc)).max()) for t in ts)
                        if gas else 1.0)
            du = dp = dbulk = 0.0
            where = set()
            for t in ts:
                ua, ub = pc.field(serial, t, "U", nc), pc.field(path, t, "U", nc)
                e = np.abs(ua - ub).max(axis=1) / max(np.abs(ua).max(), 1e-300)
                du = max(du, float(e.max()))
                al = pc.field(serial, t, "alpha", nc)
                dbulk = max(dbulk, float(e[al > 0.999].max()))
                where |= {(int(i), round(float(al[i]), 2)) for i in np.where(e > 1e-10)[0]}
                if gas:
                    dp = max(dp, float(np.abs(pc.field(serial, t, "p", nc)
                                              - pc.field(path, t, "p", nc)).max()) / pmax)
            drift = sv = 0.0
            if gas:
                M, N, s, _ = pc.mass_log(path)
                drift = float(np.abs(M - N - (M[0] - N[0])).max() / M[0])
                sv = float(s.mean())
            good = du <= 1e-10 and dp <= 1e-10 and drift <= 1e-12
            ok = ok and good
            rows.append((("gas" if gas else "dry"), n, nst[gas, n], secs[gas, n], W[gas, n],
                         max(du, dp), drift, sv, good))
            if where:
                notes.append(f"{'gas' if gas else 'dry'} on {n} processors: U differs from "
                             f"serial by up to {du:.0e} in cells (index, alpha) "
                             f"{sorted(where)[:6]}; full cells within {dbulk:.0e}")
            print(f"{'gas' if gas else 'dry'} np {n}: rate {W[gas, n]:.3f} kg/(m s); vs serial "
                  f"{max(du, dp):.0e}; mass {drift:.0e}; {nst[gas, n]} steps {secs[gas, n]:.1f} s: "
                  f"{'PASS' if good else 'FAIL'}", flush=True)
    slow = W[True, 1] / W[False, 1]
    g_ok = slow <= 0.8
    ok = ok and g_ok
    # Crewdson's outlet gradient and Muite's N (report)
    g = f"{work}/gas_np1"
    nh, h0 = hc.outlet_row()
    G = np.mean([float(pc.field(g, t, "p", nc)[: hc.NCOL].mean()) / (h0 / 2)
                 for t in pc.times(g) if float(t) >= WINDOW[0]])
    cw = math.sqrt(max(0.0, 1 + G / (hc.RHO * hc.G_ACC)))
    ut = gr.terminal_velocity(D, pc.RHOP, gr.air_density(), gr.air_viscosity())
    N = ut / (W[True, 1] / (hc.RHO * 2 * hc.RO))
    cost = {n: (secs[True, n] / nst[True, n]) / (secs[False, n] / nst[False, n]) for n in DECOMP}
    print(f"gas on / off: {slow:.3f} (limit 0.8): {'PASS' if g_ok else 'FAIL'}; Crewdson "
          f"sqrt(1+G/rho g) {cw:.3f}; Muite N {N:.2f}; cost per step "
          + ", ".join(f"{n} proc {c:.2f}x" for n, c in cost.items()))

    with open(os.path.join(HERE, "times.csv"), "w") as f:
        f.write("gas,nprocs,steps,wall_s,rate_kg_per_m_s,vs_serial,mass_drift,solves_per_step,pass\n")
        for r in rows:
            f.write("%s,%d,%d,%.3f,%.4f,%.1e,%.1e,%.2f,%s\n" % r)
    with open(os.path.join(HERE, "results.md"), "w") as f:
        f.write("# Gas validation pack test 4: fine-powder hopper, gas on and off\n\n"
                f"Plane-flow hopper, slot {2 * hc.RO * 1e3:g} mm, walls {hc.BETA:g} deg from "
                f"vertical, {nc} body-fitted cells; quartz powder d {D * 1e3:g} mm, gas fraction "
                f"{hc.E_GAS}, plain Mohr-Coulomb phi 30, psi 1, E 1 MPa; drag ergunWenYu; outlet "
                f"opens at {hc.OPEN_TIME} s, rate from {WINDOW[0]} to {WINDOW[1]} s. Written by "
                "`fine_hopper.py`.\n\n"
                f"Gas on against gas off (serial): **{W[True, 1]:.2f} against {W[False, 1]:.2f} "
                f"kg/(m s), {100 * (slow - 1):+.1f} %** (limit: at least 20 % slower). "
                f"Crewdson's sqrt(1 + G/(rho_b g)) from the run's outlet gradient: {cw:.3f}. "
                f"Muite's N = U_t/U_out: {N:.2f} (the gas matters below about 10).\n\n"
                "| run | processors | steps | wall [s] | rate [kg/(m s)] | vs serial | gas mass | "
                "solves per step |\n| --- | --- | --- | --- | --- | --- | --- | --- |\n")
        for r in rows:
            f.write("| %s | %d | %d | %.1f | %.3f | %.0e | %.0e | %.2f |\n" % r[:8])
        f.write("\nCost per step, gas on against off: "
                + ", ".join(f"{n} processor(s) {c:.2f}x" for n, c in cost.items())
                + f" (F.9 target 1.5x).\n\n"
                + "".join(f"- {x}\n" for x in notes)
                + f"\nAll checks {'pass' if ok else 'FAIL'}.\n")
    for x in notes:
        print(x)
    print(f"\nFine hopper: {'ALL PASS' if ok else 'FAILED'}; times.csv and results.md written")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
