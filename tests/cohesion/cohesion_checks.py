"""WP 3.5 done-when: cohesion in the T0 tank. With cohesion c the material
arches over a narrow outlet and the flow stops; above a critical outlet width
it flows (decision 0054).

Tank as WP 3.3 (tests/outlet/outlet_checks.py): W = 0.4 m, 10 mm cells, walls
and floor 20 deg, outlet centred in the floor, opens at 0.3 s after a damped
settle. Fill H = 0.3 m (the arch does not depend on the head, so half of T0's
fill keeps the runs short). T0 rule set (decision 0050): Mohr-Coulomb phi 35,
psi 0, mu(I) placeholders (decision 0037), flowImpedance, normalCap, shearCap,
spreadLayers; E = 1 MPa (strength, not stiffness, sets the arch; the manual
test checks 10 MPa). Cohesion c with the default tensile strength c cot(phi).

Reference: Jenike's critical width of a plane (slot) outlet,
B_c = sigma_c / (rho g), with the Mohr-Coulomb unconfined strength
sigma_c = 2 c cos(phi) / (1 - sin(phi)): 109 mm at c = 500 Pa.

A run "stops" when its discharge rate over the last 0.15 s (of 0.5 s after
opening) is below 5 % of the Beverloo rate with C = 1 at that width,
rho sqrt(g) D^1.5; it "flows" when that rate is above 25 % of it.

Checks:
1. c = 500 Pa: D = 80 and 100 mm stop, D = 120 and 140 mm flow, so the
   critical width lies between 100 and 120 mm (Jenike 109 mm, within 10 %).
2. Control: the same tank with c = 0 flows at D = 80 mm (the stop is cohesion).
3. Every run: tank plus discharged mass constant to 1e-10, no stress in empty
   cells, and in the stopping runs some material leaves first (the material
   under the arch tears off at the tensile strength instead of hanging).

    python3 cohesion_checks.py WORKDIR [--procs N]       # CI, all runs at once
    python3 cohesion_checks.py WORKDIR --case C D [E_MPa] [T_DIS]   # one run
"""

import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "outlet"))
import outlet_checks as oc  # noqa: E402  (setup, history columns)

H = 0.30
T_DIS = 0.5
LAST = 0.15                     # rate window at the end of the run [s]
STOP, FLOW = 0.05, 0.25         # of the Beverloo rate with C = 1
T0_RULES = ("mu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650.0;\n"
            "flowImpedance yes;\nnormalCap yes;\n")


def jenike(c, phi=oc.PHI, rho=oc.RHO):
    s = math.sin(math.radians(phi))
    return 2 * c * math.cos(math.radians(phi)) / (1 - s) / (rho * oc.G_ACC)


def beverloo(D):
    return oc.RHO * math.sqrt(oc.G_ACC) * D ** 1.5


def run(work, c, D, E_mpa=1.0, t_dis=T_DIS):
    """One tank; returns its rates and checks. Not thread safe in setup
    (outlet_checks keeps its settings in module globals), so setup is done
    serially by the caller and the solver runs in parallel."""
    oc.E, oc.D, oc.H, oc.T_DIS = E_mpa * 1e6, D, H, t_dis
    path = oc.setup(work, f"c{c:g}_D{D * 1000:g}_E{E_mpa:g}")
    mp = f"{path}/constant/materialProperties"
    txt = open(mp).read().replace("c 0;", f"c {c!r};")
    if c > 0:
        txt = txt.replace("sigmaT 0;\n", "")
    open(mp, "w").write(txt + T0_RULES)
    return path


def result(path, c, D, t_dis=T_DIS):
    h = np.loadtxt(f"{path}/history.dat")
    t, m, out = h[:, oc.T], h[:, oc.MASS] / oc.CELL, h[:, oc.OUT] / oc.CELL
    i_open = int(np.where(t <= oc.T_OPEN + 1e-12)[0][-1])
    m_open = m[i_open]
    sel = t >= t[-1] - LAST - 1e-12
    rate = -np.polyfit(t[sel], m[sel], 1)[0]
    return dict(c=c, D=D, rate=rate, frac=rate / beverloo(D), out=out[-1] / m_open,
                cons=np.abs((m + out) / m[0] - 1.0).max(), s_empty=h[:, oc.SEMPTY].max(),
                t_end=t[-1])


def main():
    work = sys.argv[1]
    args = sys.argv[2:]
    t_all = time.perf_counter()
    if "--case" in args:
        i = args.index("--case")
        vals = [float(v) for v in args[i + 1:]]
        c, D = vals[0], vals[1]
        E_mpa = vals[2] if len(vals) > 2 else 1.0
        t_dis = vals[3] if len(vals) > 3 else T_DIS
        path = run(work, c, D, E_mpa, t_dis)
        sec = oc.vc.run(path)
        r = result(path, c, D, t_dis)
        print(f"c {c:g} Pa, D {D * 1e3:g} mm, E {E_mpa:g} MPa ({sec:.0f} s): rate over the "
              f"last {LAST} s {r['rate']:.2f} kg/s/m = {r['frac']:.3f} of Beverloo C = 1; "
              f"discharged {r['out']:.4f} of the fill; Jenike B_c {jenike(c) * 1e3:.0f} mm")
        return 0
    procs = int(args[args.index("--procs") + 1]) if "--procs" in args else os.cpu_count()

    runs = [(500.0, 0.08), (500.0, 0.10), (500.0, 0.12), (500.0, 0.14), (0.0, 0.08)]
    paths = [run(work, c, D) for c, D in runs]
    with ThreadPoolExecutor(max_workers=max(1, procs)) as ex:
        secs = list(ex.map(oc.vc.run, paths))
    res = [result(p, c, D) for p, (c, D) in zip(paths, runs)]

    print(f"T0 tank W = {oc.W} m, fill {H} m, 10 mm cells, E = 1 MPa, T0 rule set, "
          f"walls 20 deg; outlet opens at {oc.T_OPEN} s, run {T_DIS} s after")
    print(f"Jenike B_c = sigma_c / (rho g): {jenike(500) * 1e3:.0f} mm at c = 500 Pa")
    print(f"stop: rate over the last {LAST} s < {STOP} of Beverloo C = 1; "
          f"flow: > {FLOW}\n")
    print("     c [Pa]  D [mm]  rate [kg/s/m]  / Beverloo C=1  discharged  wall [s]")
    for r, s in zip(res, secs):
        print(f"    {r['c']:7g}  {r['D'] * 1e3:6g}  {r['rate']:13.2f}  {r['frac']:14.3f}  "
              f"{r['out']:10.4f}  {s:7.0f}")
    print()
    ok = []

    def report(name, good, detail):
        ok.append(good)
        print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}")

    st = [r for r in res if r["c"] > 0 and r["D"] <= 0.1001]
    fl = [r for r in res if r["c"] > 0 and r["D"] >= 0.1199]
    report("1 critical width", all(r["frac"] < STOP for r in st)
           and all(r["frac"] > FLOW for r in fl),
           "c = 500 Pa stops at " + ", ".join(f"{r['D'] * 1e3:g}" for r in st)
           + " mm and flows at " + ", ".join(f"{r['D'] * 1e3:g}" for r in fl)
           + f" mm (B_c between 100 and 120 mm, Jenike {jenike(500) * 1e3:.0f} mm)")
    c0 = res[-1]
    report("2 cohesionless control", c0["frac"] > FLOW,
           f"c = 0 at 80 mm flows at {c0['frac']:.3f} of Beverloo C = 1 (> {FLOW})")
    first = ", ".join(f"{r['out']:.4f}" for r in st)
    report("3 mass, void, tearing",
           all(r["cons"] <= 1e-10 and r["s_empty"] == 0.0 for r in res)
           and all(r["out"] > 0 for r in st),
           f"mass constant to {max(r['cons'] for r in res):.0e} (<= 1e-10), largest stress "
           f"in empty cells {max(r['s_empty'] for r in res):.0e} Pa (0), stopping runs "
           f"discharged {first} of the fill first (> 0)")
    good = all(ok)
    print(f"\nCohesion checks: {'ALL PASS' if good else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if good else 1


if __name__ == "__main__":
    sys.exit(main())
