"""WP 3.4: the T0 tank discharge rate in celFoam's cycle at E = 0.1, 1 and
10 MPa, and a parallel run against serial (decision 0034).

Same tank, outlet and settings as WP 3.3 (tests/outlet/outlet_checks.py,
decision 0033): 10 mm cells, W 0.4 m, H 0.6 m, D 0.06 m, Mohr-Coulomb c = 0,
phi 35, walls and floor 20 deg, shear cap and rule 9 on (the collapse code,
decision 0032), neither rule 12 nor 13 (decision 0029). Only E changes.

Checks (1 and 3 by default, as tests/test_discharge.sh in CI; 2 with --sweep):
1. Mass balance (tank plus discharged constant to 1e-10) and no stress in
   empty cells, every run.
2. Rate at E = 1 and 10 MPa within 10 percent of the rate at 0.1 MPa (one code
   for collapse and discharge, decision 0029). NOT MET (decision 0034):
   10 mm gives 74.1 / 65.7 / 70.9 kg/s/m (1 MPa -11.3 %), 5 mm on four
   processors 78.0 / 55.0 at 0.1 / 1 MPa (-29 %). This check is the done-when
   of the follow-up package and runs only by hand (tests/manual/).
3. Two processors (split at x = 0.2 m, through the middle of the outlet) at
   E = 0.1 MPa: mass at opening equal to serial to 1e-10, rate within 2 percent
   of serial (a plastic flow parts by round-off growth, decision 0032, so the
   flow is compared by its rate, not cell by cell).

Run after celVoidCheck is built:

    python3 discharge_checks.py WORKDIR                 # checks 1 and 3 (CI, ~50 s)
    python3 discharge_checks.py WORKDIR --sweep [--cell 0.005] [--procs 4] [E_MPa ...]
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "outlet"))
import outlet_checks as oc  # noqa: E402  (setup, steady_rate, history columns)

E_LIST = [0.1, 1.0, 10.0]          # MPa
T_DIS = 1.0                        # as WP 3.3; at 5 mm and 1 MPa the window ends at 0.9 s
TOL_E = 0.10
TOL_PAR = 0.02
# prototype without rules 12 and 13 (decision 0027, 10 mm): 81 / 73 / 54 kg/s/m
PROTO = {0.1: 81.35, 1.0: 73.0, 10.0: 54.0}


def case(work, E_mpa, nProcs=1):
    oc.E = E_mpa * 1e6
    oc.T_DIS = T_DIS
    name = f"E{E_mpa:g}" + (f"p{nProcs}" if nProcs > 1 else "")
    path = oc.setup(work, name)
    if nProcs > 1:
        oc.cs.write(f"{path}/system/decomposeParDict",
                    oc.cs.HEADER.format(cls="dictionary", obj="decomposeParDict")
                    + f"numberOfSubdomains {nProcs};\nmethod simple;\n"
                    f"coeffs {{ n ({nProcs} 1 1); }}\n")
    sec = oc.vc.run(path, nProcs)
    h = np.loadtxt(f"{path}/history.dat")
    t, m, out = h[:, oc.T], h[:, oc.MASS] / oc.CELL, h[:, oc.OUT] / oc.CELL
    i_open = int(np.where(t <= oc.T_OPEN + 1e-12)[0][-1])
    m_open = m[i_open]
    rate, win = oc.steady_rate(t[i_open:] - oc.T_OPEN, m[i_open:], m_open)
    C = rate / (oc.RHO * math.sqrt(oc.G_ACC) * oc.D ** 1.5)
    return dict(E=E_mpa, rate=rate, C=C, win=win, m_open=m_open,
                cons=np.abs((m + out) / m[0] - 1.0).max(), s_empty=h[:, oc.SEMPTY].max(),
                steps=int(h[-1, oc.STEP]), sec=sec, m_end=m[-1] / m_open)


def main():
    work = sys.argv[1]
    args = sys.argv[2:]
    sweep = "--sweep" in args
    procs, cell = 1, oc.CELL
    if "--procs" in args:
        procs = int(args[args.index("--procs") + 1])
    if "--cell" in args:
        cell = float(args[args.index("--cell") + 1])
        oc.CELL = oc.vc.CELL = cell
    nums = [a for i, a in enumerate(args) if not a.startswith("--")
            and (i == 0 or args[i - 1] not in ("--procs", "--cell"))]
    elist = ([float(a) for a in nums] or E_LIST) if sweep else [0.1]
    t_all = time.perf_counter()
    print(f"T0 tank W = {oc.W} m, H = {oc.H} m, outlet D = {oc.D} m, "
          f"{1000 * cell:g} mm cells, Mohr-Coulomb c = 0, phi {oc.PHI:g}; walls {oc.DELTA_W:g} deg; opens at "
          f"{oc.T_OPEN} s;\nshear cap and rule 9 on, neither rule 12 nor 13\n")
    res, runs = [], {}

    def report(name, good, detail):
        res.append(good)
        print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}")

    for E_mpa in elist:
        r = case(work, E_mpa, procs)
        runs[E_mpa] = r
        print(f"  E = {E_mpa:g} MPa: {r['rate']:.1f} kg/s/m, C = {r['C']:.3f} over t = "
              f"{r['win'][0]:.3f} to {r['win'][1]:.3f} s after opening (prototype "
              f"{PROTO.get(E_mpa, float('nan')):g}); {r['steps']} steps, {r['sec']:.0f} s"
              + (f" on {procs} processors" if procs > 1 else ""))
        sys.stdout.flush()
    print()
    report("1 mass and void", all(r["cons"] <= 1e-10 and r["s_empty"] == 0.0
                                  for r in runs.values()),
           "tank plus discharged to " + ", ".join(f"{r['cons']:.0e}" for r in runs.values())
           + "; stress in empty cells " + ", ".join(f"{r['s_empty']:.0e}"
                                                   for r in runs.values()))
    if len(runs) > 1:
        r0 = runs[elist[0]]["rate"]
        dev = {e: r["rate"] / r0 - 1.0 for e, r in runs.items()}
        report("2 rate against E", all(np.isfinite(d) and abs(d) <= TOL_E for d in dev.values()),
               ", ".join(f"{e:g} MPa {100 * d:+.1f} %" for e, d in dev.items())
               + f" of the {elist[0]:g} MPa rate (within {100 * TOL_E:g} %)")
    if not sweep:
        p = case(work, 0.1, 2)
        s = runs[0.1]
        dmo = abs(p["m_open"] / s["m_open"] - 1.0)
        dr = p["rate"] / s["rate"] - 1.0
        report("3 parallel", dmo <= 1e-10 and abs(dr) <= TOL_PAR and p["cons"] <= 1e-10,
               f"2 processors: mass at opening {dmo:.0e} from serial (<= 1e-10), rate "
               f"{p['rate']:.2f} against {s['rate']:.2f} kg/s/m ({100 * dr:+.2f} %, within "
               f"{100 * TOL_PAR:g} %), mass balance {p['cons']:.0e}; {p['sec']:.0f} s")
    ok = all(res)
    print(f"\nDischarge (C++) checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
