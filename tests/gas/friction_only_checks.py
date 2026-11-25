"""WP F.8g: pore pressure that only lowers the friction (gasProperties
frictionOnly; decision 0133; Gueugneau et al. 2017, Aravena et al. 2021).

Coarse dam break (the F.8 column, 20 mm cells, 150 x 23, to 0.4 s):

I  Full support kept (diffusivity 1e-6, frictionOnly): the grains have no
   friction, so the column spreads like the dry case with friction removed
   (phi 0.5, mu2 tan 1 deg, floor 0): front from the gate within 15 % of it
   at 0.2, 0.3 and 0.4 s. (Before F.8g the pore pressure acted as a force:
   the interior stayed still and only the edge cell was blown off.)
S  The same with the support lost (diffusivity 1e-6, initial p = 0): the
   front is behind I's at 0.4 s (friction acts).
P  2 processors (cut at x = 0.1 m) against serial (report): a cohesionless
   flow at the Mohr-Coulomb apex parts serial from parallel from round-off;
   the dry frictionless case parts the same way (known issue, decisions 0132,
   0133; F.9).

tests/test_friction_only.sh, about 1 minute. Result 2026-11-25 (F.8h: the
support ratio carried; flow impedance and the wall predictor see the
lowered friction): I within 5.6 % (before F.8h 25-29 % behind), S passes.

    python3 friction_only_checks.py WORKDIR
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import collapse_cases as cc  # noqa: E402
import collapse_checks as ck  # noqa: E402
import gas_flow_checks as gf  # noqa: E402
import pore_cases as pc  # noqa: E402

cc.CELL, cc.NX, cc.NY = 0.02, 150, 23
cc.N = cc.NX * cc.NY
END = 0.4
TIMES = (0.2, 0.3, 0.4)
res = []


def report(name, good, detail):
    res.append(bool(good))
    print(f"  {name}: {detail}: {'PASS' if good else 'FAIL'}", flush=True)


def fronts(path):
    h = cc.deposit_history(path)
    return np.array([np.interp(t, h[:, 0], h[:, 1]) - cc.X0 for t in TIMES])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    base_mc = cc.MC
    # I: full support, friction only
    wi = f"{work}/support"
    cc.make(wi, kappa=1e-6, friction_only=True, end=END, extra_gas="pTolerance 1e-9;\n")
    cc.run(wi)
    # reference: dry, friction removed
    wd = f"{work}/frictionless"
    cc.MC = base_mc.replace(f"phi {cc.PHI!r};", "phi 0.5;").replace(
        f"mu2 {cc.MU2!r};", f"mu2 {math.tan(math.radians(1.0))!r};")
    floor = cc.FLOOR
    cc.FLOOR = 0.0
    cc.make(wd, gas=False, cap=False, end=END)
    cc.FLOOR = floor
    cc.MC = base_mc
    cc.run(wd)
    fi, fd = fronts(wi), fronts(wd)
    dev = float(np.max(np.abs(fi / fd - 1)))
    report("I full support spreads as frictionless", dev <= 0.15,
           "front from the gate " + ", ".join(f"{a:.3f}/{b:.3f} m" for a, b in zip(fi, fd))
           + f" at {TIMES} s (support/frictionless), largest {100 * dev:.1f} % (<= 15 %)")
    # S: no support
    ws = f"{work}/nosupport"
    cc.make(ws, kappa=1e-6, friction_only=True, end=END, extra_gas="pTolerance 1e-9;\n")
    txt = open(f"{ws}/0/p").read()
    pc.write(f"{ws}/0/p", re.sub(r"internalField nonuniform List<scalar> \d+\n\(.*?\n\);",
                                 "internalField uniform 0;", txt, count=1, flags=re.S))
    cc.run(ws)
    fs = fronts(ws)
    report("S support lost slows the front", fs[-1] < 0.8 * fi[-1],
           f"front at {TIMES[-1]} s {fs[-1]:.3f} m against {fi[-1]:.3f} m with support")
    # P: 2 processors
    w2 = f"{work}/support2"
    cc.make(w2, kappa=1e-6, friction_only=True, end=END, extra_gas="pTolerance 1e-9;\n")
    ck.write_manual(w2, 2)
    cc.run(w2, 2)
    d = max(gf.rel(gf.vfield(wi, t, f, cc.N), gf.vfield(w2, t, f, cc.N))
            for t in pc.times(wi) for f in ("p", "U"))
    print(f"  P 2 processors (report, apex round-off, F.9): p, U within {d:.1e} of their largest")
    print(f"\ntest_friction_only: {'PASS' if all(res) else 'FAIL'} ({sum(res)}/{len(res)}), "
          f"{time.perf_counter() - t0:.0f} s")
    return 0 if all(res) else 1


if __name__ == "__main__":
    sys.exit(main())
