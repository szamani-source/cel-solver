"""WP F.8k done-when (decision 0138): the dry dam break of the fluidised
dam-break case (tests/gas/collapse_cases.py: column 0.2 x 0.4 m, floor
19 deg, E 1 MPa) with plain Mohr-Coulomb (phi 24 deg, c 0, no mu(I)),
yieldShearJump no and slowFlowJump yes converges with the mesh at the
default Courant number 0.4.

Run-out at 1.4 s (from the gate, in h0): the front is the end of the
continuous deposit at least 20 mm thick, two cells of the coarsest mesh
(Sia, 2026-11-27). Each mesh ends the deposit's thin tail where it is one
cell thick, so a thinner threshold measures the mesh, not the scheme
(WP F.8k: at 5 mm thickness 2.30 / 2.46 / 2.84 h0). Passes when the 10, 5
and 2.5 mm run-outs lie within 5 % of the 2.5 mm one (WP F.8k: 1.58 / 1.69 /
1.62 h0; without slowFlowJump 1.40 / 1.67 / 1.75).


    python3 collapse_convergence.py WORKDIR [cell_mm ...] [--var slowFlowJump]
"""

import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "gas"))
import collapse_cases as cc  # noqa: E402

T_END, THICK, TOL = 1.4, 0.02, 0.05
PROCS = {10: 1, 5: 2, 2.5: 4}


def setup(cell):
    cc.CELL = cell
    cc.NX, cc.NY = int(round(3 / cell)), int(round(0.45 / cell))
    cc.N = cc.NX * cc.NY


def material(extra):
    m = cc.MC.replace("flowImpedance yes;", "flowImpedance no;")
    m = re.sub(r"mu2 [^\n]*\n|I0 [^\n]*\n", "", m)
    return m + "\nyieldShearJump no;\n" + extra


def runout(path):
    rows = []
    for t in cc.pc.times(path):
        hc, _ = cc.columns(path, t)
        k = int(round(cc.X0 / cc.CELL))
        while k < cc.NX and hc[k] >= THICK:
            k += 1
        rows.append((float(t), (k * cc.CELL - cc.X0) / cc.H0))
    return rows


def main():
    cc.FLOOR = 19.0     # the runs of WP F.8k (floor 22 deg came later)
    args = sys.argv[1:]
    work = args.pop(0)
    extra = "slowFlowJump yes;"
    if "--base" in args:
        args.remove("--base")
        extra = ""
    cells = [float(a) for a in args] or [10, 5, 2.5]
    base = cc.MC
    res = {}
    t0 = time.perf_counter()
    for mm in cells:
        setup(mm / 1000)
        cc.MC = material(extra)
        path = f"{work}/conv_{mm:g}"
        if not os.path.isdir(f"{path}/system"):
            cc.make(path, gas=False, decomp=(PROCS.get(mm, 1), 1, 1), end=T_END, slow=False)
        secs = cc.run(path, PROCS.get(mm, 1))
        cc.MC = base
        h = runout(path)
        res[mm] = h[-1][1]
        print(f"  {mm:g} mm: run-out {res[mm]:.3f} h0 at {h[-1][0]:g} s; "
              + ", ".join(f"{x:.2f}" for t, x in h[9::10]) + f"; {secs:.0f} s", flush=True)
    ref = res[min(res)]
    worst = max(abs(v / ref - 1) for v in res.values())
    ok = worst <= TOL
    print(f"\ncollapse_convergence: largest difference from {min(res):g} mm {100 * worst:.1f} % "
          f"(<= {100 * TOL:g} %): {'PASS' if ok else 'FAIL'}, {time.perf_counter() - t0:.0f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
