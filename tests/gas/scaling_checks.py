"""WP F.9: cost of the gas and parallel speed-up (decision 0141).

The F.5 plane-flow hopper refined to 64 columns (8,704 body-fitted cells),
0.2 mm grains (the airImpededHopper tutorial), T0 rule set, outlet open from
the start so that the grains flow and the gas works from the first step.
STEPS steps each (the end time is STEPS times the first CFL step, which
stays fixed while the material is near rest), wall time per step from
log.celFoam (ExecutionTime of the last step less that of the first).

Checks (F.9 done-when, Sia 2026-11-24, decision 0125 item 8):
C  gas on costs at most 1.5 times gas off per step, serial;
S  4 processors run the gas case at least 2.5 times faster than serial
   (v1: 2.6 at 620 cells per processor, WP 6.3a). 2 processors reported.
Needs 4 free cores: run it on an idle machine, one run at a time.

    python3 scaling_checks.py WORKDIR [STEPS]
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hopper_gas_cases as hc  # noqa: E402
import pore_cases as pc  # noqa: E402

NCOL, D = 64, 2e-4
failed = []


def report(name, good, detail):
    print(f"  {'PASS' if good else 'FAIL'} {name}: {detail}", flush=True)
    if not good:
        failed.append(name)


def per_step(path):
    log = open(f"{path}/log.celFoam").read()
    ex = [float(x) for x in re.findall(r"^ExecutionTime = (\S+) s", log, re.M)]
    return (ex[-1] - ex[0]) / (len(ex) - 1), len(ex)


def case(path, *, gas, n, end):
    dc = {1: (1, 1, 1), 2: (2, 1, 1), 4: (2, 2, 1)}[n]
    # the user default pTolerance (1e-4 Pa, decision 0141), not the tests' 1e-14 p0
    hc.make_hopper(path, d=D, end=end, write_interval=end, gas=gas, decomp=dc, ncol=NCOL,
                   open_time=0.0, extra_gas="pTolerance 0.0001;\n")
    pc.run(path, n, set_fields=True)
    return per_step(path)


def main(work, steps):
    os.makedirs(work, exist_ok=True)
    print(f"Hopper {hc.ncells(NCOL)} cells, d {D * 1e3:g} mm, T0 rules, outlet open, {steps} steps")
    # first CFL step
    probe = f"{work}/probe"
    hc.make_hopper(probe, d=D, end=1e-9, write_interval=1.0, gas=False, ncol=NCOL, open_time=0.0)
    pc.run(probe, set_fields=True)
    dt = float(re.findall(r"deltaT = (\S+)", open(f"{probe}/log.celFoam").read())[0])
    end = dt * (steps + 0.5)
    t = {}
    for name, gas, n in (("dry1", False, 1), ("gas1", True, 1), ("gas2", True, 2), ("gas4", True, 4)):
        t[name], k = case(f"{work}/{name}", gas=gas, n=n, end=end)
        print(f"  {name}: {t[name] * 1e3:.2f} ms per step ({k} steps)", flush=True)
    r = t["gas1"] / t["dry1"]
    report("gas on at most 1.5 times gas off per step", r <= 1.5,
           f"{t['gas1'] * 1e3:.2f} ms against {t['dry1'] * 1e3:.2f} ms ({r:.2f} times)")
    print(f"  2 processors: {t['gas1'] / t['gas2']:.2f} times serial")
    s = t["gas1"] / t["gas4"]
    report("4 processors at least 2.5 times serial", s >= 2.5,
           f"{s:.2f} times ({hc.ncells(NCOL) // 4} cells per processor)")
    if failed:
        print("FAILED: " + ", ".join(failed))
        return 1
    print("test_gas_scaling: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 200))
