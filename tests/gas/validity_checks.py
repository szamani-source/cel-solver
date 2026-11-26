"""Range of validity of the gas coupling: the lift check (decision 0135).

The gas coupling holds the grains at their packed density, so it is valid
while the grains' own stress carries their weight less the gas's share.
The lift check reports the share of the grain mass whose weight the gas
carries and warns once it stays above liftWarn (1 %) for liftTime (0.01 s).

B  Fluidisation bed (pack test 2 case, 8 x 40, air ramped to 1.4 U_mf over
   0.28 s): no warning before 0.9 U_mf; the warning comes between 0.9 and
   1.1 U_mf (the time at which the inflow passes them); share above 50 %
   at 0.26 s.
H  Fine-powder hopper with gas (pack test 4 case, outlet opens at 0.15 s,
   to 0.5 s): no warning, share at most 1 % at every write.

    python3 validity_checks.py WORKDIR
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import fluid_cases as fl  # noqa: E402
import hopper_gas_cases as hc  # noqa: E402
import pore_cases as pc  # noqa: E402

SHARE = re.compile(r"Gas lift: the gas carries the weight of (\S+) %")
WARN = re.compile(r"since t = (\S+) s the gas has carried")


def log(path):
    return open(f"{path}/log.celFoam").read()


def main():
    work = sys.argv[1]
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + msg)
        ok = ok and cond

    # B: bed
    path = f"{work}/bed"
    fl.make_bed(path, end=0.26, write_interval=0.01, nx=8)
    pc.run(path, 1)
    txt = log(path)
    shares = [float(s) for s in SHARE.findall(txt)]
    warn = WARN.findall(txt)
    u_wy = fl.U_WY
    t_lo = fl.T_RAMP * 0.9 * u_wy / fl.U_END
    t_hi = fl.T_RAMP * 1.1 * u_wy / fl.U_END
    print(f"B  bed: U_mf (Wen-Yu) {u_wy:.5f} m/s reached at {fl.T_RAMP * u_wy / fl.U_END:.4f} s; "
          f"0.9 / 1.1 U_mf at {t_lo:.4f} / {t_hi:.4f} s")
    print("   share at writes [%]: " + " ".join(f"{s:.3g}" for s in shares))
    check(len(warn) == 1, f"bed warns once (got {len(warn)})")
    if warn:
        ts = float(warn[0])
        check(t_lo <= ts <= t_hi, f"bed carried from {ts:.4f} s, within 0.9-1.1 U_mf")
    check(bool(shares) and shares[-1] > 50, f"bed share at 0.26 s {shares[-1] if shares else 0:.1f} % > 50 %")

    # H: hopper
    path = f"{work}/hopper"
    hc.make_hopper(path, d=1e-4, end=0.5, write_interval=0.025, gas=True, decomp=(1, 1, 1),
                   rules=False, extra_gas="pTolerance 1e-12;\n")
    pc.run(path, 1, set_fields=True)
    txt = log(path)
    shares = [float(s) for s in SHARE.findall(txt)]
    print("H  hopper: share at writes [%]: " + " ".join(f"{s:.3g}" for s in shares))
    check(not WARN.findall(txt), "hopper does not warn")
    check(bool(shares) and max(shares) <= 1, f"hopper share at most 1 % (max {max(shares) if shares else 0:.3g} %)")

    print("All checks PASS." if ok else "Some checks FAIL.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
