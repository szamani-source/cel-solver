"""WP F.9 done-when: the gas runs restart and run in parallel (decision 0141).

Cases from earlier packages:
H  the F.5 air-impeded hopper (0.1 mm grains, outlet opens at 0.15 s);
D  the F.8 fluidised dam break (route B, frictionOnly, slowFlowJump, 10 mm).

A  Restart: each case run again from a copy holding only a middle write time
   (startFrom latestTime; ascii at 17 digits round-trips every double)
   matches the continuous run at the end to 1e-12 of each field's largest
   value, every field written (H: T0 rule set, restarted at 0.2 s after the
   outlet opened, compared at 0.25 s; D: restarted at 0.1 s, compared at
   0.14 s). celFoam writes the gas solver's state with the fields
   (<time>/uniform/poreGas: gas totals, last step's pressure and length,
   frictionOnly's support pressures, lift clock); without it the dam break
   parted by order 1 after a restart.
B  Processors: H (plain Mohr-Coulomb, whose rules have no switch that
   round-off flips, decision 0080) on 2 (2 x 1) and 4 (2 x 2) processors
   to 0.25 s: p, U within 1e-10 of their largest, gas mass to 1e-12.
   D on 2 and 4 processors with the cuts through the column (x = 0.1 m;
   4: also y = 0.1 m): the time up to which p, U, sigma, alpha stay within
   1e-10 is REPORTED, not checked. The dam break starts at zero effective
   stress, where round-off decides yield at the free face (a serial run
   with E changed by 1e-13 departs as fast), and dry its free-surface rules
   (flowImpedance, fullBlend, overfillStress) flip on round-off, so the
   done-when's 1e-10 over the whole run is not met for it (decision 0141); the full-run comparison stays
   in tests/manual/test_fluidised_collapse.sh (1e-10, fails).

    python3 restart_mpi_checks.py WORKDIR
"""

import os
import re
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import collapse_cases as cc  # noqa: E402
import collapse_checks as ck  # noqa: E402
import hopper_gas_cases as hc  # noqa: E402
import pore_cases as pc  # noqa: E402

failed = []


def report(name, good, detail):
    print(f"  {'PASS' if good else 'FAIL'} {name}: {detail}", flush=True)
    if not good:
        failed.append(name)


def restart_copy(src, dst, t_restart):
    """Copy of src with only the write times up to t_restart, resumed with
    startFrom latestTime (serial)."""
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    for d in os.listdir(dst):
        if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > t_restart + 1e-12:
            shutil.rmtree(os.path.join(dst, d))
        elif d.startswith("processor") or d.startswith("log.celFoam"):
            p = os.path.join(dst, d)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    cd = f"{dst}/system/controlDict"
    txt = open(cd).read()
    open(cd, "w").write(re.sub(r"startFrom\s+\w+;", "startFrom       latestTime;", txt))
    pc.sh(dst, ["celFoam"], "celFoam")


def fields(path, t):
    return sorted(f for f in os.listdir(f"{path}/{t}") if os.path.isfile(f"{path}/{t}/{f}"))


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.abs(a - b).max() / max(np.abs(a).max(), 1e-300))


def compare(a, b, t, names=None):
    """Largest difference of each field relative to its largest value."""
    out = {}
    for f in names or fields(a, t):
        try:
            out[f] = rel(pc.field(a, t, f), pc.field(b, t, f))
        except (ValueError, AttributeError, IndexError):
            continue        # not a plain internal field
    return out


def check_restart(work):
    print("A: restart against the continuous run", flush=True)
    h = f"{work}/H"
    hc.make_hopper(h, d=1e-4, end=0.25, write_interval=0.05)
    pc.run(h, set_fields=True)
    restart_copy(h, f"{h}r", 0.2)
    d = compare(h, f"{h}r", "0.25")
    worst = max(d, key=d.get)
    report("hopper restarted at 0.2 s (outlet open)", d[worst] <= 1e-12 and "p" in d,
           f"{len(d)} fields, largest difference {d[worst]:.1e} ({worst})")

    w = f"{work}/D"
    cc.make(w, end=0.14)
    cc.run(w)
    restart_copy(w, f"{w}r", 0.1)
    d = compare(w, f"{w}r", "0.14")
    worst = max(d, key=d.get)
    report("dam break restarted at 0.1 s", d[worst] <= 1e-12 and "p" in d,
           f"{len(d)} fields, largest difference {d[worst]:.1e} ({worst})")
    log = open(f"{w}r/log.celFoam").read()
    report("the restart reads the gas state", "restart state read from 0.1/uniform/poreGas" in log,
           "0.1/uniform/poreGas")


def check_processors(work):
    print("B: 2 and 4 processors against serial", flush=True)
    secs = {}
    for n, dc in ((1, (1, 1, 1)), (2, (2, 1, 1)), (4, (2, 2, 1))):
        p = f"{work}/HP{n}"
        hc.make_hopper(p, d=1e-4, end=0.25, write_interval=0.05, decomp=dc, rules=False,
                       extra_gas="pTolerance 1e-12;\n")
        secs[n] = pc.run(p, n, set_fields=True)
    nc = hc.ncells()
    for n in (2, 4):
        a, b = f"{work}/HP1", f"{work}/HP{n}"
        worst = max(max(rel(pc.field(a, t, f, nc), pc.field(b, t, f, nc)) for f in ("p", "U"))
                    for t in pc.times(a))
        m, ni, _, _ = pc.mass_log(b)
        drift = float(np.abs((m - ni) - (m[0] - ni[0])).max() / m[0])
        report(f"hopper on {n} processors", worst <= 1e-10 and drift <= 1e-12,
               f"p, U within {worst:.1e}, gas mass {drift:.1e}")

    end, dt_w = 0.02, 0.002
    for n in (1, 2, 4):
        p = f"{work}/DP{n}"
        cc.make(p, end=end)
        cd = f"{p}/system/controlDict"
        txt = open(cd).read().replace(f"writeInterval   {cc.WRITE!r};", f"writeInterval   {dt_w!r};")
        open(cd, "w").write(txt)
        if n > 1:
            ck.write_manual(p, n)
        cc.run(p, n)
    for n in (2, 4):
        a, b = f"{work}/DP1", f"{work}/DP{n}"
        last, first = 0.0, None
        for t in pc.times(a):
            d = max(compare(a, b, t, ("p", "U", "sigma", "alpha")).values())
            if d > 1e-10:
                first = (float(t), d)
                break
            last = float(t)
        # reported, not a pass condition: the done-when's 1e-10 over the whole
        # run is not met for the dam break (decision 0141; Sia to decide)
        print(f"  REPORT dam break on {n} processors (cuts through the column): within 1e-10 up "
              f"to {last:g} s"
              + (f", then {first[1]:.1e} at {first[0]:g} s (round-off flips a free-surface rule)"
                 if first else f" (to {end:g} s)"), flush=True)


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check_restart(work)
    check_processors(work)
    print(f"gas restart and MPI checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_gas_restart_mpi: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
