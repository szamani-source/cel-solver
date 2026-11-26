"""WP T.8 done-when: temperature restarts and runs in parallel (decision 0136).

Cases, both from earlier packages:
S  the T.4 Schumann column (gas, twoTemperature: T and Tg), to 1 s;
M  the T.6 moving bed past a heated wall (equilibrium, inlet, outlet, flux
   wall), to 0.7 s.

A  Restart: each case run again from a copy holding only a middle write time
   (startFrom latestTime; ascii at 17 digits round-trips every double) matches
   the continuous run at the end to 1e-10 of its temperature range (S: T, Tg
   restarted at 0.5 s; M: T restarted at 0.3 s).
B  Processors: S on 4 (cuts along the column) and M on 2 and 4 (cuts along the
   flow) match serial to 1e-10 of the range.
C  Cost: the S column with temperature against the same column without
   thermalProperties (the gas solved in both): wall time per step at most
   1.2 times (best of two runs each).

    python3 thermal_restart_mpi_checks.py WORKDIR
"""

import os
import re
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import conduction_checks as cc  # noqa: E402
import gas_heat_checks as gh  # noqa: E402
import moving_bed_checks as mb  # noqa: E402

report = cc.report
failed = cc.failed


def restart_copy(src, dst, t_restart):
    shutil.copytree(src, dst)
    for d in os.listdir(dst):
        if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > t_restart + 1e-12:
            shutil.rmtree(os.path.join(dst, d))
    for d in os.listdir(dst):
        if d.startswith("processor"):
            shutil.rmtree(os.path.join(dst, d))
    cd = f"{dst}/system/controlDict"
    txt = open(cd).read()
    assert "startFrom       startTime;" in txt
    open(cd, "w").write(txt.replace("startFrom       startTime;", "startFrom       latestTime;"))
    cc.sh(dst, ["celFoam"], "celFoam")


def diff(a, b, t, fields, scale):
    return max(np.abs(cc.read_internal(f"{a}/{t:g}/{f}") - cc.read_internal(f"{b}/{t:g}/{f}")).max()
               for f in fields) / scale


def steps_and_time(path):
    log = open(f"{path}/log.celFoam").read()
    n = len(re.findall(r"^Time = ", log, re.M))
    secs = float(re.findall(r"^ExecutionTime = (\S+) s", log, re.M)[-1])
    return n, secs


def check_schumann(work):
    print("[S] Schumann column (T.4)", flush=True)
    s1 = os.path.join(work, "S1")
    t0 = time.perf_counter()
    gh.column(s1, "twoTemperature", gh.T0 + gh.DT0, end=1.0, write=0.5)
    print(f"    serial {time.perf_counter() - t0:.0f} s", flush=True)
    sr = os.path.join(work, "S1r")
    restart_copy(s1, sr, 0.5)
    d = diff(s1, sr, 1.0, ("T", "Tg"), gh.DT0)
    report("Schumann column restarted at 0.5 s", d < 1e-10, f"largest difference {d:.1e} of the step")
    s4 = os.path.join(work, "S4")
    gh.column(s4, "twoTemperature", gh.T0 + gh.DT0, decomp=(4, 1, 1), end=1.0, write=0.5)
    d = diff(s1, s4, 1.0, ("T", "Tg"), gh.DT0)
    report("Schumann column on 4 processors", d < 1e-10, f"largest difference {d:.1e} of the step")
    return s1


def check_bed(work):
    print("[M] Moving bed past a heated wall (T.6)", flush=True)
    m1 = os.path.join(work, "M1")
    mb.make(m1)
    cc.run(m1)
    rng = np.ptp(cc.read_internal(f"{m1}/{mb.END:g}/T"))
    mr = os.path.join(work, "M1r")
    restart_copy(m1, mr, 0.3)
    d = diff(m1, mr, mb.END, ("T",), rng)
    report("moving bed restarted at 0.3 s", d < 1e-10, f"largest difference {d:.1e} of the range")
    for n in (2, 4):
        mn = os.path.join(work, f"M{n}")
        mb.make(mn, decomp=(n, 1, 1))
        cc.run(mn, n)
        d = diff(m1, mn, mb.END, ("T",), rng)
        report(f"moving bed on {n} processors", d < 1e-10, f"largest difference {d:.1e} of the range")


def check_cost(work):
    print("[C] Cost of temperature with the gas", flush=True)
    # Best of two runs each (wall time on a shared machine varies by 10 %)
    runs = {}
    for model in ("twoTemperature", "none"):
        best = None
        for k in range(2):
            p = os.path.join(work, f"cost_{model}_{k}")
            gh.column(p, model, gh.T0 + gh.DT0, end=1.0, write=0.5)
            n, t = steps_and_time(p)
            best = (n, t) if best is None or t < best[1] else best
        runs[model] = best
    n1, t1 = runs["twoTemperature"]
    n0, t0 = runs["none"]
    r = (t1 / n1) / (t0 / n0)
    report("temperature costs at most 1.2 times a gas-on step", r <= 1.2,
           f"{t1 / n1 * 1e3:.3f} ms against {t0 / n0 * 1e3:.3f} ms per step ({r:.3f} times), "
           f"{n1} and {n0} steps")


def main(work):
    os.makedirs(work, exist_ok=True)
    t0 = time.perf_counter()
    check_schumann(work)
    check_bed(work)
    check_cost(work)
    print(f"thermal restart and MPI checks took {time.perf_counter() - t0:.0f} s")
    if failed:
        print("FAILED: " + ", ".join(failed))
        sys.exit(1)
    print("test_thermal_restart_mpi: PASS")


if __name__ == "__main__":
    main(sys.argv[1])
