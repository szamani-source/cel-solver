"""Checkpoint and restart of the test driver celVoidCheck (decision 0116).

The Couto silo on the coarse wall-aligned STL mesh (buffered_case,
gravity-on start, nonlocal fluidity, fullBlend, outlet opening at 3 ms) on
2 ranks, run to 6 ms with checkpoints every 2 ms: (A) straight through,
(B) to 5 ms, then restarted from the latest checkpoint (4 ms) to 6 ms.
Pass: the cells at 6 ms and the history agree to 1e-9 (relative to the
largest value of each column).

    python3 voidcheck_restart.py WORKDIR
"""
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "hopper"))
import hopper_pressure_checks as hp  # noqa: E402
import buffered_case as bc  # noqa: E402

vc = hp.vc
NP = 2
os.environ["CEL_HOPPER_RULES"] = hp.MU_I + hp.RULES + "fullBlend yes;\nnonlocalA 0.48;\n"


def case(work, name, t_end):
    check = dict(tEnd=t_end, maxCo=0.4, writeTimes=[0.006], tOpen=0.003,
                 checkpointEvery=0.002)
    p, c = bc.aligned_buffered_case(work, name, "mc", check, nProcs=NP)
    dp = f"{p}/system/decomposeParDict"
    txt = open(dp).read().replace(f"n ({NP} 1 1)", f"n (1 {NP} 1)")
    open(dp, "w").write(txt)
    open(f"{p}/system/controlDict", "a").write(
        "gravityRamp 0.1;\ndamping 100;\nwriteFormat binary;\n")
    return p


def resume(p, t_end):
    """Restart p from its latest written time and run to t_end."""
    d = f"{p}/system/voidCheckDict"
    txt = re.sub(r"tEnd [^;]*;", f"tEnd {t_end!r};", open(d).read())
    open(d, "w").write(txt)
    open(f"{p}/system/controlDict", "a").write("startFrom latestTime;\n")
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    with open(f"{p}/log.restart", "w") as f:
        subprocess.run(["mpirun", "-np", str(NP), "--oversubscribe", "celVoidCheck",
                        "-parallel"], cwd=p, stdout=f, stderr=subprocess.STDOUT,
                       env=env, check=True)


def main():
    work = os.path.abspath(sys.argv[1])
    shutil.rmtree(work, ignore_errors=True)
    a = case(work, "straight", 0.006)
    vc.run(a, NP)
    b = case(work, "restarted", 0.005)
    vc.run(b, NP)
    last = max(float(t) for t in os.listdir(f"{b}/processor0")
               if re.fullmatch(r"[0-9.e-]+", t) and float(t) > 0)
    resume(b, 0.006)
    ok = "Restart from t = " in open(f"{b}/log.restart").read()
    ca, cb = hp.read_cells(a, "cells_0.006"), hp.read_cells(b, "cells_0.006")
    ha, hb = vc.history(a), vc.history(b)
    ok &= ha.shape == hb.shape
    cell_err = np.max(np.abs(ca - cb) / np.maximum(np.abs(ca).max(0), 1e-300))
    hist_err = (np.max(np.abs(ha - hb)[:, 2:12] / np.maximum(np.abs(ha[:, 2:12]).max(0), 1e-300))
                if ok else np.inf)
    print(f"restart from t = {last:g} s; history {ha.shape[0]} / {hb.shape[0]} rows; "
          f"largest relative difference: cells {cell_err:.2e}, history {hist_err:.2e}")
    ok &= cell_err < 1e-9 and hist_err < 1e-9
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
