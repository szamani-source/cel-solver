"""WP 6.3a: parallel scaling of celFoam on the T0 tank (tutorials/siloDischarge,
E 10 MPa, T0 rule set). Run by tests/manual/test_scaling.sh:

    python3 scaling_run.py WORKDIR [--cells 10,5] [--procs 1,2,4]

For each cell size it runs the settling tank (gravity ramp, damping, outlet
closed) for a fixed span of time on 1, 2 and 4 processors (simple split
along x), and prints wall time, steps, cost per cell per step and speed-up
against 1 processor. Checks only that every run finished with the same number
of steps and the same mass (the runs are short, the numbers are the result).
"""

import argparse
import os
import re
import shutil
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SPAN = {10: 0.05, 5: 0.02}          # simulated time [s]


def sh(path, cmd, log):
    env = dict(os.environ, PWD=os.path.abspath(path), OMPI_ALLOW_RUN_AS_ROOT="1",
               OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    with open(f"{path}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed in {path}")


def sub(fn, pat, rep):
    s = open(fn).read()
    s, n = re.subn(pat, rep, s, flags=re.M)
    assert n > 0, (fn, pat)
    open(fn, "w").write(s)


def case(work, mm, n):
    path = f"{work}/tank{mm}mm_np{n}"
    shutil.copytree(f"{ROOT}/tutorials/siloDischarge", path,
                    ignore=shutil.ignore_patterns("*.png", "*.py", "README.md"))
    f = 10 // mm
    bm = f"{path}/system/blockMeshDict"
    sub(bm, r"\(17 62 1\)", f"({17 * f} {62 * f} 1)")
    sub(bm, r"\(6 62 1\)", f"({6 * f} {62 * f} 1)")
    if mm != 10:
        # one cell thick: the slab is as thick as the cells (stableDeltaT uses cbrt V)
        sub(bm, r" 0\.01\)", f" {mm / 1000!r})")
    cd = f"{path}/system/controlDict"
    sub(cd, r"^endTime\s+[^;]*;", f"endTime {SPAN[mm]!r};")
    sub(cd, r"^writeInterval\s+[^;]*;", "writeInterval 100;")
    dp = f"{path}/system/decomposeParDict"
    sub(dp, r"^numberOfSubdomains\s+[^;]*;", f"numberOfSubdomains {n};")
    sub(dp, r"n \(2 1 1\)", f"n ({n} 1 1)")
    sh(path, ["blockMesh"], "blockMesh")
    shutil.copytree(f"{path}/0.orig", f"{path}/0")
    sh(path, ["setFields"], "setFields")
    if n > 1:
        sh(path, ["decomposePar", "-force"], "decomposePar")
    t0 = time.perf_counter()
    if n > 1:
        sh(path, ["mpirun", "-np", str(n), "--oversubscribe", "celFoam", "-parallel"], "celFoam")
    else:
        sh(path, ["celFoam"], "celFoam")
    wall = time.perf_counter() - t0
    log = open(f"{path}/log.celFoam").read()
    steps = len(re.findall(r"^Time = ", log, flags=re.M))
    mass = re.findall(r"^Total mass\s*=\s*(\S+)", log, flags=re.M)[-1]
    ncell = 40 * 62 * f * f
    return dict(wall=wall, steps=steps, mass=mass, cells=ncell)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("work")
    ap.add_argument("--cells", default="10,5")
    ap.add_argument("--procs", default="1,2,4")
    a = ap.parse_args()
    ok = True
    print(f"{'cells':>6} {'procs':>5} {'cells':>7} {'steps':>6} {'wall s':>8} "
          f"{'us/cell/step':>12} {'speed-up':>8} {'efficiency':>10}")
    for mm in (int(x) for x in a.cells.split(",")):
        base = None
        for n in (int(x) for x in a.procs.split(",")):
            r = case(a.work, mm, n)
            if base is None:
                base = r
            same = r["steps"] == base["steps"] and r["mass"] == base["mass"]
            ok = ok and r["steps"] == base["steps"]
            cost = r["wall"] * n / (r["cells"] * r["steps"]) * 1e6
            sp = base["wall"] / r["wall"]
            print(f"{mm:>4}mm {n:>5} {r['cells']:>7} {r['steps']:>6} {r['wall']:>8.1f} "
                  f"{cost:>12.2f} {sp:>8.2f} {100 * sp / n:>9.0f}%"
                  + ("" if same else f"  (mass {r['mass']} vs {base['mass']})"))
    print(f"\nWP 6.3a scaling: {'DONE' if ok else 'STEP COUNTS DIFFER'}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
