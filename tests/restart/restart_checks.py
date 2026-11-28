"""WP 6.2a done-when check: a restarted celFoam run matches the continuous
run (decision 0081). Run by tests/test_restart.sh after celFoam is built:

    python3 restart_checks.py WORKDIR [case ...]

Each case runs celFoam from 0 to T_END, writing at T_RESTART and T_END
(writeFormat binary, so fields round-trip exactly), then again from a copy
that holds only the T_RESTART time (startFrom latestTime) to T_END. The
fields at T_END (rho, U, sigma, epsP, epsPEq, alpha), the body state and the
mass counters must be identical (bit for bit), and so must the appended
lines of postProcessing/bodyForces.dat after the restart time.

Cases (the T0 tank of tutorials/solid/siloDischarge, E 0.1 MPa, 10 mm cells):
  silo     outlet opens at 0.3 s, restart at 0.4 s while it discharges, end
           0.5 s (flow impedance, outlet, damping switched off at opening)
  bodies   tank filled to 0.3 m, outlet closed; a free steel block (rigid,
           STL) drops onto the surface on the left, a steel footing of hex8
           elements is pushed in on the right; restart at 0.05 s (both in
           contact), end 0.1 s
  parallel the bodies case on 2 processors (split between the bodies), run
           and restarted in parallel

Prints PASS/FAIL per case, exits non-zero on any failure.
"""

import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
for d in ("lagrange", "stlwalls"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import stl_wall_checks as sw  # noqa: E402

FIELDS = ["rho", "U", "sigma", "epsP", "epsPEq", "alpha"]


def sh(path, cmd, log):
    env = dict(os.environ, PWD=os.path.abspath(path), OMPI_ALLOW_RUN_AS_ROOT="1",
               OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    with open(f"{path}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=path, stdout=f, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed in {path}:\n"
                           + open(f"{path}/log.{log}").read()[-3000:])


def set_entry(path, rel, key, value):
    fn = f"{path}/{rel}"
    s = open(fn).read()
    s, n = re.subn(rf"^(\s*){key}\s+[^;]*;", rf"\g<1>{key} {value};", s, flags=re.M)
    if n == 0:
        s += f"{key} {value};\n"
    open(fn, "w").write(s)


def make_silo(path, t_restart, t_end, fill=0.6, open_time=0.3):
    shutil.copytree(f"{ROOT}/tutorials/solid/siloDischarge", path,
                    ignore=shutil.ignore_patterns("*.png", "*.py", "README.md"))
    set_entry(path, "constant/materialProperties", "E", "1e5")
    s = open(f"{path}/constant/materialProperties").read()
    s = re.sub(r"openTime\s+[^;]*;", f"openTime {open_time!r};", s)
    open(f"{path}/constant/materialProperties", "w").write(s)
    s = open(f"{path}/system/setFieldsDict").read()
    s = s.replace("box (-1 -1 -1) (1 0.6 1);", f"box (-1 -1 -1) (1 {fill!r} 1);")
    open(f"{path}/system/setFieldsDict", "w").write(s)
    cd = "system/controlDict"
    set_entry(path, cd, "endTime", repr(t_end))
    set_entry(path, cd, "writeInterval", repr(round(t_end - t_restart, 12)))
    set_entry(path, cd, "writeFormat", "binary")
    set_entry(path, cd, "writePrecision", "17")
    set_entry(path, cd, "timePrecision", "12")
    return path


def add_bodies(path):
    """Free steel block on the left, hex8 steel footing on the right."""
    lo, hi = (0.05, 0.302, -0.01), (0.12, 0.372, 0.02)
    cs.write(f"{path}/constant/triSurface/block.stl", sw.box_stl(lo, hi))
    m = 7850.0 * (hi[0] - lo[0]) * (hi[1] - lo[1]) * 0.01
    v = 0.25
    table = f"table ((0 (0 {-v!r} 0)) (10 (0 {-v!r} 0)))"
    flo, fhi = (0.27, 0.305, -0.01), (0.34, 0.375, 0.02)
    bodies = (
        f"    block\n    {{\n        surface block.stl;\n        friction 20;\n"
        f"        mass {m!r};\n        velocity (0 0 0);\n        omega (0 0 0);\n    }}\n"
        f"    footing\n    {{\n        deformable\n        {{\n"
        f"            min ({flo[0]!r} {flo[1]!r} {flo[2]!r});\n"
        f"            max ({fhi[0]!r} {fhi[1]!r} {fhi[2]!r});\n"
        f"            cells (7 7 1);\n            density 7850;\n"
        f"            material {{ type linearElastic; E 210e9; nu 0.3; }}\n"
        f"            planeStrain yes;\n            damping 0.2;\n"
        f"            drive ({flo[0] - 0.01!r} {fhi[1] - 1e-6!r} -1) "
        f"({fhi[0] + 0.01!r} {fhi[1] + 1e-6!r} 1);\n        }}\n"
        f"        velocity {table};\n        friction 20;\n    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + bodies + "}\n")


def run_celfoam(path, nProcs, log):
    if nProcs > 1:
        sh(path, ["mpirun", "-np", str(nProcs), "--oversubscribe", "celFoam", "-parallel"], log)
    else:
        sh(path, ["celFoam"], log)


def time_dirs(path, nProcs):
    base = f"{path}/processor0" if nProcs > 1 else path
    out = []
    for d in os.listdir(base):
        try:
            out.append((float(d), d))
        except ValueError:
            pass
    return sorted(out)


def compare_dirs(a, b, files):
    """Names of files that differ (or are missing) between two directories."""
    bad = []
    for f in files:
        fa, fb = f"{a}/{f}", f"{b}/{f}"
        if not os.path.isfile(fa) and not os.path.isfile(fb):
            continue
        if not (os.path.isfile(fa) and os.path.isfile(fb)):
            bad.append(f + " (missing)")
            continue
        if open(fa, "rb").read() != open(fb, "rb").read():
            bad.append(f)
    return bad


def log_value(path, log, key):
    vals = re.findall(rf"^{key}\s*=\s*(\S+)", open(f"{path}/log.{log}").read(), flags=re.M)
    return vals[-1] if vals else None


def forces_after(path, t0):
    fn = f"{path}/postProcessing/bodyForces.dat"
    if not os.path.isfile(fn):
        return []
    rows = [line for line in open(fn) if not line.startswith("#")]
    return [r for r in rows if float(r.split()[0]) > t0 + 1e-12]


def check(work, name, build, t_restart, t_end, nProcs=1):
    print(f"[{name}] continuous run 0 to {t_end} s against a restart from {t_restart} s"
          + (f" ({nProcs} processors)" if nProcs > 1 else ""))
    t0 = time.perf_counter()
    full = build(f"{work}/{name}_full", t_restart, t_end)
    sh(full, ["blockMesh"], "blockMesh")
    shutil.copytree(f"{full}/0.orig", f"{full}/0")
    sh(full, ["setFields"], "setFields")
    if nProcs > 1:
        set_entry(full, "system/decomposeParDict", "numberOfSubdomains", str(nProcs))
        sh(full, ["decomposePar", "-force"], "decomposePar")
    run_celfoam(full, nProcs, "celFoam")

    times = time_dirs(full, nProcs)
    names = {round(t, 9): d for t, d in times}
    tr, te = names.get(round(t_restart, 9)), names.get(round(t_end, 9))
    if tr is None or te is None:
        print(f"    FAIL: written times {[d for _, d in times]}, need {t_restart} and {t_end}")
        return False

    # The restart copy: everything but the times after the restart time
    rest = f"{work}/{name}_restart"
    shutil.copytree(full, rest, ignore=shutil.ignore_patterns(te, "log.*"))
    # bodyForces.dat as the first run left it at the restart time
    fn = f"{rest}/postProcessing/bodyForces.dat"
    if os.path.isfile(fn):
        keep = [r for r in open(fn) if r.startswith("#") or float(r.split()[0]) <= t_restart + 1e-12]
        open(fn, "w").writelines(keep)
    set_entry(rest, "system/controlDict", "startFrom", "latestTime")
    run_celfoam(rest, nProcs, "celFoam")

    files = FIELDS + ["uniform/time", "uniform/rigidBodies", "uniform/eulerianSolid"]
    bad = []
    procs = [f"processor{p}/" for p in range(nProcs)] if nProcs > 1 else [""]
    for pr in procs:
        bad += [pr + f for f in compare_dirs(f"{full}/{pr}{te}", f"{rest}/{pr}{te}", files)]
    ok = not bad
    det = [f"fields and state at {te} s: " + ("identical" if not bad else "differ: " + ", ".join(bad))]
    for key in ("Total mass", "Discharged mass"):
        a, b = log_value(full, "celFoam", key), log_value(rest, "celFoam", key)
        if a is None:
            continue
        ok = ok and a == b
        det.append(f"{key.lower()} at the end: continuous {a}, restarted {b}")
    fa, fr = forces_after(full, t_restart), forces_after(rest, t_restart)
    if fa or fr:
        same = fa == fr
        ok = ok and same and len(fa) > 0
        det.append(f"bodyForces.dat after {t_restart} s: {len(fa)} lines, "
                   + ("identical" if same else "differ"))
    det.append(f"{time.perf_counter() - t0:.0f} s")
    print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    for line in det:
        print(f"    {line}")
    return ok


def build_bodies(path, t_restart, t_end):
    make_silo(path, t_restart, t_end, fill=0.3, open_time=100.0)
    add_bodies(path)
    return path


def main():
    work = sys.argv[1]
    cases = {
        "silo": lambda: check(work, "silo", make_silo, 0.4, 0.5),
        "bodies": lambda: check(work, "bodies", build_bodies, 0.05, 0.1),
        "parallel": lambda: check(work, "parallel", build_bodies, 0.05, 0.1, nProcs=2),
    }
    only = sys.argv[2:] or list(cases)
    t0 = time.perf_counter()
    ok = True
    for name in only:
        ok = cases[name]() and ok
    print(f"\nWP 6.2a checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
