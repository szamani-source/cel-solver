"""WP G.6 done-when checks (first part): the optional faster-time-step mesh
options are off by default and lengthen the step when on.

1. Defaults: the Couto mesh and its time step are unchanged (5.85003e-6 s).
2. conicalHopper outletColumns even + axisColumn 2: same cell count, the step
   at least 1.5 times longer (celCheck), Mesh OK, the same mass placed, and
   celFoam runs 50 steps keeping the mass.
3. flatBottom (axisymmetric) axisColumn 2: the step at least 1.3 times longer.
4. Asking Gmsh for these options is an error that says why.

Whether wall pressures and the discharge stay the same with them is not
checked here (it needs a settling run); they stay off until it is.

    python3 tests/celsetup/fastdt_checks.py WORKDIR
"""
import copy
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from celsetup import case, miniyaml  # noqa: E402

EX = os.path.join(ROOT, "tools", "celsetup", "examples")
fails = []


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg, flush=True)
    if not ok:
        fails.append(msg)


def make(raw, d):
    case.write_case(case.resolve(raw), d)
    ok, msg = case.run_mesh_and_fill(d)
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "celCheck"), d, "--json",
                    d + ".json"], capture_output=True)
    return ok, case.fill(d), json.load(open(d + ".json"))


def main(work):
    os.makedirs(work, exist_ok=True)
    couto = miniyaml.load(os.path.join(EX, "coutoHopper.yaml"))
    ok, f0, r0 = make(couto, os.path.join(work, "default"))
    check(ok and abs(r0["dt"] / 5.85003e-6 - 1) < 1e-4, f"defaults: Couto step {r0['dt']:.6g} s")
    fast = copy.deepcopy(couto)
    fast["mesh"].update(outletColumns="even", axisColumn=2)
    d = os.path.join(work, "fast")
    ok, f1, r1 = make(fast, d)
    check(ok and r1["nCells"] == r0["nCells"] and r1["dt"] / r0["dt"] > 1.5
          and abs(f1["mass"] / f0["mass"] - 1) < 2e-3,
          f"outletColumns even, axisColumn 2: step {r1['dt']:.4g} s, "
          f"{r1['dt'] / r0['dt']:.2f} times longer; {f1['mass']:.1f} kg placed")
    with open(os.path.join(d, "log.checkMesh"), "w") as f:
        subprocess.run(["checkMesh"], cwd=d, stdout=f, stderr=subprocess.STDOUT,
                       env=dict(os.environ, PWD=d))
    check("Mesh OK." in open(os.path.join(d, "log.checkMesh")).read(), "checkMesh: Mesh OK")
    for k in ("endTime", "writeInterval"):
        subprocess.run(["foamDictionary", "-entry", k, "-set", repr(50.5 * r1["dt"]),
                        "system/controlDict"], cwd=d, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    txt = re.sub(r"^NP=\d+", "NP=1", open(f"{d}/Allrun").read(), flags=re.M)
    open(f"{d}/Allrun", "w").write(txt)
    rc = subprocess.run(["./Allrun"], cwd=d, env=dict(os.environ, PWD=d)).returncode
    log = open(os.path.join(d, "log.celFoam")).read()
    mass = [float(x) for x in re.findall(r"Total mass\s+=\s+(\S+)", log)]
    sd = re.search(r"Stable time step\s+deltaT = (\S+)", log)
    steps = len(re.findall(r"^Time = ", log, flags=re.M))
    check(rc == 0 and steps >= 45 and len(mass) >= 2 and abs(mass[-1] / mass[0] - 1) < 1e-9
          and sd and abs(float(sd.group(1)) / r1["dt"] - 1) < 1e-3,
          f"celFoam runs {steps} steps at {sd.group(1) if sd else '?'} s, mass kept")

    silo = miniyaml.load(os.path.join(EX, "planeTank.yaml"))
    silo["geometry"]["shape"] = "axisymmetric"
    _, _, s0 = make(silo, os.path.join(work, "silo"))
    silo["mesh"]["axisColumn"] = 2
    _, _, s1 = make(silo, os.path.join(work, "siloFast"))
    check(s1["dt"] / s0["dt"] > 1.3, f"flatBottom axisColumn 2: step {s1['dt']:.4g} s, "
          f"{s1['dt'] / s0['dt']:.2f} times longer")
    g = copy.deepcopy(fast)
    g["mesh"]["mesher"] = "gmsh"
    try:
        case.write_case(case.resolve(g), os.path.join(work, "g"))
        err = ""
    except ValueError as e:
        err = str(e)
    check("need mesher: blockMesh" in err, f"Gmsh with these options: {err}")
    print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
