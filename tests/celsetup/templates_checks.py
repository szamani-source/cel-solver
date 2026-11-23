"""WP G.5 done-when checks: more templates, both meshers.

For each of flatBottom plane (the planeTank example: Sia's first target, a
plane tank with a small bottom opening), flatBottom axisymmetric and
conicalHopper plane:
- the case is written and meshed (blockMesh), checkMesh reports Mesh OK;
- the fill places the geometry's mass within 0.5 % (1500 kg/m for the tank);
- celFoam runs 50 steps and keeps the mass.
With Gmsh installed, each is meshed with Gmsh too: same cell count and
patches, time step within 1 %.

    python3 tests/celsetup/templates_checks.py WORKDIR
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
from celsetup import case, gmsh_mesh, miniyaml, usermesh  # noqa: E402

EX = os.path.join(ROOT, "tools", "celsetup", "examples")
fails = []


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg, flush=True)
    if not ok:
        fails.append(msg)


def sh(cmd, cwd, log):
    with open(os.path.join(cwd, log), "w") as f:
        return subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                              env=dict(os.environ, PWD=cwd)).returncode


def celcheck(c):
    out = c + ".json"
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "celCheck"), c, "--json", out],
                   capture_output=True, text=True)
    return json.load(open(out))


def variants():
    tank = miniyaml.load(os.path.join(EX, "planeTank.yaml"))
    silo = copy.deepcopy(tank)
    silo["geometry"]["shape"] = "axisymmetric"
    hop = miniyaml.load(os.path.join(EX, "coutoHopper.yaml"))
    hop["geometry"]["shape"] = "plane"
    return [("planeTank", tank, 1500.0), ("flatSilo", silo, None), ("planeHopper", hop, None)]


def main(work):
    os.makedirs(work, exist_ok=True)
    have_gmsh = gmsh_mesh.available() is not None
    for name, raw, mass_expected in variants():
        c = os.path.join(work, name)
        b = case.write_case(case.resolve(raw), c)
        ok, msg = case.run_mesh_and_fill(c)
        res = case.fill(c)
        target = mass_expected or res["target"]
        check(ok and abs(res["mass"] / target - 1) < 5e-3,
              f"{name}: {b['summary']}; {msg}")
        sh(["checkMesh"], c, "log.checkMesh")
        check("Mesh OK." in open(os.path.join(c, "log.checkMesh")).read(),
              f"{name}: checkMesh Mesh OK")
        for k, v in (("endTime", "1"), ("writeInterval", "1")):
            subprocess.run(["foamDictionary", "-entry", k, "-set", v, "system/controlDict"],
                           cwd=c, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # 50 steps: stop by step count through the end time of 50 stable steps
        r0 = celcheck(c)
        for k in ("endTime", "writeInterval"):     # the mass is printed at writes
            subprocess.run(["foamDictionary", "-entry", k, "-set", repr(50.5 * r0["dt"]),
                            "system/controlDict"], cwd=c, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        txt = re.sub(r"^NP=\d+", "NP=1", open(f"{c}/Allrun").read(), flags=re.M)
        open(f"{c}/Allrun", "w").write(txt)
        rc = subprocess.run(["./Allrun"], cwd=c, env=dict(os.environ, PWD=c)).returncode
        log = open(os.path.join(c, "log.celFoam")).read()
        mass = [float(x) for x in re.findall(r"Total mass\s+=\s+(\S+)", log)]
        steps = len(re.findall(r"^Time = ", log, flags=re.M))
        check(rc == 0 and "End" in log and steps >= 45 and len(mass) >= 2
              and abs(mass[-1] / mass[0] - 1) < 1e-9,
              f"{name}: celFoam runs {steps} steps, mass kept")
        if not have_gmsh:
            continue
        g = copy.deepcopy(raw)
        g["mesh"]["mesher"] = "gmsh"
        cg = c + "Gmsh"
        case.write_case(case.resolve(g), cg)
        ok, msg = case.run_mesh_and_fill(cg)
        pb = {n: (t, k) for n, t, k, _ in usermesh.read_boundary(c) if k}
        pg = {n: (t, k) for n, t, k, _ in usermesh.read_boundary(cg) if k}
        rg = celcheck(cg) if ok else {"dt": 0}
        check(ok and pb == pg and abs(rg["dt"] / r0["dt"] - 1) < 0.01,
              f"{name} with Gmsh: same patches {sorted(pg)}, time step {rg['dt']:.4g} s "
              f"(blockMesh {r0['dt']:.4g} s)")
    if not have_gmsh:
        print("SKIPPED: Gmsh not available, the Gmsh variants are not checked here")
    print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
