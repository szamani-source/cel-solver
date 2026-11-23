"""WP S.1 done-when checks: celSetup writes a complete, runnable Couto case.

1. The YAML reader reads the example (numbers, yes/no, comments) and a
   missing required entry gives a readable error.
2. With the Couto inputs and the validated mesh numbers, the mesh points
   and the STL shell equal those of tests/hopper/buffered_case.py
   (aligned_buffered_case, decision 0107) exactly.
3. checkMesh reports Mesh OK.
4. The fill places the Couto mass (1503 kg) within 0.2 %.
5. celFoam runs 0.3 ms (about 50 steps) from the written case, keeps the
   mass to 1e-9, and a second Allrun resumes from the last write.

    python3 tests/celsetup/celsetup_checks.py WORKDIR
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
sys.path.insert(0, os.path.join(ROOT, "tests", "hopper"))
from celsetup import case, foam, miniyaml  # noqa: E402

fails = []


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg, flush=True)
    if not ok:
        fails.append(msg)


def sh(cmd, cwd, log):
    with open(os.path.join(cwd, log), "w") as f:
        return subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                              env=dict(os.environ, PWD=cwd)).returncode


def main(work):
    os.makedirs(work, exist_ok=True)
    ex = os.path.join(ROOT, "tools", "celsetup", "examples", "coutoHopper.yaml")
    raw = miniyaml.load(ex)
    check(raw["material"]["E"] == 10.674e6 and raw["geometry"]["fillHeight"] == "top"
          and raw["walls"]["frictionCoefficient"] == 0.2, "YAML reader: numbers, words, comments")
    bad = dict(raw, material=dict(raw["material"]))
    del bad["material"]["d"]
    try:
        case.resolve(bad)
        check(False, "missing material d is reported")
    except ValueError as e:
        check("'d' is required" in str(e), f"missing material d is reported: {e}")

    # S.2: the template file lists every entry; filled in, it resolves
    from celsetup import schema, conical
    tpl = schema.template_yaml("conicalHopper", conical.SCHEMA)
    allkeys = [e.key for ents in schema.sections(conical.SCHEMA).values() for e in ents]
    missing = [k for k in allkeys if not any(l.lstrip(" #").startswith(k + ":")
                                            for l in tpl.splitlines())]
    check(not missing, f"--template lists all {len(allkeys)} entries with descriptions"
          + (f" (missing {missing})" if missing else ""))
    t = miniyaml.loads(tpl)
    t["geometry"].update(raw["geometry"])
    t["walls"]["frictionCoefficient"] = 0.2
    for k in ("rho", "E", "phi", "d"):
        t["material"][k] = raw["material"][k]
    cfg = case.resolve(t)
    check(cfg["mesh"]["outletCell"] == 0.004 and cfg["material"]["fullBlend"] == "auto"
          and abs(cfg["walls"]["friction"] - 11.3099) < 1e-3,
          "filled-in template resolves (defaults, preset values, friction from mu_w)")
    try:
        case.resolve(dict(raw, mesh=dict(raw["mesh"], outletcell=0.003)))
        check(False, "a misspelt entry is reported")
    except ValueError as e:
        check("outletcell" in str(e), f"a misspelt entry is reported: {e}")
    t2 = dict(t, material=dict(t["material"], fullBlend="no", model="druckerPrager"))
    mt = case.material_text(case.resolve(t2)["material"])
    check("fullBlend       no;" in mt and "shearCap        yes;" in mt
          and "fit             plane_strain;" in mt, "a single switch can be turned off; DP fit")

    # 2. same mesh and STL as the validated script
    p = os.path.join(work, "couto")
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "celSetup"), ex, p],
                       capture_output=True, text=True)
    print(r.stdout + r.stderr)
    check(r.returncode == 0, "celSetup writes, meshes and fills the Couto case")
    import buffered_case as bc
    ref, _ = bc.aligned_buffered_case(work, "reference", "mc",
                                      dict(tEnd=1e-6, maxCo=0.4, writeTimes=[]), {}, 1)
    sh(["blockMesh"], ref, "log.blockMesh")
    a = foam.read_points(f"{ref}/constant/polyMesh/points")
    b = foam.read_points(f"{p}/constant/polyMesh/points")
    same = a.shape == b.shape and float(np.abs(a - b).max()) == 0.0
    check(same, f"mesh points equal the validated mesh ({b.shape[0]} points)")
    check(open(f"{ref}/constant/triSurface/shell.stl").read()
          == open(f"{p}/constant/triSurface/shell.stl").read(), "STL shell equals the validated one")

    # 3. checkMesh
    sh(["checkMesh"], p, "log.checkMesh")
    log = open(f"{p}/log.checkMesh").read()
    nc = re.search(r"cells:\s+(\d+)", log).group(1)
    check("Mesh OK." in log, f"checkMesh: Mesh OK ({nc} cells)")

    # 4. mass placed
    m = re.search(r"([\d.]+) kg placed \(geometry ([\d.]+) kg\)", r.stdout)
    placed, geo = float(m.group(1)), float(m.group(2))
    check(abs(placed / geo - 1) < 2e-3 and abs(geo - 1503) < 2,
          f"fill places {placed:.1f} kg, geometry {geo:.1f} kg (Couto 1503 kg)")

    # 5. short run, then resume
    subprocess.run(["foamDictionary", "-entry", "endTime", "-set", "3e-4", "system/controlDict"],
                   cwd=p, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["foamDictionary", "-entry", "writeInterval", "-set", "3e-4",
                    "system/controlDict"], cwd=p, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    txt = open(f"{p}/Allrun").read().replace("NP=4", "NP=1")
    open(f"{p}/Allrun", "w").write(txt)
    env = dict(os.environ, PWD=p)
    rc = subprocess.run(["./Allrun"], cwd=p, env=env).returncode
    log = open(f"{p}/log.celFoam").read()
    mass = [float(x) for x in re.findall(r"Total mass\s+=\s+(\S+)", log)]
    steps = len(re.findall(r"^Time = ", log, flags=re.M))
    ok = rc == 0 and "End" in log and steps > 20 and len(mass) >= 2
    check(ok, f"celFoam runs to 0.3 ms ({steps} steps)")
    if mass:
        drift = abs(mass[-1] / mass[0] - 1)
        check(drift < 1e-9, f"mass kept to {drift:.1e}")
    subprocess.run(["foamDictionary", "-entry", "endTime", "-set", "4e-4", "system/controlDict"],
                   cwd=p, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["./Allrun"], cwd=p, env=env)
    log = open(f"{p}/log.celFoam").read()
    check("restart state read from 0.0003" in log, "second Allrun resumes from the last write")
    print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
