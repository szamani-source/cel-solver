"""WP G.3 done-when checks: a user-provided mesh, its patches and their
boundary conditions.

1. A Gmsh .msh file (tests/celsetup/data/tank2d.msh, a 2D half tank) is
   converted with gmshToFoam; the patches get their conditions (wall,
   outlet, symmetry, free, empty), the fill box places 0.4 x 0.5 x 0.01 m of
   material (3.0 kg at 1500 kg/m^3) and celFoam runs 50 steps keeping the
   mass.
2. The patch suggestions: frontAndBack empty, axis symmetry.
3. A polyMesh folder imports, with scale.
4. The Couto wedge (written by celSetup with blockMesh) imports as a user
   mesh: wedge angle 2.5 deg found, the same time step from celCheck.
5. Missing and unknown patches give readable errors; the YAML file
   written with the boundaries reads back to the same case.
6. The GUI server loads the mesh and lists the patches.

    python3 tests/celsetup/usermesh_checks.py WORKDIR
"""
import json
import math
import os
import re
import socket
import subprocess
import sys
import time
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from celsetup import case, foam, miniyaml, schema, usermesh  # noqa: E402

fails = []
MSH = os.path.join(HERE, "data", "tank2d.msh")
EX = os.path.join(ROOT, "tools", "celsetup", "examples")


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg, flush=True)
    if not ok:
        fails.append(msg)


def tank_cfg(file=MSH, **extra):
    raw = miniyaml.load(os.path.join(EX, "userTank.yaml"))
    raw["mesh"]["file"] = file
    raw["mesh"].update(extra)
    return raw


def celcheck(c, work, name):
    out = os.path.join(work, name)
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "celCheck"), c, "--json", out],
                   capture_output=True, text=True)
    return json.load(open(out))


def main(work):
    os.makedirs(work, exist_ok=True)
    # 1. Gmsh .msh, 2D
    t = os.path.join(work, "tank")
    case.write_case(case.resolve(tank_cfg()), t)
    ok, msg = case.run_mesh_and_fill(t)
    check(ok, f".msh imported and filled: {msg}")
    if not ok:
        print("FAILED: " + "; ".join(fails))
        return 1
    m = foam.PolyMesh(t)
    types = {n: ty for n, ty, _, _ in usermesh.read_boundary(t)}
    check(types == {"outlet": "patch", "bottom": "patch", "wall": "patch", "top": "patch",
                    "axis": "symmetryPlane", "frontAndBack": "empty"},
          f"patch types set: {types}")
    mp = open(os.path.join(t, "constant", "materialProperties")).read()
    check(re.search(r"walls\s*\{\s*bottom\s*\{ friction 20.0; \}\s*wall\s*\{ friction 20.0; \}", mp)
          and "outlet { friction 20.0; openTime 0.5; }" in mp,
          "materialProperties: walls bottom and wall, outlet opens at run openTime")
    res = case.fill(t)
    check(abs(res["mass"] / 3.0 - 1) < 0.005 and m.nCells == 1080,
          f"fill: {res['mass']:.4f} kg placed (0.4 x 0.5 x 0.01 m x 1500 = 3.0; 8 x 8 samples per cell, within 0.5 %), {m.nCells} cells")
    for k, v in (("endTime", "2e-3"), ("writeInterval", "2e-3")):
        subprocess.run(["foamDictionary", "-entry", k, "-set", v, "system/controlDict"], cwd=t,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rc = subprocess.run(["./Allrun"], cwd=t, env=dict(os.environ, PWD=t)).returncode
    log = open(os.path.join(t, "log.celFoam")).read()
    mass = [float(x) for x in re.findall(r"Total mass\s+=\s+(\S+)", log)]
    steps = len(re.findall(r"^Time = ", log, flags=re.M))
    check(rc == 0 and "End" in log and steps >= 30 and len(mass) >= 2
          and abs(mass[-1] / mass[0] - 1) < 1e-9,
          f"celFoam runs {steps} steps on the imported mesh, mass kept")

    # 2. suggestions
    d = usermesh.describe_mesh(os.path.join(work, "tank"))
    sug = {p["name"]: p["suggest"] for p in d["patches"]}
    check(sug.get("frontAndBack") == "empty" and sug.get("axis") == "symmetry" and d["kind"] == "2D",
          f"suggested conditions: {sug}")

    # 3. polyMesh folder, scaled
    p = os.path.join(work, "poly")
    case.write_case(case.resolve(tank_cfg(os.path.join(t, "constant", "polyMesh"), scale=0.5)), p)
    ok, msg = case.mesh(p)
    hi = foam.PolyMesh(p).points.max(0)
    check(ok and abs(hi[0] - 0.2) < 1e-9 and abs(hi[1] - 0.3) < 1e-9,
          f"polyMesh folder imported at scale 0.5: top corner ({hi[0]:g}, {hi[1]:g}) m")

    # 4. the Couto wedge as a user mesh
    c = os.path.join(work, "couto")
    case.write_case(case.resolve(miniyaml.load(os.path.join(EX, "coutoHopper.yaml"))), c)
    case.mesh(c)
    w = os.path.join(work, "wedge")
    raw = dict(template="userMesh", case="w", mesh=dict(file=c),
               geometry=dict(fillYMax=2.0),
               boundaries=dict(outlet=dict(type="outlet", friction=11.3),
                               outer=dict(type="wall", friction=11.3),
                               bottom=dict(type="wall", friction=11.3), top="free",
                               front="wedge", back="wedge", axis="empty"),
               material=dict(rho=856, E=10.674e6, phi=30, d=0.004))
    case.write_case(case.resolve(raw), w)
    ok, msg = case.run_mesh_and_fill(w)
    th = json.load(open(os.path.join(w, "system", "celSetup.json"))).get("theta") or 0
    rw, rc_ = celcheck(w, work, "w.json"), celcheck(c, work, "c.json")
    check(ok and abs(math.degrees(th) - 2.5) < 1e-6 and abs(rw["dt"] / rc_["dt"] - 1) < 1e-9,
          f"Couto wedge as a user mesh: {msg}; wedge half angle {math.degrees(th):.4f} deg, "
          f"time step {rw['dt']:.6g} s as celSetup's own case")

    # 5. errors and the YAML round trip
    bad = tank_cfg()
    del bad["boundaries"]["top"]
    b = os.path.join(work, "bad")
    case.write_case(case.resolve(bad), b)
    ok, msg = case.mesh(b)
    check(not ok and "give a condition for top" in msg, f"missing patch: {msg}")
    bad["boundaries"]["roof"] = {"type": "free"}
    case.write_case(case.resolve(bad), b)
    ok, msg = case.mesh(b)
    check(not ok and "no patch roof" in msg, f"unknown patch: {msg}")
    try:
        case.resolve(dict(tank_cfg(), boundaries={"wall": {"type": "wall"}}))
        case.write_case(case.resolve(dict(tank_cfg(), boundaries={"wall": {"type": "wall"}})), b)
        err = ""
    except ValueError as e:
        err = str(e)
    check("'friction' is required for a wall" in err, f"wall without friction: {err}")
    cfg = tank_cfg()
    y = schema.template_yaml("userMesh", usermesh.SCHEMA, cfg, boundary_types=usermesh.BTYPES)
    back = miniyaml.loads(y)
    check(case.resolve(back)["boundaries"] == case.resolve(cfg)["boundaries"],
          "YAML written with the boundaries reads back the same")

    # 6. GUI: load the mesh, list patches
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    srv = subprocess.Popen([sys.executable, os.path.join(ROOT, "tools", "celSetup"), "--gui",
                            "--no-browser", "--port", str(port)], cwd=work)
    base = f"http://127.0.0.1:{port}"

    def post(path, body):
        req = urllib.request.Request(base + path, json.dumps(body).encode(),
                                     {"Content-Type": "application/json"})
        return json.loads(urllib.request.urlopen(req, timeout=120).read())
    try:
        for _ in range(50):
            try:
                init = json.loads(urllib.request.urlopen(base + "/api/init", timeout=5).read())
                break
            except OSError:
                time.sleep(0.2)
        check("userMesh" in init["schema"] and "wall" in init["btypes"],
              "GUI offers the userMesh template and the boundary conditions")
        form = tank_cfg()
        form["boundaries"] = {}
        r = post("/api/loadmesh", {"config": form})
        names = sorted(p["name"] for p in r.get("patches", []))
        check(r.get("nCells") == 1080 and len(r.get("polys", [])) == 1080
              and names == ["axis", "bottom", "frontAndBack", "outlet", "top", "wall"],
              f"GUI loads the mesh: {r.get('nCells')} cells, patches {names}")
        m = post("/api/mesh", {"config": tank_cfg()})
        check(m.get("nCells") == 1080 and set(m.get("patches") or {}) ==
              {"axis", "bottom", "outlet", "top", "wall"},
              f"GUI preview with every patch given: verdict {m.get('verdict')}, "
              f"patches drawn {sorted(m.get('patches') or {})}")
    finally:
        srv.terminate()
        srv.wait()
    print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
