"""WP G.1 done-when checks: celSetup meshes the Couto hopper with Gmsh.

1. With `mesh: mesher: gmsh` celSetup writes system/mesh.geo (no
   blockMeshDict), and Allrun meshes through `celSetup --mesh`.
2. (Gmsh installed) Gmsh plus gmshToFoam give a mesh with the same cell
   count and the same patches and patch types as blockMesh; checkMesh
   reports Mesh OK.
3. The fill places the Couto mass within 0.2 %; celCheck's time step is
   within 1 % of the blockMesh mesh's.
4. celFoam runs 0.3 ms on the Gmsh mesh and keeps the mass to 1e-9.

Without Gmsh (the CI image) steps 2 to 4 are skipped and said so.

    python3 tests/celsetup/gmsh_checks.py WORKDIR
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from celsetup import case, gmsh_mesh, miniyaml  # noqa: E402

fails = []


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg, flush=True)
    if not ok:
        fails.append(msg)


def patches(c):
    txt = open(os.path.join(c, "constant", "polyMesh", "boundary")).read()
    return {n: (t, int(f)) for n, t, f in re.findall(
        r"\n\s*(\w+)\s*\{[^}]*?type\s+(\w+);[^}]*?nFaces\s+(\d+);", txt)}


def check_json(c, work, name):
    out = os.path.join(work, name)
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "celCheck"), c, "--json", out],
                   capture_output=True, text=True)
    return json.load(open(out))


def main(work):
    os.makedirs(work, exist_ok=True)
    ex = os.path.join(ROOT, "tools", "celsetup", "examples", "coutoHopper.yaml")
    raw = miniyaml.load(ex)
    raw["mesh"]["mesher"] = "gmsh"
    gm = os.path.join(work, "gmsh")
    case.write_case(case.resolve(raw), gm)
    geo = open(os.path.join(gm, "system", "mesh.geo")).read()
    names = set(re.findall(r'Physical \w+\("(\w+)"\)', geo))
    check(not os.path.exists(os.path.join(gm, "system", "blockMeshDict"))
          and names == {"outlet", "bottom", "outer", "top", "front", "back", "internal"}
          and geo.count("Transfinite Surface") == 6,
          f"mesher gmsh writes system/mesh.geo: 6 transfinite blocks, groups {sorted(names)}")
    check('"$CELSETUP" --mesh .' in open(os.path.join(gm, "Allrun")).read(),
          "Allrun meshes through celSetup --mesh")
    if gmsh_mesh.available() is None or os.environ.get("WM_PROJECT_DIR") is None:
        print("SKIPPED: Gmsh or OpenFOAM not available, meshing not checked here")
        print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
        return 1 if fails else 0

    bm = os.path.join(work, "blockMesh")
    case.write_case(case.resolve(miniyaml.load(ex)), bm)
    for c in (bm, gm):
        ok, msg = case.run_mesh_and_fill(c)
        check(ok, f"{os.path.basename(c)}: {msg}")
    pb, pg = patches(bm), patches(gm)
    pb.pop("axis", None)                          # blockMesh's collapsed axis, no faces
    check(pb == pg, f"same patches, types and face counts as blockMesh: {pg}")
    with open(os.path.join(gm, "log.checkMesh"), "w") as f:
        subprocess.run(["checkMesh"], cwd=gm, stdout=f, stderr=subprocess.STDOUT,
                       env=dict(os.environ, PWD=gm))
    log = open(os.path.join(gm, "log.checkMesh")).read()
    nc = re.search(r"cells:\s+(\d+)", log).group(1)
    check("Mesh OK." in log and nc == "7584", f"checkMesh: Mesh OK, {nc} cells (blockMesh 7584)")
    rb, rg = check_json(bm, work, "bm.json"), check_json(gm, work, "gm.json")
    check(abs(rg["dt"] / rb["dt"] - 1) < 0.01,
          f"celCheck time step {rg['dt']:.4g} s (blockMesh {rb['dt']:.4g} s)")
    check(rg["nonOrth"]["max"] < 45, f"non-orthogonality max {rg['nonOrth']['max']:.1f} deg "
          f"(blockMesh {rb['nonOrth']['max']:.1f})")

    for k, v in (("endTime", "3e-4"), ("writeInterval", "3e-4")):
        subprocess.run(["foamDictionary", "-entry", k, "-set", v, "system/controlDict"],
                       cwd=gm, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    txt = open(f"{gm}/Allrun").read().replace("NP=4", "NP=1")
    open(f"{gm}/Allrun", "w").write(txt)
    rc = subprocess.run(["./Allrun"], cwd=gm, env=dict(os.environ, PWD=gm)).returncode
    log = open(f"{gm}/log.celFoam").read()
    mass = [float(x) for x in re.findall(r"Total mass\s+=\s+(\S+)", log)]
    steps = len(re.findall(r"^Time = ", log, flags=re.M))
    check(rc == 0 and "End" in log and steps > 20 and len(mass) >= 2,
          f"celFoam runs to 0.3 ms on the Gmsh mesh ({steps} steps)")
    if len(mass) >= 2:
        drift = abs(mass[-1] / mass[0] - 1)
        check(drift < 1e-9, f"mass kept to {drift:.1e}")
    print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
