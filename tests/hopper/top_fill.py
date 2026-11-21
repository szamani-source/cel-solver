"""Couto silo filled through a narrow inlet at the top, outlet closed, then
left to settle; the at-rest wall pressures are read after settling (Sia,
2026-11-20). MC with the T0 rules, on the wall-aligned mesh with the STL
shell inside it (buffered_case.aligned_buffered_case, WP G.5).

- The silo starts empty. The inlet is the part of the mesh top with
  r < R_IN (a disc narrower than the silo), split off the top patch with
  topoSet/createPatch, R_HEAD of headroom above the silo top so the heap
  never reaches the inlet (an inlet touching the heap would press on it).
- The inlet feeds V_IN m/s at the bulk density until the Couto fill mass
  (the silo full to its top at rho0) has gone in. Gravity is on from the
  start (0.1 s ramp) and there is no damping while the material falls, so
  the stream lands at its free-fall speed.
- Damping (100 1/s) starts T_FALL after the inlet closes, when the last
  material has landed (controlDict dampingStart), and the silo settles to
  T_SETTLE after the inlet closes. Pressures at three late writes show
  whether they are steady.

    CEL_NPROCS=4 python3 top_fill.py WORKDIR

Environment: CEL_RIN (m, 0.25), CEL_VIN (m/s, 3), CEL_HEAD (m, 0.40),
CEL_SETTLE (s, 1.0), CEL_FILL_MASS (kg, the Couto mass): the mass the inlet
feeds. The inlet time is set from the inlet patch's real radius on the mesh
(the faces whose centre lies inside R_IN, so the outer edge of the last
column): the nominal pi R_IN^2 put in 1541 kg for 1503 (decision 0114).
"""
import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hopper_case as hc  # noqa: E402
import hopper_pressure_checks as hp  # noqa: E402
import buffered_case as bc  # noqa: E402

vc = hp.vc
cs = bc.cs
NP = int(os.environ.get("CEL_NPROCS", "1"))
R_IN = float(os.environ.get("CEL_RIN", 0.25))
V_IN = float(os.environ.get("CEL_VIN", 3.0))
R_HEAD = float(os.environ.get("CEL_HEAD", 0.40))
H_TOP = 0.025
T_FALL = 0.3
T_SETTLE = float(os.environ.get("CEL_SETTLE", 1.0))   # s after the inlet closes
CELLS = ["GA0", "GA1", "GA2", "GA3", "GA4", "GA5", "GA6", "GA7"]


def split_inlet(p, r_in):
    """Move the top faces with r < r_in into a new patch `inlet`."""
    xr = r_in * math.cos(math.radians(hc.THETA))
    cs.write(f"{p}/system/topoSetDict", cs.HEADER.format(cls="dictionary", obj="topoSetDict")
             + "actions\n(\n"
             "    { name inletFaces; type faceSet; action new; source patchToFace; patch top; }\n"
             "    { name inletFaces; type faceSet; action subset; source boxToFace;\n"
             f"      box (-1 -100 -1) ({xr!r} 100 1); }}\n);\n")
    cs.write(f"{p}/system/createPatchDict", cs.HEADER.format(cls="dictionary", obj="createPatchDict")
             + "pointSync false;\npatches\n(\n    {\n        name inlet;\n"
             "        patchInfo { type patch; }\n        constructFrom set;\n"
             "        set inletFaces;\n    }\n);\n")


def sh(p, cmd, log, env=None):
    with open(f"{p}/log.{log}", "w") as f:
        r = subprocess.run(cmd, cwd=p, stdout=f, stderr=subprocess.STDOUT,
                           env=dict(env or os.environ, PWD=p))
    if r.returncode:
        raise RuntimeError(f"{cmd[0]} failed:\n" + open(f"{p}/log.{log}").read()[-3000:])


def mesh(p):
    sh(p, ["blockMesh"], "blockMesh")
    sh(p, ["topoSet"], "topoSet")
    sh(p, ["createPatch", "-overwrite"], "createPatch")


def inlet_radius(p):
    """Outer radius of the inlet patch on the mesh (largest r of its points)."""
    import re
    sys.path.insert(0, os.path.join(HERE, "..", "angled"))
    import angled_checks as ac
    m = f"{p}/constant/polyMesh"
    b = open(f"{m}/boundary").read()
    blk = b[b.index("inlet"):]
    n = int(re.search(r"nFaces\s+(\d+)", blk).group(1))
    s0 = int(re.search(r"startFace\s+(\d+)", blk).group(1))
    pts = ac.foam_list(f"{m}/points", "point")
    faces = ac.foam_list(f"{m}/faces", "face")
    ids = sorted({i for f in faces[s0:s0 + n] for i in f})
    return float(pts[ids, 0].max() / math.cos(math.radians(hc.THETA)))


def run(p):
    def sh_(cmd, log, env=None):
        sh(p, cmd, log, env)
    if NP > 1:
        sh_(["decomposePar", "-force"], "decomposePar")
        env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
        sh_(["mpirun", "-np", str(NP), "--oversubscribe", "celVoidCheck", "-parallel"],
            "driver", env)
    else:
        sh_(["celVoidCheck"], "driver")


def main():
    t0 = time.perf_counter()
    work = os.path.abspath(sys.argv[1])
    c = dict(hc.COUTO)
    b, ht = hc.hopper_geometry(c)
    H = ht + c["Hc"]
    vol = (math.pi * c["R"] ** 2 * c["Hc"]
           + math.pi * ht / 3 * (c["R"] ** 2 + c["R"] * c["Rout"] + c["Rout"] ** 2))
    # the spout below the lip (Y0 to 0) is part of the Couto fill too
    vol += math.pi * c["Rout"] ** 2 * (-bc.Y0)
    mass = float(os.environ.get("CEL_FILL_MASS", vol * c["rho"]))
    # mesh once with a dummy time to read the inlet's real radius
    q, _ = bc.aligned_buffered_case(work, "inlet_probe", "mc",
                                    dict(tEnd=1e-6, fill=[-1.0, -1.0, 1.0, bc.Y0 - 1.0]),
                                    dict(c), head_rows=int(round(R_HEAD / H_TOP)))
    split_inlet(q, R_IN)
    mesh(q)
    r_mesh = inlet_radius(q)
    t_fill = mass / (c["rho"] * V_IN * math.pi * r_mesh ** 2)
    t_damp = round(t_fill + T_FALL, 3)
    t_end = round(t_fill + T_SETTLE, 2)
    writes = [round(t_fill + 0.2, 2), round(t_end - 0.2, 2), t_end]
    os.environ.setdefault("CEL_HOPPER_RULES", hp.MU_I + hp.RULES + "fullBlend yes;\n")
    check = dict(tEnd=t_end, maxCo=0.4, writeTimes=writes,
                 fill=[-1.0, -1.0, 1.0, bc.Y0 - 1.0])
    p, c = bc.aligned_buffered_case(work, "couto_mc_topfill", "mc", check, c, nProcs=NP,
                                    head_rows=int(round(R_HEAD / H_TOP)))
    mp = f"{p}/constant/materialProperties"
    txt = open(mp).read()
    txt += (f"inlets\n{{\n    inlet\n    {{\n        speed {V_IN!r};\n        startTime 0;\n"
            f"        endTime {t_fill!r};\n        density {c['rho']!r};\n"
            "        pressure 0;\n        friction 0;\n    }\n}\n")
    open(mp, "w").write(txt)
    dp = f"{p}/system/decomposeParDict"
    txt = open(dp).read().replace(f"n ({NP} 1 1)", f"n (1 {NP} 1)")
    open(dp, "w").write(txt)
    open(f"{p}/system/controlDict", "a").write(
        f"gravityRamp 0.1;\ndamping 100;\ndampingStart {t_damp!r};\n")
    split_inlet(p, R_IN)
    mesh(p)
    print(f"inlet r < {R_IN} m (on the mesh {r_mesh:.4f} m) at {V_IN} m/s, closes at {t_fill:.3f} s "
          f"(fill mass {mass:.0f} kg); damping from {t_damp} s; end {t_end} s", flush=True)
    run(p)
    h = vc.history(p)
    f = math.pi / (math.sin(math.radians(hc.THETA)) * math.cos(math.radians(hc.THETA)))
    print(f"run {time.perf_counter() - t0:.0f} s; mass in the silo {h[-1, vc.MASS] * f:.0f} kg "
          f"(target {mass:.0f}); largest speed at the end {h[-1, vc.UMAX]:.3f} m/s")
    hp.CELL = 0.012
    rows = hp.table()
    for w in writes:
        r = hp.pressures(p, c, f"cells_{w:g}")
        print(f"  t {w:g}: " + " ".join(f"{n} {r[n]:.2f}" for n in CELLS), flush=True)
        a = hp.read_cells(p, f"cells_{w:g}")
        out, rr, yy = bc.leak_report(c, a)
        print(f"    leak: {out.sum()} cells, {hp.leak_volume(p, a, out) * c['rho']:.4f} kg")
    print("  experiment at rest: " + " ".join(
        "{} {:.2f}-{:.2f}".format(n, *hp.exp_range(rows, "static", n)) for n in CELLS))
    a = hp.read_cells(p, f"cells_{t_end:g}")
    top = a[(a[:, 2] > 0.5) & (a[:, 10] > 0)]
    print(f"  heap: highest half-full cell {top[:, 1].max():.3f} m (silo top {H:.3f} m)")


if __name__ == "__main__":
    main()
