"""Write a complete celFoam case from a celSetup configuration, and fill
the material (0/alpha, 0/rho) once the mesh exists."""

import json
import math
import os
import shutil
import subprocess

import numpy as np

from . import conical, foam, schema

TEMPLATES = {"conicalHopper": conical}
TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve(cfg):
    """Defaults, presets and derived values filled in; raises ValueError
    with a readable message on a missing or wrong entry."""
    t = cfg.get("template")
    if t not in TEMPLATES:
        raise ValueError(f"template must be one of {', '.join(TEMPLATES)} (got {t!r})")
    mod = TEMPLATES[t]
    given = {k: dict(v) if isinstance(v, dict) else v for k, v in cfg.items()}
    mesh = dict(given.get("mesh") or {})
    preset = mesh.get("preset", "standard")
    if preset not in mod.PRESETS:
        raise ValueError(f"mesh: preset must be one of {', '.join(mod.PRESETS)}")
    given["mesh"] = dict(mod.PRESETS[preset],
                         **{k: v for k, v in mesh.items() if v not in ("preset", None)})
    out = {"template": t, "case": cfg.get("case", "celCase")}
    for sec, entries in schema.sections(mod.SCHEMA).items():
        known = {e.key for e in entries}
        extra = set(given.get(sec) or {}) - known
        if extra and sec != "material":       # material passes others verbatim
            raise ValueError(f"{sec}: unknown entr{'y' if len(extra) == 1 else 'ies'} "
                             f"{', '.join(sorted(extra))} (celSetup --describe {t} lists them)")
        out[sec] = schema.apply_defaults(given, sec, entries)
    w = out["walls"]
    if w.get("friction") is None:
        if w.get("frictionCoefficient") is None:
            raise ValueError("walls: give friction (deg) or frictionCoefficient")
        w["friction"] = math.degrees(math.atan(float(w["frictionCoefficient"])))
    return out


def material_text(mat):
    phi = float(mat["phi"])
    lines = [("type", mat["model"]), ("E", float(mat["E"])), ("nu", float(mat["nu"])),
             ("phi", phi), ("c", float(mat["c"])), ("psi", float(mat["psi"]))]
    if mat["model"] == "mohrCoulomb":
        st = mat["sigmaT"]
        if st == "auto":
            c = float(mat["c"])
            st = 0.0 if c == 0 else c / math.tan(math.radians(phi))
        lines.append(("sigmaT", float(st)))
    else:
        lines.append(("fit", mat["fit"]))
    mu2 = mat["mu2"]
    mu2 = math.tan(math.radians(phi)) + 0.26 if mu2 == "auto" else float(mu2)
    rhoS = mat["rhoS"]
    rhoS = float(mat["rho"]) / 0.62 if rhoS == "auto" else float(rhoS)
    lines += [("mu2", round(mu2, 6)), ("I0", float(mat["I0"])), ("d", float(mat["d"])),
              ("rhoS", round(rhoS, 3))]
    for k in ("muIAlpha", "nonlocalA", "capLambda"):
        if mat.get(k) is not None:
            lines.append((k, float(mat[k])))
    if mat.get("capLambda") is not None:
        lines.append(("capP0", float(mat["capP0"])))
    known = {e.key for e in schema.MATERIAL}
    for k, v in mat.items():                       # anything else, verbatim
        if k not in known:
            lines.append((k, v))
    on = mat["rules"] == "recommended"
    for k in schema.SWITCHES:
        v = mat[k]
        lines.append((k, on if v == "auto" else v in (True, "yes")))
    if mat["voidPlasticFlow"] in (True, "yes"):
        lines.append(("voidPlasticFlow", True))
    return "".join(f"{k:<16}{foam.fmt(v)};\n" for k, v in lines)


def write_case(cfg, out_dir, src="celSetup"):
    """Write every file of the case. Returns the template's build result."""
    t = TEMPLATES[cfg["template"]]
    b = t.build(cfg, src)
    run, mat, w = cfg["run"], cfg["material"], cfg["walls"]
    os.makedirs(out_dir, exist_ok=True)
    for rel, txt in b["files"].items():
        foam.write(os.path.join(out_dir, rel), txt)
    H = lambda cls, obj: foam.header(cls, obj, src)  # noqa: E731
    fr = float(w["friction"])
    walls = "".join(f"    {p:<10}{{ friction 0; }}\n" for p in b["boundaryWalls"])
    outlets = "".join(f"    {p} {{ friction {fr!r}; openTime {float(run['openTime'])!r}; }}\n"
                      for p in b["outlets"])
    foam.write(f"{out_dir}/constant/materialProperties",
               H("dictionary", "materialProperties")
               + f"rho             {float(mat['rho'])!r};\nfreeSurface     yes;\n\n"
               + material_text(mat)
               + "\n// domain sides outside the STL wall: frictionless\nwalls\n{\n" + walls + "}\n"
               + "outlets\n{\n" + outlets + "}\n")
    bodies = "".join(f"    {n}\n    {{\n        surface {s};\n        friction {fr!r};\n    }}\n"
                     for n, s in b["bodies"].items())
    foam.write(f"{out_dir}/constant/rigidBodies",
               H("dictionary", "rigidBodies") + "bodies\n{\n" + bodies + "}\n")
    foam.write(f"{out_dir}/constant/g", H("uniformDimensionedVectorField", "g")
               + f"dimensions      [0 1 -2 0 0 0 0];\nvalue           (0 {-float(run['gravity'])!r} 0);\n")
    foam.write(f"{out_dir}/system/controlDict", H("dictionary", "controlDict") + f"""
application     celFoam;
startFrom       latestTime;     // a second Allrun resumes from the last write
startTime       0;
stopAt          endTime;
endTime         {float(run['endTime'])!r};
deltaT          1e-6;           // celFoam sets the stable step itself
writeControl    adjustableRunTime;
writeInterval   {float(run['writeInterval'])!r};
purgeWrite      0;
writeFormat     ascii;
writePrecision  10;
timeFormat      general;
timePrecision   8;
runTimeModifiable false;
maxCo           {float(run['maxCo'])!r};
gravityRamp     {float(run['gravityRamp'])!r};
damping         {float(run['damping'])!r};
dampingStart    {float(run['dampingStart'])!r};
""")
    foam.write(f"{out_dir}/system/fvSchemes", H("dictionary", "fvSchemes")
               + "ddtSchemes { default Euler; }\ngradSchemes { default Gauss linear; }\n"
               "divSchemes { default none; }\nlaplacianSchemes { default none; }\n"
               "interpolationSchemes { default linear; }\nsnGradSchemes { default corrected; }\n")
    foam.write(f"{out_dir}/system/fvSolution", H("dictionary", "fvSolution"))
    n = int(run["nProcs"])
    foam.write(f"{out_dir}/system/decomposeParDict", H("dictionary", "decomposeParDict")
               + f"numberOfSubdomains {n};\nmethod simple;\ncoeffs {{ n {b['decompose'].format(n=n)}; }}\n")
    for f in ("alpha", "rho"):
        foam.write(f"{out_dir}/0.orig/{f}", H("volScalarField", f)
                   + "dimensions      " + ("[0 0 0 0 0 0 0]" if f == "alpha" else "[1 -3 0 0 0 0 0]")
                   + ";\ninternalField   uniform 0;\nboundaryField\n{\n"
                   "    #includeEtc \"caseDicts/setConstraintTypes\"\n"
                   "    \".*\" { type zeroGradient; }\n}\n")
    meta = dict(b["notes"], wedge=b["wedge"], theta=b["theta"], fill=b["fill"],
                solids=b["solids"], volume=b["volume"], rho=float(mat["rho"]),
                mass=b["volume"] * float(mat["rho"]), d=float(mat["d"]),
                nonlocalA=mat.get("nonlocalA"), E=float(mat["E"]), nu=float(mat["nu"]),
                maxCo=float(run["maxCo"]), endTime=float(run["endTime"]), nProcs=n,
                summary=b["summary"])
    foam.write(f"{out_dir}/system/celSetup.json", json.dumps(meta, indent=1))
    setup = os.path.join(TOOLS, "celSetup")
    foam.write(f"{out_dir}/Allrun", f"""#!/bin/sh
# Written by celSetup. Runs the case, and resumes from the last write when
# run again (checkpoints every writeInterval; startFrom latestTime).
cd "${{0%/*}}" || exit
CELSETUP="${{CELSETUP:-{setup}}}"
command -v celSetup >/dev/null 2>&1 && [ ! -x "$CELSETUP" ] && CELSETUP=celSetup
NP={n}
if [ ! -f constant/polyMesh/faces ]; then
    blockMesh > log.blockMesh 2>&1 || {{ echo "blockMesh failed, see log.blockMesh"; exit 1; }}
fi
if [ ! -f 0/alpha ]; then
    python3 "$CELSETUP" --fill . || exit 1
fi
if [ "$NP" -gt 1 ]; then
    [ -d processor0 ] || decomposePar -force > log.decomposePar 2>&1 || exit 1
    mpirun -np "$NP" celFoam -parallel >> log.celFoam 2>&1
else
    celFoam >> log.celFoam 2>&1
fi
""", 0o755)
    foam.write(f"{out_dir}/Allclean", """#!/bin/sh
# Removes the mesh, results and logs; keeps the setup.
cd "${0%/*}" || exit
. ${WM_PROJECT_DIR:?}/bin/tools/CleanFunctions
cleanCase
rm -rf 0
""", 0o755)
    return b


# --- material fill ----------------------------------------------------------

def _inside(pts, poly):
    """Even-odd point in polygon, vectorised: pts (n, 2), poly (m, 2)."""
    x, y = pts[:, 0][:, None], pts[:, 1][:, None]
    a = poly
    b = np.roll(poly, -1, axis=0)
    ay, by = a[:, 1][None, :], b[:, 1][None, :]
    ax, bx = a[:, 0][None, :], b[:, 0][None, :]
    cond = (ay > y) != (by > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = ax + (y - ay) * (bx - ax) / (by - ay)
    return (cond & (x < xi)).sum(1) % 2 == 1


def fill(case, nsub=8):
    """Write 0/alpha and 0/rho: alpha is the material fraction of the open
    part of each cell (the part outside the STL bodies), from nsub x nsub
    sample points per cell (weighted by radius on a wedge)."""
    meta = json.load(open(os.path.join(case, "system", "celSetup.json")))
    mesh = foam.PolyMesh(case)
    cells, polys, wedge = mesh.front_polygons()
    fillp = np.array(meta["fill"])
    solids = [np.array(s) for s in meta["solids"]]
    u = (np.arange(nsub) + 0.5) / nsub
    U, Vv = np.meshgrid(u, u)
    U, Vv = U.ravel(), Vv.ravel()
    alpha = np.zeros(mesh.nCells)
    openf = np.ones(mesh.nCells)
    fb = (fillp.min(0), fillp.max(0))
    for c, q in zip(cells, polys):
        if len(q) == 3:                     # wedge axis prism seen as a triangle
            q = np.vstack([q, q[-1]])
        lo, hi = q.min(0), q.max(0)
        touches_fill = not (hi[0] < fb[0][0] or lo[0] > fb[1][0] or hi[1] < fb[0][1] or lo[1] > fb[1][1])
        sb = [s for s in solids if not (hi[0] < s[:, 0].min() or lo[0] > s[:, 0].max()
                                        or hi[1] < s[:, 1].min() or lo[1] > s[:, 1].max())]
        if not touches_fill and not sb:
            continue
        # bilinear map of the quad (corner order as in the face)
        p = ((1 - U) * (1 - Vv))[:, None] * q[0] + (U * (1 - Vv))[:, None] * q[1] \
            + (U * Vv)[:, None] * q[2] + ((1 - U) * Vv)[:, None] * q[3]
        w = np.maximum(p[:, 0], 0) if wedge else np.ones(len(p))
        if w.sum() <= 0:
            w = np.ones(len(p))
        solid = np.zeros(len(p), bool)
        for s in sb:
            solid |= _inside(p, s)
        mat = _inside(p, fillp) & ~solid if touches_fill else np.zeros(len(p), bool)
        of = (w * ~solid).sum() / w.sum()
        openf[c] = of
        alpha[c] = (w * mat).sum() / (w * ~solid).sum() if of > 0 else 0.0
    alpha = np.where(alpha > 1 - 1e-9, 1.0, np.where(alpha < 1e-9, 0.0, alpha))
    rho = meta["rho"] * alpha
    os.makedirs(os.path.join(case, "0"), exist_ok=True)
    for name, val in (("alpha", alpha), ("rho", rho)):
        src = open(os.path.join(case, "0.orig", name)).read()
        body = "nonuniform List<scalar>\n%d\n(\n%s\n)" % (len(val), "\n".join(f"{v:.10g}" for v in val))
        open(os.path.join(case, "0", name), "w").write(
            src.replace("internalField   uniform 0", "internalField   " + body))
    # full-geometry mass placed (wedge: scale by 2 pi / wedge angle)
    vol = (alpha * openf * np.abs(mesh.V)).sum()
    if wedge:
        # flat-faced wedge cell: r dr dy 2 sin(t) cos(t)
        vol *= math.pi / (math.sin(meta["theta"]) * math.cos(meta["theta"]))
    return dict(mass=vol * meta["rho"], target=meta["mass"], cells=int((alpha > 0).sum()))


def run_mesh_and_fill(case):
    """blockMesh, then the fill. Returns (ok, message)."""
    if shutil.which("blockMesh") is None:
        return False, "blockMesh not found: source OpenFOAM's bashrc, then run ./Allrun"
    with open(os.path.join(case, "log.blockMesh"), "w") as f:
        r = subprocess.run(["blockMesh"], cwd=case, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        return False, "blockMesh failed, see " + os.path.join(case, "log.blockMesh")
    res = fill(case)
    return True, (f"filled {res['cells']} cells: {res['mass']:.1f} kg placed "
                  f"(geometry {res['target']:.1f} kg)")
