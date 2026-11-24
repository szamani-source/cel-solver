"""User-provided meshes (WP G.3): template `userMesh`.

The mesh comes from the user: an OpenFOAM case or polyMesh folder, or a
file in a format OpenFOAM converts (.msh from Gmsh or Fluent, .geo meshed
with Gmsh, .unv, .neu, .ans). celSetup imports it into the case, lists its
patches, and the user gives each patch a boundary condition:

    boundaries:
      hopperWall: {type: wall, friction: 20}      # or as an indented block
      outlet:     {type: outlet, friction: 20}
      axis:       {type: symmetry}
      front:      {type: wedge}

The material fills the cells inside a box (fillXMin ... fillYMax; on a
wedge x is the distance from the axis). The mesh must be 2D (one cell thick,
front and back `empty`) or an axisymmetric wedge about the y axis: the fill
and the checker work on its (x, y) cross-section.
"""

import json
import math
import os
import re
import shutil
import subprocess

import numpy as np

from . import foam
from .schema import AUTO, REQUIRED, E

SCHEMA = {
    "mesh": [
        E("file", REQUIRED, "", "the mesh: an OpenFOAM case or polyMesh folder, or a file: "
          ".msh (Gmsh or Fluent), .geo (meshed with Gmsh), .unv (I-DEAS), .neu (Gambit), "
          ".ans (ANSYS)"),
        E("format", "auto", "", "mesh format; auto: from the file name and contents",
          ["auto", "polyMesh", "gmsh", "geo", "fluent", "unv", "gambit", "ansys"]),
        E("scale", 1.0, "", "multiply the coordinates by this (0.001 for a mesh drawn in mm)"),
    ],
    "geometry": [
        E("fillXMin", None, "m", "material fills the cells inside this box; empty: no limit "
          "(x is the distance from the axis on a wedge)"),
        E("fillXMax", None, "m", "right side of the fill box; empty: no limit"),
        E("fillYMin", None, "m", "bottom of the fill box; empty: no limit"),
        E("fillYMax", None, "m", "top of the fill box (the fill level); empty: no limit"),
    ],
}
PRESETS = {"standard": {}}
BOUNDARIES = True

# Boundary conditions a patch can take, with their entries
BTYPES = {
    "wall": ("a wall: the material slides along it with Coulomb friction",
             [E("friction", REQUIRED, "deg", "wall friction angle (0: frictionless)")]),
    "outlet": ("closed (a wall) until openTime, then open: material leaves and is counted",
               [E("friction", REQUIRED, "deg", "wall friction angle while closed"),
                E("openTime", AUTO, "s", "opens at this time; auto: run openTime")]),
    "inlet": ("feeds fresh material into the domain at a set speed",
              [E("speed", REQUIRED, "m/s", "speed of the material into the domain"),
               E("startTime", 0.0, "s", "feeding starts"),
               E("endTime", None, "s", "feeding stops; empty: never"),
               E("density", AUTO, "kg/m^3", "density of the fed material; auto: material rho"),
               E("pressure", 0.0, "Pa", "pressure of the fed material"),
               E("friction", 0.0, "deg", "wall friction angle while closed")]),
    "symmetry": ("a symmetry plane: no normal velocity, no shear", []),
    "wedge": ("a side of an axisymmetric wedge (front or back)", []),
    "empty": ("a side of a 2D mesh (front or back)", []),
    "free": ("open and traction-free: nothing holds the material there", []),
}
FOAM_TYPE = {"wall": "patch", "outlet": "patch", "inlet": "patch", "symmetry": "symmetryPlane",
             "wedge": "wedge", "empty": "empty", "free": "patch"}
CONVERTERS = {"gmsh": "gmshToFoam", "fluent": "fluentMeshToFoam", "unv": "ideasUnvToFoam",
              "gambit": "gambitToFoam", "ansys": "ansysToFoam"}


def detect_format(path):
    if os.path.isdir(path):
        return "polyMesh"
    ext = os.path.splitext(path)[1].lower()
    if ext == ".msh":
        with open(path, "rb") as f:
            head = f.read(200).decode(errors="replace")
        return "gmsh" if "$MeshFormat" in head else "fluent"
    fmt = {".geo": "geo", ".unv": "unv", ".neu": "gambit", ".ans": "ansys"}.get(ext)
    if fmt is None:
        raise ValueError(f"mesh: cannot tell the format of {path}; give mesh: format")
    return fmt


def polymesh_dir(path):
    for p in (path, os.path.join(path, "polyMesh"), os.path.join(path, "constant", "polyMesh")):
        if os.path.exists(os.path.join(p, "faces")) and os.path.exists(os.path.join(p, "owner")):
            return p
    raise ValueError(f"mesh: no polyMesh (points, faces, owner, ...) in {path}")


def _run(cmd, case, log):
    with open(os.path.join(case, log), "w") as f:
        r = subprocess.run(cmd, cwd=case, stdout=f, stderr=subprocess.STDOUT,
                           env=dict(os.environ, PWD=os.path.abspath(case)))
    if r.returncode != 0:
        raise ValueError(f"{cmd[0]} failed, see {os.path.join(case, log)}")


def import_mesh(case, spec):
    """Copy or convert the user's mesh into case/constant/polyMesh and scale
    it. The case needs system/controlDict (write_case has written it)."""
    src = os.path.abspath(os.path.expanduser(spec["file"]))
    if not os.path.exists(src):
        raise ValueError(f"mesh: {src} not found")
    fmt = spec.get("format") or "auto"
    fmt = detect_format(src) if fmt == "auto" else fmt
    dst = os.path.join(case, "constant", "polyMesh")
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    if fmt == "polyMesh":
        shutil.copytree(polymesh_dir(src), dst, ignore=shutil.ignore_patterns("sets", "*.gz"))
    else:
        if fmt == "geo":
            from . import gmsh_mesh
            if gmsh_mesh.available() is None:
                raise ValueError("mesh: a .geo file needs Gmsh: pip install gmsh")
            shutil.copy(src, os.path.join(case, "system", "mesh.geo"))
            msh = os.path.join(case, "mesh.msh")
            r = subprocess.run(["gmsh", "-3", "system/mesh.geo", "-o", "mesh.msh", "-format",
                                "msh22"], cwd=case, capture_output=True) \
                if gmsh_mesh.available() == "command" else None
            if r is None:
                import gmsh
                gmsh.initialize(["gmsh"], interruptible=False)
                try:
                    gmsh.option.setNumber("General.Terminal", 0)
                    gmsh.open(os.path.join(case, "system", "mesh.geo"))
                    gmsh.model.mesh.generate(3)
                    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
                    gmsh.write(msh)
                except Exception as e:  # gmsh raises plain Exceptions
                    raise ValueError(f"Gmsh failed on {src}: {e}")
                finally:
                    gmsh.finalize()
            elif r.returncode != 0:
                raise ValueError(f"Gmsh failed on {src}")
            src, fmt = msh, "gmsh"
        tool = CONVERTERS[fmt]
        if shutil.which(tool) is None:
            raise ValueError(f"mesh: {tool} not found: source OpenFOAM's bashrc")
        _run([tool, src], case, f"log.{tool}")
        if os.path.exists(os.path.join(case, "mesh.msh")):
            os.remove(os.path.join(case, "mesh.msh"))
    for d in ("sets",):
        shutil.rmtree(os.path.join(dst, d), ignore_errors=True)
    scale = float(spec.get("scale") or 1.0)
    if scale != 1.0:
        _run(["transformPoints", "-scale", repr(scale)], case, "log.transformPoints")
    return fmt


def read_boundary(case):
    """[(name, type, nFaces, startFace)] from constant/polyMesh/boundary."""
    txt = open(os.path.join(case, "constant", "polyMesh", "boundary")).read()
    out = []
    for name, body in re.findall(r"\n\s*([\w.:-]+)\s*\{([^}]*)\}", txt):
        t = re.search(r"\btype\s+(\w+);", body)
        n = re.search(r"\bnFaces\s+(\d+);", body)
        s = re.search(r"\bstartFace\s+(\d+);", body)
        if t and n and s:
            out.append((name, t.group(1), int(n.group(1)), int(s.group(1))))
    return out


def suggest(m, name, ftype):
    """A likely boundary condition for a patch from its type and normals:
    wedge, empty or symmetry when it says so or the faces show it."""
    if ftype in ("wedge", "empty"):
        return ftype
    if ftype in ("symmetryPlane", "symmetry"):
        return "symmetry"
    f = m.patch_faces(name)
    if not len(f):
        return "free"
    n = m.Sf[f] / np.maximum(m.magSf[f], 1e-300)[:, None]
    if np.all(np.abs(n[:, 2]) > 0.999):
        return "empty"
    if np.all(np.abs(n[:, 1]) < 1e-6) and np.all(np.abs(n[:, 2]) > 0.5):
        th = np.abs(np.arctan2(n[:, 0], np.abs(n[:, 2])))
        if th.max() - th.min() < 1e-6 and th.max() > 1e-6:
            return "wedge"
    if ftype == "wall":
        return "wall"
    if np.all(np.abs(n[:, 0]) > 0.999) and np.allclose(m.Cf[f, 0], 0, atol=1e-9):
        return "symmetry"
    return None


def describe_mesh(case):
    """What the GUI shows after a load: patches with a suggested condition,
    their edges in the (x, y) section, cell polygons, bounds and kind."""
    m = foam.PolyMesh(case)
    pts = m.points
    lo, hi = pts.min(0), pts.max(0)
    patches = []
    wedge_like = any(suggest(m, n, t) == "wedge" for n, t, _, _ in read_boundary(case))
    for name, ftype, nf, _ in read_boundary(case):
        s = suggest(m, name, ftype)
        segs = []
        if s not in ("wedge", "empty"):
            for fi in m.patch_faces(name):
                q = pts[m.faces[fi]]
                r = np.hypot(q[:, 0], q[:, 2]) if wedge_like else q[:, 0]
                xy = np.unique(np.round(np.column_stack([r, q[:, 1]]), 9), axis=0)
                if len(xy) >= 2:
                    segs.append(xy[[0, -1]].tolist() if len(xy) == 2 else xy.tolist())
        patches.append(dict(name=name, foamType=ftype, nFaces=nf, suggest=s, segments=segs))
    kinds = {p["suggest"] for p in patches}
    kind = "wedge" if "wedge" in kinds else ("2D" if "empty" in kinds else "3D")
    return dict(nCells=m.nCells, nPoints=len(pts), bounds=[lo.tolist(), hi.tolist()],
                patches=patches, kind=kind)


def build(cfg, src):
    """Template build for a user mesh: no mesh files (the import is the mesh
    step); boundary conditions from cfg['boundaries'], the fill box."""
    m, g, run, mat = cfg["mesh"], cfg["geometry"], cfg["run"], cfg["material"]
    bnd = cfg.get("boundaries") or {}
    walls, outlets, inlets, types = {}, {}, {}, {}
    for p, d in bnd.items():
        d = d if isinstance(d, dict) else {"type": d}
        t = d.get("type")
        if t not in BTYPES:
            raise ValueError(f"boundaries: {p}: type must be one of {', '.join(BTYPES)} "
                             f"(got {t!r})")
        vals = {}
        for e in BTYPES[t][1]:
            v = d.get(e.key)
            if v in (None, ""):
                if e.default == REQUIRED:
                    raise ValueError(f"boundaries: {p}: '{e.key}' is required for a {t} "
                                     f"({e.text}, {e.unit})")
                v = e.default
            vals[e.key] = v
        extra = set(d) - {"type"} - set(vals)
        if extra:
            raise ValueError(f"boundaries: {p}: unknown entr{'y' if len(extra) == 1 else 'ies'} "
                             f"{', '.join(sorted(extra))} for a {t}")
        types[p] = FOAM_TYPE[t]
        if t == "wall":
            walls[p] = float(vals["friction"])
        elif t == "outlet":
            ot = run["openTime"] if vals["openTime"] == AUTO else vals["openTime"]
            outlets[p] = dict(friction=float(vals["friction"]), openTime=float(ot))
        elif t == "inlet":
            inl = dict(speed=float(vals["speed"]), startTime=float(vals["startTime"]),
                       density=float(mat["rho"]) if vals["density"] == AUTO
                       else float(vals["density"]),
                       pressure=float(vals["pressure"]), friction=float(vals["friction"]))
            if vals["endTime"] is not None:
                inl["endTime"] = float(vals["endTime"])
            inlets[p] = inl
    wedge = "wedge" in types.values()
    big = 1e6

    def lim(k, d):
        return d if g.get(k) is None else float(g[k])
    x0, x1 = lim("fillXMin", -big), lim("fillXMax", big)
    y0, y1 = lim("fillYMin", -big), lim("fillYMax", big)
    if x1 <= x0 or y1 <= y0:
        raise ValueError("geometry: the fill box is empty (fillXMax <= fillXMin or "
                         "fillYMax <= fillYMin)")
    fill = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    return dict(
        files={},
        wedge=wedge, theta=None, mesher="import",
        importSpec=dict(file=m["file"], format=m["format"], scale=m["scale"]),
        patchTypes=types,
        bodies={}, walls=walls, outlets=outlets, inlets=inlets,
        wallsNote="walls: mesh patches with their friction (deg), from celSetup boundaries",
        fill=fill, solids=[], volume=None,
        decompose="(1 {n} 1)",
        notes=dict(template="userMesh", fillTop=y1 if y1 < big else None),
        summary=(f"user mesh {os.path.basename(str(m['file']).rstrip('/'))}; boundaries: "
                 + ", ".join(f"{p} {d.get('type') if isinstance(d, dict) else d}"
                             for p, d in bnd.items())),
    )


def finish_import(case):
    """After import: check every patch has a condition, set the patch types,
    measure the wedge angle and the fill volume target. Returns a message."""
    from . import gmsh_mesh
    meta_p = os.path.join(case, "system", "celSetup.json")
    meta = json.load(open(meta_p))
    types = meta.get("patchTypes", {})
    names = [b[0] for b in read_boundary(case)]
    missing = [n for n in names if n not in types]
    unknown = [n for n in types if n not in names]
    if unknown:
        raise ValueError(f"boundaries: no patch {', '.join(unknown)} in the mesh "
                         f"(it has {', '.join(names)})")
    if missing:
        raise ValueError(f"boundaries: give a condition for {', '.join(missing)}")
    gmsh_mesh.set_patch_types(case, types)
    m = foam.PolyMesh(case)
    if meta["wedge"]:
        f = m.patch_faces([n for n, t in types.items() if t == "wedge"][0])
        n = m.Sf[f[0]] / m.magSf[f[0]]
        if abs(n[1]) > 1e-6:
            raise ValueError("mesh: the wedge must be about the y axis")
        meta["theta"] = float(math.atan2(abs(n[0]), abs(n[2])))
    elif "empty" not in types.values():
        raise ValueError("mesh: 3D meshes are not supported here yet; give the front and "
                         "back patches of a 2D mesh type empty, or use a wedge")
    cells, polys, _ = m.front_polygons()
    if len(cells) != m.nCells:
        raise ValueError(f"mesh: not one cell thick ({len(cells)} front faces, "
                         f"{m.nCells} cells)")
    if max(len(q) for q in polys) > 4:
        raise ValueError("mesh: cells must be quadrilaterals or triangles in the section")
    json.dump(meta, open(meta_p, "w"), indent=1)
    return f"imported {m.nCells} cells, {len(names)} patches"
