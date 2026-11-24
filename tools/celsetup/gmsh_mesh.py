"""Gmsh mesher for the templates (WP G.1): the same wall-aligned block
layout as the blockMesh dictionary, built in Gmsh as structured
(transfinite, recombined) quadrilateral surfaces in the (r, y) plane, then
revolved into a one-cell wedge (or extruded into a one-cell plane slice),
written as a version 2.2 .msh file and converted with OpenFOAM's
gmshToFoam. Gmsh is GPL and installs with `apt install gmsh` or
`pip install gmsh`; celSetup uses its Python module when it is there, else
the `gmsh` command on a .geo file it writes.

Differences from blockMesh: Gmsh places interior nodes by transfinite
interpolation of the four sides, blockMesh by its own interpolation of the
graded edges, so the meshes agree in layout and cell counts and differ a
little inside the blocks (celCheck compares them).
"""

import os
import re
import shutil
import subprocess

from . import foam


def available():
    try:
        import gmsh  # noqa: F401
        return "python"
    except (ImportError, OSError):
        pass
    return "command" if shutil.which("gmsh") else None


def _progression(n, ratio):
    """Gmsh Progression coefficient: each cell r times the previous one, so
    that last/first = ratio over n cells."""
    return 1.0 if n <= 1 or ratio == 1 else ratio ** (1.0 / (n - 1))


def geo_script(L, G, depth=None):
    """A .geo script (built-in kernel) for the conical hopper layout L
    (conical.layout). Curves split where the grading changes, so each side
    of each surface has a single progression."""
    lines = ["// written by celSetup (WP G.1)", "SetFactory(\"Built-in\");",
             "Mesh.MshFileVersion = 2.2;", "Mesh.RecombineAll = 1;", "Mesh.Algorithm = 8;"]
    pid = [0]

    def pt(r, y):
        pid[0] += 1
        lines.append(f"Point({pid[0]}) = {{{r!r}, {y!r}, 0}};")
        return pid[0]
    cid = [0]

    def line(a, b, n, coef=1.0):
        cid[0] += 1
        lines.append(f"Line({cid[0]}) = {{{a}, {b}}};")
        lines.append(f"Transfinite Curve {{{cid[0]}}} = {n + 1} Using Progression {coef!r};")
        return cid[0]

    Y0, y_t, Hm, R, D, TS = G.Y0, L["y_t"], L["Hm"], L["R"], L["D"], G.TS
    M, O = L["M"], L["O"]
    ncols, nb = L["ncols"], L["nb"]
    # row breaks: hopper block one segment, cylinder one or two
    ys = [Y0, y_t]
    y = y_t
    for s in L["rows2"]:
        y += s[0]
        ys.append(y)
    ys[-1] = Hm
    rows = [L["rows1"][0]] + list(L["rows2"])
    # points: axis, shell middle M, outer O at each break
    ax = [pt(0.0, yy) for yy in ys]
    m_r = [float(M[0][0]), float(M[-1][0])] + [R + 0.5 * TS] * (len(ys) - 2)
    o_r = [float(O[0][0]), float(O[-1][0])] + [R + D] * (len(ys) - 2)
    mp = [pt(m_r[i], ys[i]) for i in range(len(ys))]
    op = [pt(o_r[i], ys[i]) for i in range(len(ys))]
    cw = _progression(ncols, L["wallRatio"])
    # horizontal curves (axis -> M, M -> O) at each break
    hA = [line(ax[i], mp[i], ncols, cw) for i in range(len(ys))]
    hB = [line(mp[i], op[i], nb) for i in range(len(ys))]
    # vertical curves per segment
    vA, vM, vO = [], [], []
    for i, (ln, n, ratio) in enumerate(rows):
        c = _progression(n, ratio)
        vA.append(line(ax[i], ax[i + 1], n, c))
        if i == 0:
            # the curved lower part of M and O: splines through the paths
            mids = [pt(float(r), float(yy)) for r, yy in M[1:-1]]
            cid[0] += 1
            lines.append(f"Spline({cid[0]}) = {{{', '.join(map(str, [mp[0]] + mids + [mp[1]]))}}};")
            lines.append(f"Transfinite Curve {{{cid[0]}}} = {n + 1} Using Progression {c!r};")
            vM.append(cid[0])
            mids = [pt(float(r), float(yy)) for r, yy in O[1:-1]]
            cid[0] += 1
            lines.append(f"BSpline({cid[0]}) = {{{', '.join(map(str, [op[0]] + mids + [op[1]]))}}};")
            lines.append(f"Transfinite Curve {{{cid[0]}}} = {n + 1} Using Progression {c!r};")
            vO.append(cid[0])
        else:
            vM.append(line(mp[i], mp[i + 1], n, c))
            vO.append(line(op[i], op[i + 1], n, c))
    # surfaces: inner (axis-M) and buffer (M-O) per row segment
    sid = 0
    inner, buf = [], []
    for i in range(len(rows)):
        for kind, (lo, hi, left, right) in (("in", (hA[i], hA[i + 1], vA[i], vM[i])),
                                             ("buf", (hB[i], hB[i + 1], vM[i], vO[i]))):
            sid += 1
            lines.append(f"Curve Loop({sid}) = {{{lo}, {right}, -{hi}, -{left}}};")
            lines.append(f"Plane Surface({sid}) = {{{sid}}};")
            lines.append(f"Transfinite Surface {{{sid}}};")
            lines.append(f"Recombine Surface {{{sid}}};")
            (inner if kind == "in" else buf).append(sid)
    surfs = inner + buf
    sl = ", ".join(map(str, surfs))
    if depth:
        lines.append(f"Translate {{0, 0, {-0.5 * depth!r}}} {{ Surface{{{sl}}}; }}")
        lines.append(f"ex[] = Extrude {{0, 0, {depth!r}}} {{ Surface{{{sl}}}; Layers{{1}}; Recombine; }};")
    else:
        th = G.theta
        lines.append(f"Rotate {{{{0, 1, 0}}, {{0, 0, 0}}, {th!r}}} {{ Surface{{{sl}}}; }}")
        lines.append(f"ex[] = Extrude {{{{0, 1, 0}}, {{0, 0, 0}}, {-2 * th!r}}} "
                     f"{{ Surface{{{sl}}}; Layers{{1}}; Recombine; }};")
    # Extrude returns, per surface: top, volume, then the side surfaces in
    # the order of the curve loop (lo, right, hi, left). Revolving drops the
    # inner surfaces' left side (on the axis, it collapses), the last one.
    lines.append("// physical groups from the extrusion result")
    start, k = {}, 0
    for s in surfs:
        start[s] = k
        k += 5 if (s in inner and not depth) else 6

    def side(s, j):
        return f"ex[{start[s] + 2 + j}]"
    outlet = [side(inner[0], 0)]
    bottom = [side(buf[0], 0)]
    outer = [side(b, 1) for b in buf]
    top = [side(inner[-1], 2), side(buf[-1], 2)]
    axis = [side(s, 3) for s in inner]
    back = [f"ex[{start[s]}]" for s in surfs]
    lines.append(f"Physical Surface(\"outlet\") = {{{', '.join(outlet)}}};")
    lines.append(f"Physical Surface(\"bottom\") = {{{', '.join(bottom)}}};")
    lines.append(f"Physical Surface(\"outer\") = {{{', '.join(outer)}}};")
    lines.append(f"Physical Surface(\"top\") = {{{', '.join(top)}}};")
    if depth:
        lines.append(f"Physical Surface(\"axis\") = {{{', '.join(axis)}}};")
    lines.append(f"Physical Surface(\"front\") = {{{sl}}};")
    lines.append(f"Physical Surface(\"back\") = {{{', '.join(back)}}};")
    lines.append(f"Physical Volume(\"internal\") = {{{', '.join(f'ex[{start[s] + 1}]' for s in surfs)}}};")
    return "\n".join(lines) + "\n"


def mesh_case(case, types):
    """Mesh system/mesh.geo with Gmsh (3D), convert with gmshToFoam and set
    the patch types {patch: type}. Returns (ok, message)."""
    how = available()
    if how is None:
        return False, "Gmsh not found: apt install gmsh, or pip install gmsh"
    if shutil.which("gmshToFoam") is None:
        return False, "gmshToFoam not found: source OpenFOAM's bashrc"
    gp = os.path.join(case, "system", "mesh.geo")
    msh = os.path.join(case, "mesh.msh")
    with open(os.path.join(case, "log.gmsh"), "w") as log:
        if how == "python":
            import gmsh
            # not interruptible: the GUI calls this from a worker thread
            gmsh.initialize(["gmsh"], interruptible=False)
            try:
                gmsh.option.setNumber("General.Terminal", 0)
                gmsh.open(gp)
                gmsh.model.mesh.generate(3)
                gmsh.write(msh)
            except Exception as e:  # gmsh raises plain Exceptions
                log.write(str(e) + "\n")
                return False, f"Gmsh failed: {e}"
            finally:
                gmsh.finalize()
        else:
            r = subprocess.run(["gmsh", "-3", gp, "-o", msh, "-format", "msh22"],
                               stdout=log, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                return False, "Gmsh failed, see " + os.path.join(case, "log.gmsh")
    with open(os.path.join(case, "log.gmshToFoam"), "w") as f:
        r = subprocess.run(["gmshToFoam", "mesh.msh"], cwd=case, stdout=f,
                           stderr=subprocess.STDOUT,
                           env=dict(os.environ, PWD=os.path.abspath(case)))
    if r.returncode != 0:
        return False, "gmshToFoam failed, see " + os.path.join(case, "log.gmshToFoam")
    set_patch_types(case, types)
    os.remove(msh)
    return True, "Gmsh mesh written"


def set_patch_types(case, types):
    """Change patch types in constant/polyMesh/boundary (gmshToFoam writes
    every patch as 'patch'); drop the physicalType entries."""
    p = os.path.join(case, "constant", "polyMesh", "boundary")
    txt = open(p).read()
    for name, ty in types.items():
        txt = re.sub(r"(\n\s*" + name + r"\s*\{[^}]*?type\s+)\w+;", r"\g<1>" + ty + ";", txt)
    txt = re.sub(r"\n\s*physicalType\s+\w+;", "", txt)
    txt = re.sub(r"\n\s*inGroups\s+[^;]*;", "", txt)
    open(p, "w").write(txt)
