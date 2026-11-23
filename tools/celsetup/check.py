"""celCheck: a mesh check before a run (WP S.3).

Reads a case (constant/polyMesh, materialProperties, controlDict and, when
celSetup wrote the case, system/celSetup.json) and reports:

- the stable time step, computed as celFoam does (eulerianSolid::cellSize
  and stableDeltaT), the cell that sets it and how isolated it is;
- cell sizes where they matter: across the outlet, along the walls, against
  the grain size d (shear bands 10 to 20 d wide) and the nonlocal length;
- non-orthogonality and the size jump between neighbours;
- the estimated run time on 1 to 8 cores;
- optionally a picture coloured by cell size (needs matplotlib) and a JSON
  dump for the GUI.

Rules are those of the manual's meshing guide (decisions 0038, 0105, 0111,
0115). Each is PASS, WARN or FAIL with the numbers behind it.
"""

import json
import math
import os
import re

import numpy as np

from . import foam

# Run-time cost: Couto run 1 of decision 0115 (12,336 cells, dt 5.85e-6 s,
# 2 s simulated, 3.7 h on 4 cores, T0 rules + fullBlend + nonlocal):
# 3.7 h x 4 cores / (12,336 x 341,880 cell-steps) = 12.6 us per cell-step-core.
US_PER_CELL_STEP = 12.6
CORE_EFFICIENCY = {1: 1.0, 2: 0.95, 4: 0.9, 8: 0.75}


def read_dict_entries(path):
    """Top-level 'key value;' entries of an OpenFOAM dictionary (simple)."""
    if not os.path.exists(path):
        return {}
    txt = re.sub(r"/\*.*?\*/", "", open(path).read(), flags=re.S)
    txt = re.sub(r"//[^\n]*", "", txt)
    out, depth = {}, 0
    for m in re.finditer(r"[{}]|(\w+)\s+([^;{}]+);", txt):
        s = m.group(0)
        if s == "{":
            depth += 1
        elif s == "}":
            depth -= 1
        elif depth == 0:
            out.setdefault(m.group(1), m.group(2).strip())
    return out


def subdict_names(path, name):
    """Entry names inside a top-level sub-dictionary (walls, outlets)."""
    if not os.path.exists(path):
        return []
    txt = re.sub(r"//[^\n]*", "", open(path).read())
    m = re.search(r"\b" + name + r"\s*\{", txt)
    if not m:
        return []
    i, depth, names = m.end(), 1, []
    while i < len(txt) and depth > 0:
        mm = re.compile(r"(\w+)\s*\{|\}").search(txt, i)
        if not mm:
            break
        if mm.group(0) == "}":
            depth -= 1
        else:
            if depth == 1:
                names.append(mm.group(1))
            depth += 1
        i = mm.end()
    return names


def _f(d, k, default=None):
    try:
        return float(d[k])
    except (KeyError, ValueError):
        return default


class Mesh(foam.PolyMesh):
    """PolyMesh with celFoam's stable cell length."""

    def boundary_kind(self):
        kinds = np.full(len(self.faces), "internal", dtype=object)
        for n, p in self.patches.items():
            kinds[p["startFace"]:p["startFace"] + p["nFaces"]] = p["type"]
        return kinds

    def angled(self):
        kinds = self.boundary_kind()
        S = self.Sf
        mag = np.maximum(self.magSf, 1e-300)
        off = np.abs(S).max(1) < (1 - 1e-6) * mag
        use = (kinds != "empty") & (kinds != "wedge") & (self.magSf > 0)
        return bool((off & use).any())

    def stable_length(self):
        """celFoam's cell length for the time step (eulerianSolid::cellSize)."""
        nc, ni = self.nCells, self.nInternal
        kinds = self.boundary_kind()
        wedge = bool(self.patches_of_type("wedge"))
        angled = self.angled()
        own, nei = self.owner, self.neighbour
        V = self.V
        skip = (kinds == "empty") | (kinds == "wedge")
        # V over the largest face (angled)
        Smax = np.zeros(nc)
        np.maximum.at(Smax, own[:ni], self.magSf[:ni])
        np.maximum.at(Smax, nei, self.magSf[:ni])
        b = np.arange(ni, len(self.faces))
        b = b[~skip[b]]
        np.maximum.at(Smax, own[b], self.magSf[b])
        angled_len = V / np.maximum(Smax, 1e-300)
        if not wedge:
            return (angled_len if angled else np.cbrt(V)), dict(wedge=False, angled=angled)
        A = np.zeros((nc, 3))
        np.add.at(A, own[:ni], 0.5 * np.abs(self.Sf[:ni]))
        np.add.at(A, nei, 0.5 * np.abs(self.Sf[:ni]))
        notempty = np.arange(ni, len(self.faces))
        notempty = notempty[(kinds[notempty] != "empty") & (kinds[notempty] != "wedge")]
        np.add.at(A, own[notempty], 0.5 * np.abs(self.Sf[notempty]))
        Sw = np.zeros((nc, 3))
        wf = np.arange(ni, len(self.faces))
        wf = wf[kinds[wf] == "wedge"]
        np.add.at(Sw, own[wf], self.Sf[wf])
        if angled:
            d = angled_len.copy()
        else:
            # geometric directions: those with extent (empty direction left out)
            span = self.points.max(0) - self.points.min(0)
            d = np.full(nc, np.inf)
            for a in range(3):
                if span[a] > 1e-12:
                    ok = A[:, a] > 1e-300
                    d[ok] = np.minimum(d[ok], V[ok] / A[ok, a])
        sw = np.linalg.norm(Sw, axis=1)
        r = np.where(sw > 1e-300, V / np.maximum(sw, 1e-300), np.inf)
        hoop = np.isfinite(r)
        d_in = d.copy()
        d[hoop] = 1 / np.sqrt(1 / d[hoop] ** 2 + 1 / r[hoop] ** 2)
        return d, dict(wedge=True, angled=angled, inPlane=d_in, r=r)


def inplane_sizes(polys):
    """Smaller and larger side span of each (r, y) cell polygon."""
    small, large = [], []
    for q in polys:
        if len(q) == 3:
            q = np.vstack([q, q[-1]])
        e = [np.linalg.norm(q[(i + 1) % 4] - q[i]) for i in range(4)]
        a, b = max(e[0], e[2]), max(e[1], e[3])
        small.append(min(a, b))
        large.append(max(a, b))
    return np.array(small), np.array(large)


def non_orthogonality(m):
    ni = m.nInternal
    dvec = m.C[m.neighbour] - m.C[m.owner[:ni]]
    n = m.Sf[:ni] / np.maximum(m.magSf[:ni], 1e-300)[:, None]
    c = np.einsum("ij,ij->i", dvec, n) / np.maximum(np.linalg.norm(dvec, axis=1), 1e-300)
    ok = m.magSf[:ni] > 0
    return np.degrees(np.arccos(np.clip(c[ok], -1, 1)))


def run_check(case):
    case = os.path.abspath(case)
    mp = read_dict_entries(os.path.join(case, "constant", "materialProperties"))
    cd = read_dict_entries(os.path.join(case, "system", "controlDict"))
    meta_p = os.path.join(case, "system", "celSetup.json")
    meta = json.load(open(meta_p)) if os.path.exists(meta_p) else {}
    m = Mesh(case)
    res = dict(case=case, nCells=m.nCells, checks=[])

    def check(name, status, text):
        res["checks"].append(dict(name=name, status=status, text=text))

    # material and time step
    rho = _f(mp, "rho")
    E, nu = _f(mp, "E"), _f(mp, "nu", 0.3)
    if rho is None or E is None:
        raise ValueError("materialProperties needs rho and E for the time step")
    M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    c = math.sqrt(M / rho)
    maxCo = _f(cd, "maxCo", 0.4)
    L, info = m.stable_length()
    dts = maxCo * L / c
    k = int(np.argmin(dts))
    dt = float(dts[k])
    res.update(dt=dt, waveSpeed=c, maxCo=maxCo, wedge=info["wedge"], angled=info["angled"])
    lim = dict(cell=k, length=float(L[k]), centre=m.C[k].tolist())
    if info["wedge"]:
        lim.update(r=float(np.hypot(m.C[k, 0], m.C[k, 2])), y=float(m.C[k, 1]),
                   inPlane=float(info["inPlane"][k]), hoopRadius=float(info["r"][k]))
    res["limit"] = lim
    p1 = float(np.percentile(L, 1))
    res["dt1pct"] = maxCo * p1 / c
    n_near = int((L < 1.1 * L[k]).sum())
    where = (f"at r {lim['r'] * 1000:.1f} mm, y {lim['y'] * 1000:.1f} mm" if info["wedge"]
             else "at ({:.4f} {:.4f} {:.4f}) m".format(*lim["centre"]))
    why = ""
    if info["wedge"]:
        why = (f"; in-plane height {lim['inPlane'] * 1000:.2f} mm, hoop length (distance from the axis) "
               f"{lim['hoopRadius'] * 1000:.2f} mm")
    gain = res["dt1pct"] / dt
    status = "PASS" if gain < 1.25 else "WARN"
    check("time step", status,
          f"dt {dt:.4g} s (maxCo {maxCo:g}, P-wave {c:.0f} m/s), set by a {L[k] * 1000:.2f} mm cell "
          f"{where}{why}. {n_near} cells within 10 % of it. Without the smallest 1 % of cells "
          f"the step would be {gain:.2f} times longer")

    # in-plane sizes (2D and wedge) or cube roots (3D)
    cells, polys, wedge = m.front_polygons()
    if len(cells):
        small, large = inplane_sizes(polys)
        size = np.zeros(m.nCells)
        size[cells] = small
        rr = np.array([q[:, 0].mean() for q in polys])
        yy = np.array([q[:, 1].mean() for q in polys])
        res["polys"] = dict(cells=cells.tolist(), small=small.tolist())
    else:
        size = np.cbrt(m.V)
        cells = np.arange(m.nCells)
        rr, yy = m.C[:, 0], m.C[:, 1]
    # cells whose centre lies inside an STL wall (celSetup cases) hold no
    # material: left out of the size-jump rule
    solid = np.zeros(m.nCells, bool)
    if meta.get("solids") and len(cells):
        from .case import _inside
        P = np.column_stack([rr, yy])
        for sp in meta["solids"]:
            solid[cells[_inside(P, np.array(sp))]] = True
    res["size"] = dict(min=float(size[cells].min()), median=float(np.median(size[cells])),
                       max=float(size[cells].max()))
    d = _f(mp, "d")
    A = _f(mp, "nonlocalA")

    # outlet
    outlets = subdict_names(os.path.join(case, "constant", "materialProperties"), "outlets")
    for o in outlets:
        if o not in m.patches:
            continue
        f = m.patch_faces(o)
        oc = m.owner[f]
        osz = size[oc]
        if wedge:
            r = np.hypot(m.Cf[f, 0], m.Cf[f, 2])
            half = float(r.max() + 0.5 * osz[np.argmax(r)])
        else:
            x = m.Cf[f, 0]
            half = 0.5 * float(x.max() - x.min() + osz.mean())
        n_across = half / float(osz.mean())
        st = "PASS" if n_across >= 35 else ("WARN" if n_across >= 20 else "FAIL")
        what = "radius" if wedge else "half-width"
        check(f"outlet '{o}'", st,
              f"cells {osz.min() * 1000:.1f} to {osz.max() * 1000:.1f} mm, {n_across:.0f} across the "
              f"outlet {what} ({half * 1000:.0f} mm). Rule: about 40 (3 to 5 mm on the 0.35 m Couto "
              f"outlet); the outlet cells set the discharge rate")
        res.setdefault("outlets", {})[o] = dict(min=float(osz.min()), max=float(osz.max()),
                                                across=n_across)

    # walls: cells near the wall line (celSetup case) or next to wall patches
    wall_cells = None
    if meta.get("wall") and len(cells):
        W = np.array(meta["wall"])
        P = np.column_stack([rr, yy])
        dist = np.full(len(P), np.inf)
        for a, b in zip(W[:-1], W[1:]):
            ab = b - a
            t = np.clip(((P - a) @ ab) / max(ab @ ab, 1e-300), 0, 1)
            dist = np.minimum(dist, np.linalg.norm(P - (a + t[:, None] * ab), axis=1))
        sel = (dist < 1.5 * size[cells]) & (yy > 0) & (yy < meta.get("fillTop", np.inf))
        wall_cells = cells[sel]
    else:
        names = [n for n in subdict_names(os.path.join(case, "constant", "materialProperties"),
                                          "walls") if n in m.patches]
        if names:
            wall_cells = np.unique(np.concatenate([m.owner[m.patch_faces(n)] for n in names]))
    if wall_cells is not None and len(wall_cells):
        ws = size[wall_cells]
        txt = f"cells along the wall {ws.min() * 1000:.1f} to {ws.max() * 1000:.1f} mm " \
              f"(median {np.median(ws) * 1000:.1f})"
        st = "PASS"
        if d:
            band = 10 * d
            n_band = band / float(np.median(ws))
            txt += f"; {n_band:.1f} cells across a 10 d = {band * 1000:.0f} mm shear band"
            st = "PASS" if n_band >= 2 else ("WARN" if n_band >= 1 else "FAIL")
        check("walls", st, txt + ". Rule: finer columns at walls, a band of 10 to 20 d "
              "needs a few cells")
        res["wallCells"] = wall_cells.tolist()

    # grain size and nonlocal length
    if d:
        core = float(np.median(size[cells]))
        check("grain size", "INFO",
              f"d {d * 1000:g} mm: smallest cell {size[cells].min() / d:.2f} d, median {core / d:.1f} d, "
              f"largest {size[cells].max() / d:.1f} d. A local friction law puts a shear band in 1 to 4 "
              f"cells (mesh dependent); cells near d converge it with the nonlocal law")
        if A:
            xi = 3 * d
            n_xi = xi / core
            st = "PASS" if n_xi >= 1 else "WARN"
            check("nonlocal length", st,
                  f"nonlocalA {A:g}: cooperative length about 3 d = {xi * 1000:.0f} mm, "
                  f"{n_xi:.2f} median cells per length (decision 0115: 12 mm rows still did not "
                  f"resolve it; cells near d do)")

    # mesh quality
    no = non_orthogonality(m)
    st = "PASS" if no.max() < 65 else ("WARN" if no.max() < 70 else "FAIL")
    check("non-orthogonality", st, f"max {no.max():.1f} deg, mean {no.mean():.1f} deg "
          "(checkMesh limit 70)")
    ni = m.nInternal
    s_o, s_n = size[m.owner[:ni]], size[m.neighbour]
    ok = (s_o > 0) & (s_n > 0) & ~solid[m.owner[:ni]] & ~solid[m.neighbour]
    jump = np.maximum(s_o[ok] / s_n[ok], s_n[ok] / s_o[ok])
    st = "PASS" if jump.max() <= 1.3 else "WARN"
    fj = np.arange(ni)[ok][np.argmax(jump)]
    cj = m.Cf[fj]
    rj = np.hypot(cj[0], cj[2]) if wedge else cj[0]
    check("size jump", st, f"largest neighbour size ratio {jump.max():.2f} at "
          f"{'r' if wedge else 'x'} {rj * 1000:.0f} mm, y {cj[1] * 1000:.0f} mm; "
          f"{int((jump > 1.3).sum())} faces above 1.3 (rule: about 1.2 or less; cells inside "
          "an STL wall left out)")
    res["nonOrth"] = dict(max=float(no.max()), mean=float(no.mean()))
    res["jump"] = float(jump.max())

    # run time
    T = _f(cd, "endTime", meta.get("endTime", 1.0))
    t0 = _f(cd, "startTime", 0.0) or 0.0
    steps = (T - t0) / dt
    hours = {n: m.nCells * steps * US_PER_CELL_STEP * 1e-6 / (n * e) / 3600
             for n, e in CORE_EFFICIENCY.items()}
    res.update(endTime=T, steps=steps, hours=hours)
    check("run time", "INFO",
          f"{T - t0:g} s simulated = {steps:,.0f} steps x {m.nCells:,} cells: about "
          + ", ".join((f"{h:.1f} h" if h >= 1 else f"{h * 60:.0f} min") + f" on {n}"
                      for n, h in hours.items())
          + " cores (rough, from the Couto runs)")
    if meta.get("mass"):
        res["mass"] = meta["mass"]
    if meta.get("solids"):
        res["solids"] = meta["solids"]
    return res, m


def report(res):
    lines = [f"celCheck: {res['case']}", f"  {res['nCells']:,} cells, "
             + ("axisymmetric wedge" if res["wedge"] else "plane/3D")
             + (", angled cells" if res["angled"] else "")]
    for c in res["checks"]:
        lines.append(f"  {c['status']:<5} {c['name']}: {c['text']}")
    worst = "FAIL" if any(c["status"] == "FAIL" for c in res["checks"]) else (
        "WARN" if any(c["status"] == "WARN" for c in res["checks"]) else "PASS")
    lines.append(f"  verdict: {worst}")
    return "\n".join(lines), worst


def picture(res, m, png, views=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    cells, polys, wedge = m.front_polygons()
    small, _ = inplane_sizes(polys)
    val = small * 1000
    lim = res["limit"]
    allp = np.vstack(polys)
    lo, hi = allp.min(0), allp.max(0)

    def box(cx, cy, w):
        x0 = min(max(cx - w / 2, lo[0] - 0.02 * w), hi[0] - w)
        y0 = min(max(cy - w / 2, lo[1] - 0.02 * w), hi[1] - w)
        return (x0, x0 + w), (y0, y0 + w)
    if views is None:
        views = [((lo[0], hi[0]), (lo[1], hi[1]), "Whole domain")]
        w = 0.15 * max(hi - lo)
        ref = (lim.get("r", lim["centre"][0]), lim.get("y", lim["centre"][1]))
        views.append(box(ref[0] + 0.3 * w, ref[1] + 0.3 * w, w)
                     + ("Around the cell that sets the time step",))
        if res.get("wallCells"):
            pos = {c: i for i, c in enumerate(cells)}
            wc = [polys[pos[c]].mean(0) for c in res["wallCells"] if c in pos]
            if wc:
                wc = np.array(wc)
                yt = lo[1] + 0.2 * (hi[1] - lo[1])
                mid = wc[np.argmin(np.abs(wc[:, 1] - yt))]
                views.append(box(mid[0], mid[1], w) + ("Along the wall",))
    fig, axs = plt.subplots(1, len(views), figsize=(4 + 5 * (len(views) - 1), 8))
    axs = np.atleast_1d(axs)
    vmax = float(np.percentile(val, 98))
    for ax, (xl, yl, title) in zip(axs, views):
        pc = PolyCollection(polys, array=val, cmap="Spectral", clim=(val.min(), vmax),
                            edgecolors="#55606b", linewidths=0.15)
        ax.add_collection(pc)
        for sp in res.get("solids", []):
            sp = np.array(sp + sp[:1])
            ax.plot(sp[:, 0], sp[:, 1], color="#000000", lw=1.0)
        ax.plot([lim.get("r", lim["centre"][0])], [lim.get("y", lim["centre"][1])], marker="o",
                mfc="none", mec="#000000", ms=10, mew=1.5, label="sets the time step")
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlabel("r (m)" if wedge else "x (m)")
    axs[0].set_ylabel("y (m)")
    axs[0].legend(loc="upper left", fontsize=8)
    fig.colorbar(pc, ax=list(axs), shrink=0.6, label="cell size, smaller side (mm)")
    fig.suptitle(f"{os.path.basename(res['case'])}: {res['nCells']:,} cells, dt {res['dt']:.3g} s, "
                 f"about {res['hours'][4]:.1f} h on 4 cores for {res['endTime']:g} s")
    fig.savefig(png, dpi=120)
    plt.close(fig)
