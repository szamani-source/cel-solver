"""OpenFOAM file helpers for celSetup and celCheck: dictionary headers,
writing files, and an ASCII polyMesh reader (numpy only)."""

import os
import re

import numpy as np

BANNER = """/*--------------------------------*- C++ -*----------------------------------*\\
  Written by celSetup ({src}). Edit the YAML file and run celSetup again
  rather than editing this file by hand.
\\*---------------------------------------------------------------------------*/
"""

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def header(cls, obj, src="celSetup"):
    return BANNER.format(src=src) + HEADER.format(cls=cls, obj=obj)


def write(path, text, mode=None):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
    if mode is not None:
        os.chmod(path, mode)


def fmt(v):
    """A Python value as an OpenFOAM dictionary value."""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "(" + " ".join(fmt(x) for x in v) + ")"
    return str(v)


# --- polyMesh reader ------------------------------------------------------

def _body(path):
    with open(path) as f:
        txt = f.read()
    if "format      binary" in txt[:2000] or "format binary" in txt[:2000]:
        raise ValueError(f"{path} is binary; set writeFormat ascii in system/controlDict "
                         "and run blockMesh again")
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    txt = re.sub(r"//[^\n]*", "", txt)
    i = txt.find("FoamFile")
    if i >= 0:
        j = txt.find("}", i)
        txt = txt[j + 1:]
    return txt


def read_points(path):
    txt = _body(path)
    n = int(re.search(r"(\d+)\s*\(", txt).group(1))
    nums = re.findall(r"\(([^()]*)\)", txt[txt.find("(") + 1:])
    a = np.array([[float(x) for x in s.split()] for s in nums[:n]])
    return a


def read_faces(path):
    txt = _body(path)
    m = re.search(r"(\d+)\s*\(", txt)
    n = int(m.group(1))
    out = [np.array(s.split(), dtype=int)
           for s in re.findall(r"\d+\s*\(([^()]*)\)", txt[m.end():])[:n]]
    return out


def read_labels(path):
    txt = _body(path)
    m = re.search(r"(\d+)\s*\(", txt)
    n = int(m.group(1))
    return np.array(txt[m.end():txt.rfind(")")].split()[:n], dtype=int)


def read_boundary(path):
    txt = _body(path)
    out = {}
    for name, body in re.findall(r"(\w+)\s*\{([^{}]*)\}", txt):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", body))
        out[name] = dict(type=d.get("type", "patch"), nFaces=int(d.get("nFaces", 0)),
                         startFace=int(d.get("startFace", 0)))
    return out


class PolyMesh:
    """Points, faces, owner, neighbour and patches of constant/polyMesh,
    with face areas, face centres, cell volumes and cell centres."""

    def __init__(self, case):
        m = os.path.join(case, "constant", "polyMesh")
        self.points = read_points(f"{m}/points")
        self.faces = read_faces(f"{m}/faces")
        self.owner = read_labels(f"{m}/owner")
        self.neighbour = read_labels(f"{m}/neighbour")
        self.patches = read_boundary(f"{m}/boundary")
        self.nCells = int(self.owner.max()) + 1
        self.nInternal = len(self.neighbour)
        self._geometry()

    def _geometry(self):
        P = self.points
        nf = len(self.faces)
        Sf = np.zeros((nf, 3))
        Cf = np.zeros((nf, 3))
        sizes = np.array([len(f) for f in self.faces])
        for n in np.unique(sizes):
            idx = np.nonzero(sizes == n)[0]
            q = P[np.array([self.faces[i] for i in idx])]          # (k, n, 3)
            c0 = q.mean(1)
            q1 = np.roll(q, -1, axis=1)
            # triangle fan about the point average (as OpenFOAM does)
            a = np.cross(q - c0[:, None], q1 - c0[:, None])
            tc = (q + q1 + c0[:, None]) / 3.0
            ma = np.linalg.norm(a, axis=2)
            Sf[idx] = 0.5 * a.sum(1)
            tot = ma.sum(1)
            Cf[idx] = np.where(tot[:, None] > 0,
                               (tc * ma[:, :, None]).sum(1) / np.where(tot > 0, tot, 1)[:, None], c0)
        self.Sf, self.Cf = Sf, Cf
        self.magSf = np.linalg.norm(Sf, axis=1)
        nc = self.nCells
        # cell volume and centre: divergence theorem about an estimate
        own, nei = self.owner, self.neighbour
        cEst = np.zeros((nc, 3))
        cnt = np.zeros(nc)
        np.add.at(cEst, own, Cf)
        np.add.at(cnt, own, 1)
        np.add.at(cEst, nei, Cf[:len(nei)])
        np.add.at(cnt, nei, 1)
        cEst /= cnt[:, None]
        V = np.zeros(nc)
        C = np.zeros((nc, 3))

        def acc(cells, sign, fidx):
            pyr = sign * np.einsum("ij,ij->i", Sf[fidx], Cf[fidx] - cEst[cells]) / 3.0
            np.add.at(V, cells, pyr)
            np.add.at(C, cells, pyr[:, None] * (0.75 * Cf[fidx] + 0.25 * cEst[cells]))
        acc(own, 1.0, np.arange(nf))
        acc(nei, -1.0, np.arange(len(nei)))
        self.V = V
        self.C = C / np.where(np.abs(V) > 0, V, 1.0)[:, None]

    def patch_faces(self, name):
        p = self.patches[name]
        return np.arange(p["startFace"], p["startFace"] + p["nFaces"])

    def patches_of_type(self, *types):
        return [n for n, p in self.patches.items() if p["type"] in types]

    def front_polygons(self):
        """2D or wedge meshes (one cell thick): for each cell, its polygon
        in the (r, y) plane (r = distance from the y axis on a wedge, x on a
        plane mesh), from the front patch of type wedge or empty with
        z < 0. Returns (cells, list of (n, 2) arrays), wedge flag."""
        wedge = bool(self.patches_of_type("wedge"))
        names = self.patches_of_type("wedge") or self.patches_of_type("empty")
        cells, polys = [], []
        for n in names:
            for fi in self.patch_faces(n):
                if self.Sf[fi][2] >= 0:      # keep the z-min side
                    continue
                q = self.points[self.faces[fi]]
                r = np.hypot(q[:, 0], q[:, 2]) if wedge else q[:, 0]
                cells.append(self.owner[fi])
                polys.append(np.column_stack([r, q[:, 1]]))
        return np.array(cells, dtype=int), polys, wedge
