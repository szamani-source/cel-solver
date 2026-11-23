"""celSetup GUI (WP S.4): a local web page, served by Python's own http.server.

    celSetup --gui [case.yaml] [--port 8765] [--no-browser]

The page (gui/index.html, app.js, style.css; no internet needed) shows a
form built from the schema (every entry with its default, unit and
description), draws the geometry at once as values change, and builds the
real mesh in the background (blockMesh or Gmsh, the fill and celCheck in a
scratch folder) to show it coloured by cell size or material, with the
checker's verdict; Compare meshes with the other mesher too and draws its
cell edges over the first (WP G.2). Save writes the YAML file and the case; Run starts ./Allrun and
shows its progress. Only 127.0.0.1 is served.
"""

import json
import os
import shutil
import subprocess
import tempfile
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

from . import case, check, foam, miniyaml, schema, usermesh

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC = os.path.join(HERE, "gui")
EXAMPLES = os.path.join(HERE, "examples")
TYPES = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
         ".css": "text/css; charset=utf-8"}


class State:
    def __init__(self, yaml_path=None):
        self.yaml_path = os.path.abspath(yaml_path) if yaml_path else None
        self.scratch = tempfile.mkdtemp(prefix="celSetupGui-")
        self.lock = threading.Lock()
        self.runs = {}
        self.user = None            # the loaded user mesh (G.3): description


def initial_config(path):
    if path and os.path.exists(path):
        return miniyaml.load(path)
    return miniyaml.load(os.path.join(EXAMPLES, "coutoHopper.yaml"))


def clean(cfg):
    """Form values arrive as strings: numbers and words as the YAML reader
    reads them; empty means not given."""
    out = {}
    for k, v in cfg.items():
        if k == "boundaries" and isinstance(v, dict):
            out[k] = {p: {kk: miniyaml._scalar(vv) if isinstance(vv, str) else vv
                          for kk, vv in d.items() if vv not in ("", None)}
                      for p, d in v.items()}
        elif isinstance(v, dict):
            out[k] = {kk: miniyaml._scalar(vv) if isinstance(vv, str) else vv
                      for kk, vv in v.items() if vv not in ("", None)}
        else:
            out[k] = v
    return out


def abs_mesh_file(state, cfg):
    """A relative user-mesh path is taken from the YAML file's folder (or
    the folder celSetup --gui started in)."""
    f = (cfg.get("mesh") or {}).get("file")
    if isinstance(f, str) and f and not os.path.isabs(os.path.expanduser(f)):
        base = os.path.dirname(state.yaml_path) if state.yaml_path else os.getcwd()
        cfg["mesh"]["file"] = os.path.join(base, f)
    elif isinstance(f, str):
        cfg["mesh"]["file"] = os.path.expanduser(f)
    return cfg


def geometry(cfg, state=None):
    """Fast outline (no mesh): fill, walls, domain."""
    c = case.resolve(cfg)
    t = case.TEMPLATES[c["template"]]
    b = t.build(c, "gui")
    fill = b["fill"]
    if c["template"] == "userMesh":
        # the fill box, clipped to the loaded mesh
        u = state.user if state else None
        if u is None:
            fill = None
        else:
            (x0, y0), (x1, y1) = u["box"]
            fx0, fy0 = max(fill[0][0], x0), max(fill[0][1], y0)
            fx1, fy1 = min(fill[2][0], x1), min(fill[2][1], y1)
            fill = [[fx0, fy0], [fx1, fy0], [fx1, fy1], [fx0, fy1]] if fx1 > fx0 and fy1 > fy0 \
                else None
    mass = None if b["volume"] is None else b["volume"] * float(c["material"]["rho"])
    return dict(fill=fill, solids=b["solids"], summary=b["summary"],
                wall=b["notes"].get("wall"), mass=mass)


def load_user_mesh(state, cfg):
    """Import the user's mesh into a scratch case and describe it: patches
    with a suggested condition and their edges, cell polygons (2D or wedge),
    bounds. Nothing about the boundaries needs to be known yet."""
    cfg = dict(cfg, boundaries={})
    c = case.resolve(cfg)
    with state.lock:
        d = os.path.join(state.scratch, "usermesh")
        if os.path.isdir(d):
            shutil.rmtree(d)
        case.write_case(c, d, src="celSetup GUI")
        fmt = usermesh.import_mesh(d, c["mesh"])
        desc = usermesh.describe_mesh(d)
        # front polygons need the 2D or wedge sides typed: the suggestions
        from . import gmsh_mesh
        gmsh_mesh.set_patch_types(d, {p["name"]: p["suggest"] for p in desc["patches"]
                                      if p["suggest"] in ("wedge", "empty")})
        m = foam.PolyMesh(d)
        cells, polys, wedge = m.front_polygons()
    if len(cells):
        allp = np.vstack(polys)
        box = [allp.min(0).tolist(), allp.max(0).tolist()]
    else:
        lo, hi = desc["bounds"]
        box = [lo[:2], hi[:2]]
    desc.update(format=fmt, file=c["mesh"]["file"], box=box, wedge=wedge,
                polys=[np.round(q, 6).tolist() for q in polys])
    state.user = desc
    return desc


def meshers():
    """Which meshers can run here, with a note for the ones that cannot."""
    from . import gmsh_mesh
    bm = shutil.which("blockMesh") is not None
    out = {"blockMesh": dict(ok=bm, note="" if bm else "blockMesh not found: source "
                             "OpenFOAM's bashrc and start celSetup --gui again")}
    g = gmsh_mesh.available()
    if g is None:
        out["gmsh"] = dict(ok=False, note="Gmsh not found: pip install gmsh (or apt install gmsh)")
    elif shutil.which("gmshToFoam") is None:
        out["gmsh"] = dict(ok=False, note="gmshToFoam not found: source OpenFOAM's bashrc")
    else:
        out["gmsh"] = dict(ok=True, note=f"Gmsh ({'Python module' if g == 'python' else 'command'})")
    return out


def other_mesher(cfg):
    """The configuration with the other mesher (for the comparison)."""
    c = {k: dict(v) if isinstance(v, dict) else v for k, v in cfg.items()}
    m = c.setdefault("mesh", {})
    m["mesher"] = "blockMesh" if m.get("mesher", "blockMesh") == "gmsh" else "gmsh"
    return c


def mesh_preview(state, cfg, sub="preview"):
    """Write the case in a scratch folder, mesh (blockMesh or Gmsh, or
    import the user's mesh), fill, celCheck."""
    c = case.resolve(cfg)
    with state.lock:
        d = os.path.join(state.scratch, sub)
        if os.path.isdir(d):
            shutil.rmtree(d)
        case.write_case(c, d, src="celSetup GUI")
        ok, msg = case.run_mesh_and_fill(d)
        if not ok:
            return dict(error=msg)
        res, m = check.run_check(d)
        text, worst = check.report(res)
        cells, polys, wedge = m.front_polygons()
        small, _ = check.inplane_sizes(polys)
        a = read_alpha(d, m.nCells)
        return dict(
            polys=[np.round(q, 6).tolist() for q in polys],
            size=np.round(small * 1000, 3).tolist(),
            alpha=np.round(a[cells], 4).tolist(),
            solids=res.get("solids", []),
            checks=res["checks"], verdict=worst, nCells=res["nCells"], dt=res["dt"],
            hours=res["hours"], limit=res["limit"], fill=msg, wedge=wedge,
            mesher=c["mesh"].get("mesher", "blockMesh"), nonOrth=res["nonOrth"]["max"],
            patches=patch_segments(d, wedge) if c["template"] == "userMesh" else None)


def patch_segments(d, wedge):
    """Edges of every non-side patch in the (x, y) section, for drawing."""
    m = foam.PolyMesh(d)
    out = {}
    for name in m.patches:
        if name in m.patches_of_type("wedge") + m.patches_of_type("empty"):
            continue
        segs = []
        for fi in m.patch_faces(name):
            q = m.points[m.faces[fi]]
            r = np.hypot(q[:, 0], q[:, 2]) if wedge else q[:, 0]
            xy = np.unique(np.round(np.column_stack([r, q[:, 1]]), 9), axis=0)
            if len(xy) >= 2:
                segs.append(xy[[0, -1]].tolist())
        out[name] = segs
    return out


def read_alpha(d, n):
    t = open(os.path.join(d, "0", "alpha")).read()
    i = t.index("(", t.index("nonuniform"))
    return np.array(t[i + 1:t.index(")", i)].split(), float)[:n]


def save(state, cfg, yaml_path, case_dir):
    c = case.resolve(cfg)
    t = case.TEMPLATES[c["template"]]
    yaml_path = os.path.abspath(os.path.expanduser(yaml_path))
    case_dir = os.path.abspath(os.path.expanduser(case_dir))
    with open(yaml_path, "w") as f:
        f.write(schema.template_yaml(c["template"], t.SCHEMA, cfg, advanced=True,
                                     boundary_types=getattr(t, "BTYPES", None)) + "\n")
    b = case.write_case(c, case_dir, src=os.path.basename(yaml_path))
    ok, msg = case.run_mesh_and_fill(case_dir)
    state.yaml_path = yaml_path
    return dict(yaml=yaml_path, case=case_dir, summary=b["summary"], mesh=msg, ok=ok)


def start_run(state, case_dir):
    case_dir = os.path.abspath(os.path.expanduser(case_dir))
    if not os.path.exists(os.path.join(case_dir, "Allrun")):
        return dict(error=f"no Allrun in {case_dir}: save the case first")
    r = state.runs.get(case_dir)
    if r and r.poll() is None:
        return dict(error="already running")
    state.runs[case_dir] = subprocess.Popen(
        ["./Allrun"], cwd=case_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        start_new_session=True, env=dict(os.environ, PWD=case_dir))
    return dict(started=case_dir)


def run_status(state, case_dir):
    case_dir = os.path.abspath(os.path.expanduser(case_dir))
    log = os.path.join(case_dir, "log.celFoam")
    out = dict(running=bool(state.runs.get(case_dir) and state.runs[case_dir].poll() is None))
    if os.path.exists(log):
        with open(log, "rb") as f:
            f.seek(max(0, os.path.getsize(log) - 20000))
            tail = f.read().decode(errors="replace")
        times = [l for l in tail.splitlines() if l.startswith("Time = ")]
        if times:
            out["time"] = float(times[-1].split()[2])
        meta = os.path.join(case_dir, "system", "celSetup.json")
        if os.path.exists(meta):
            out["endTime"] = json.load(open(meta)).get("endTime")
        out["tail"] = "\n".join(tail.splitlines()[-12:])
    return out


def make_handler(state):
    class H(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, body, ctype="application/json"):
            data = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            path = self.path.split("?")[0]
            if path in ("/", "/index.html"):
                path = "/index.html"
            if path == "/api/init":
                cfg = initial_config(state.yaml_path)
                sch = {n: schema.to_json(t.SCHEMA) for n, t in case.TEMPLATES.items()}
                presets = {n: t.PRESETS for n, t in case.TEMPLATES.items()}
                base = state.yaml_path or os.path.join(os.getcwd(), "case.yaml")
                return self._send(200, dict(config=cfg, schema=sch, presets=presets,
                                            sections=schema.SECTION_TEXT, yamlPath=base,
                                            caseDir=os.path.join(os.path.dirname(base),
                                                                 cfg.get("case", "celCase")),
                                            openfoam=shutil.which("blockMesh") is not None,
                                            meshers=meshers(),
                                            btypes={k: dict(text=v[0], entries=[
                                                e._asdict() for e in v[1]])
                                                for k, v in usermesh.BTYPES.items()}))
            if path == "/api/status":
                q = dict(p.split("=", 1) for p in self.path.split("?", 1)[1].split("&")) \
                    if "?" in self.path else {}
                from urllib.parse import unquote
                return self._send(200, run_status(state, unquote(q.get("case", ""))))
            f = os.path.normpath(os.path.join(STATIC, path.lstrip("/")))
            if f.startswith(STATIC) and os.path.isfile(f):
                ext = os.path.splitext(f)[1]
                return self._send(200, open(f, "rb").read(), TYPES.get(ext, "text/plain"))
            self._send(404, dict(error="not found"))

        def do_POST(self):
            n = int(self.headers.get("Content-Length", 0))
            try:
                body = json.loads(self.rfile.read(n) or b"{}")
                cfg = abs_mesh_file(state, clean(body.get("config", {})))
                if self.path == "/api/geometry":
                    return self._send(200, geometry(cfg, state))
                if self.path == "/api/loadmesh":
                    return self._send(200, load_user_mesh(state, cfg))
                if self.path == "/api/mesh":
                    return self._send(200, mesh_preview(state, cfg))
                if self.path == "/api/compare":
                    return self._send(200, mesh_preview(state, other_mesher(cfg), "compare"))
                if self.path == "/api/save":
                    return self._send(200, save(state, cfg, body["yamlPath"], body["caseDir"]))
                if self.path == "/api/run":
                    return self._send(200, start_run(state, body["caseDir"]))
                return self._send(404, dict(error="not found"))
            except (ValueError, KeyError, OSError) as e:
                return self._send(200, dict(error=str(e)))
    return H


def serve(yaml_path=None, port=8765, open_browser=True):
    state = State(yaml_path)
    srv = None
    for p in range(port, port + 20):
        try:
            srv = ThreadingHTTPServer(("127.0.0.1", p), make_handler(state))
            break
        except OSError:
            continue
    if srv is None:
        raise OSError(f"no free port from {port}")
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"celSetup GUI: {url}  (Ctrl+C to stop)", flush=True)
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        shutil.rmtree(state.scratch, ignore_errors=True)
    return 0
