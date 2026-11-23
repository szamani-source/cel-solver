"""WP S.4 checks: the setup GUI's server answers every call the page makes.

Starts `celSetup --gui` on a free port without a browser, loads the page and
its scripts, then: init (schema, presets, the example), geometry (fast
outline), mesh (blockMesh, fill and celCheck in a scratch folder: 7,584
cells, celCheck's verdicts), an error for a missing required entry, and save
(YAML with every entry, which celSetup reads back to the same case).

    python3 tests/celsetup/gui_checks.py WORKDIR
"""
import json
import os
import socket
import subprocess
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from celsetup import case, miniyaml  # noqa: E402

fails = []


def check(ok, msg):
    print(("PASS  " if ok else "FAIL  ") + msg, flush=True)
    if not ok:
        fails.append(msg)


def main(work):
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    srv = subprocess.Popen([sys.executable, os.path.join(ROOT, "tools", "celSetup"), "--gui",
                            "--no-browser", "--port", str(port)], cwd=work)
    base = f"http://127.0.0.1:{port}"

    def get(path):
        return urllib.request.urlopen(base + path, timeout=60).read()

    def post(path, body):
        req = urllib.request.Request(base + path, json.dumps(body).encode(),
                                     {"Content-Type": "application/json"})
        return json.loads(urllib.request.urlopen(req, timeout=120).read())
    try:
        for _ in range(50):
            try:
                get("/api/init")
                break
            except OSError:
                time.sleep(0.2)
        page = get("/").decode()
        js = get("/app.js").decode()
        css = get("/style.css").decode()
        check("app.js" in page and "function draw" in js and "prefers-color-scheme" in css
              and "http" not in page.replace("http-equiv", ""),
              "page, script and style served; no external resources")
        init = json.loads(get("/api/init"))
        sch = init["schema"]["conicalHopper"]
        n = sum(len(v) for v in sch.values())
        check(n >= 57 and all(e["text"] for v in sch.values() for e in v),
              f"init: schema with {n} entries, each with a description")
        cfg = init["config"]
        check(cfg["template"] == "conicalHopper", "init: the Couto example is loaded")
        # form values arrive as strings
        form = {k: ({kk: str(vv) for kk, vv in v.items()} if isinstance(v, dict) else v)
                for k, v in cfg.items()}
        g = post("/api/geometry", {"config": form})
        check("error" not in g and abs(g["mass"] - 1503.4) < 1, f"geometry: {g.get('summary')}")
        t0 = time.time()
        m = post("/api/mesh", {"config": form})
        dt = time.time() - t0
        ok = "error" not in m and m.get("nCells") == 7584 and len(m["polys"]) == 7584
        check(ok, f"mesh preview: {m.get('nCells')} cells, dt {m.get('dt')}, "
              f"verdict {m.get('verdict')}, {dt:.1f} s")
        if ok:
            st = {c["name"]: c["status"] for c in m["checks"]}
            check(st.get("outlet 'outlet'") == "PASS" and abs(m["dt"] / 5.85003e-6 - 1) < 1e-4,
                  "mesh preview carries celCheck's results")
        bad = json.loads(json.dumps(form))
        bad["material"]["phi"] = ""
        e = post("/api/geometry", {"config": bad})
        check("'phi' is required" in e.get("error", ""), f"missing entry reported: {e.get('error')}")
        y = os.path.join(work, "saved.yaml")
        d = os.path.join(work, "savedCase")
        r = post("/api/save", {"config": form, "yamlPath": y, "caseDir": d})
        check("error" not in r and os.path.exists(os.path.join(d, "Allrun"))
              and os.path.exists(os.path.join(d, "0", "alpha")), f"save: {r.get('mesh', r)}")
        back = case.resolve(miniyaml.load(y))
        d2 = os.path.join(work, "fromYaml")
        case.write_case(back, d2)
        def body(p):        # without the banner, which names the YAML file
            t = open(p).read()
            return t[t.index("FoamFile"):]
        same = all(body(os.path.join(d, f)) == body(os.path.join(d2, f))
                   for f in ("system/blockMeshDict", "constant/materialProperties"))
        check(same, "saved YAML reads back to the same mesh and material")
    finally:
        srv.terminate()
        srv.wait()
    print("FAILED: " + "; ".join(fails) if fails else "ALL PASS")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1])))
