// celSetup GUI (WP S.4). Plain JavaScript, no libraries, works offline.
"use strict";

const $ = (id) => document.getElementById(id);
const S = {
  init: null, template: null, values: {}, geom: null, mesh: null, colour: "size",
  view: null, meshReq: 0, geomReq: 0, runCase: null, other: null, otherReq: 0,
  user: null, hl: null,
};
// G.3: boundary condition colours on the drawing
const BCOL = { wall: "#3a3f45", outlet: "#d62728", inlet: "#1f77b4", symmetry: "#2ca02c",
  free: "#9467bd", none: "#ff9800" };
const NAMES = { blockMesh: "blockMesh", gmsh: "Gmsh" };

async function api(path, body) {
  const r = await fetch(path, body === undefined ? {} : {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body),
  });
  return r.json();
}

function toast(msg, ms = 6000) {
  const t = $("toast");
  t.textContent = msg; t.hidden = false;
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => (t.hidden = true), ms);
}

// ---------- form ----------------------------------------------------------
function fmtDefault(v) {
  if (v === null || v === undefined) return "";
  if (v === true) return "yes";
  if (v === false) return "no";
  return String(v);
}

function buildForm() {
  const sch = S.init.schema[S.template];
  const form = $("form");
  form.innerHTML = "";
  for (const [sec, entries] of Object.entries(sch)) {
    const fs = document.createElement("fieldset");
    fs.innerHTML = `<legend>${sec}</legend><div class="sectext">${S.init.sections[S.template + "." + sec] || S.init.sections[sec] || ""}</div>`;
    const adv = document.createElement("details");
    adv.innerHTML = "<summary>Advanced</summary>";
    let nAdv = 0;
    for (const e of entries) {
      const el = fieldEl(sec, e);
      if (e.level === "advanced") { adv.appendChild(el); nAdv++; } else fs.appendChild(el);
    }
    if (nAdv) fs.appendChild(adv);
    if (S.template === "userMesh" && sec === "mesh") {
      const b = document.createElement("div");
      b.className = "field";
      b.innerHTML = `<button id="loadMesh" class="primary">Load mesh</button> <span id="meshInfo" class="desc"></span>`;
      b.querySelector("button").onclick = loadMesh;
      fs.appendChild(b);
    }
    form.appendChild(fs);
    if (S.template === "userMesh" && sec === "mesh") {
      const bs = document.createElement("fieldset");
      bs.id = "bndSet";
      form.appendChild(bs);
      buildBoundaries();
    }
  }
}

// the mesher choice belongs to templates that make their own mesh
function meshToolbar() {
  const own = S.template !== "userMesh";
  $("mesher").closest("label").hidden = !own;
  $("compare").closest("label").hidden = !own;
  if (!own) $("compare").checked = false;
}

// ---------- user mesh (G.3) ------------------------------------------------
async function loadMesh() {
  $("meshInfo").textContent = "loading…";
  const r = await api("/api/loadmesh", { config: config() });
  if (r.error) { $("meshInfo").textContent = ""; showError(r.error); return; }
  S.user = r; S.mesh = null; S.other = null;
  S.values.boundaries = S.values.boundaries || {};
  for (const p of r.patches) {
    if (!S.values.boundaries[p.name] && p.suggest) S.values.boundaries[p.name] = { type: p.suggest };
  }
  $("meshInfo").textContent = `${r.format}: ${r.nCells.toLocaleString()} cells, ${r.patches.length} patches, ` +
    (r.kind === "3D" ? "3D (not supported yet: give the front and back patches type empty)" : r.kind);
  buildBoundaries();
  S.moved = false; fit();
  changed();
}

function buildBoundaries() {
  const bs = $("bndSet");
  if (!bs) return;
  bs.innerHTML = `<legend>boundaries</legend><div class="sectext">A condition for every patch of the mesh. ` +
    `Point at a patch to see it on the drawing.</div>`;
  if (!S.user) { bs.innerHTML += `<div class="desc">Give the mesh file above and press Load mesh.</div>`; return; }
  const bt = S.init.btypes;
  for (const p of S.user.patches) {
    const cur = (S.values.boundaries || {})[p.name] || {};
    const div = document.createElement("div");
    div.className = "field bnd" + (cur.type ? "" : " bad");
    div.dataset.patch = p.name;
    const sel = `<select id="b_${p.name}"><option value="">choose…</option>` + Object.keys(bt).map((t) =>
      `<option ${cur.type === t ? "selected" : ""}>${t}</option>`).join("") + "</select>";
    div.innerHTML = `<div class="row"><label for="b_${p.name}" title="${p.foamType}, ${p.nFaces} faces">` +
      `<span class="swatch" style="background:${cur.type === "wedge" || cur.type === "empty" ? "transparent" : (BCOL[cur.type] || BCOL.none)}"></span>${p.name}</label>${sel}<span class="unit">${p.nFaces}</span></div>` +
      `<div class="desc">${cur.type ? bt[cur.type].text : "no condition yet"}</div>`;
    for (const e of (cur.type ? bt[cur.type].entries : [])) {
      const req = e.default === "required";
      const ph = req ? "required" : (e.default === null ? "off" : fmtDefault(e.default));
      const row = document.createElement("div");
      row.className = "row sub";
      row.innerHTML = `<label title="${e.text}">${e.key}</label><input value="${cur[e.key] ?? ""}" placeholder="${ph}" ` +
        `title="${e.text}" spellcheck="false"><span class="unit">${e.unit || ""}</span>`;
      row.querySelector("input").addEventListener("input", (ev) => {
        const v = ev.target.value.trim();
        if (v === "") delete S.values.boundaries[p.name][e.key]; else S.values.boundaries[p.name][e.key] = v;
        changed();
      });
      div.appendChild(row);
    }
    div.querySelector("select").addEventListener("change", (ev) => {
      S.values.boundaries = S.values.boundaries || {};
      const t = ev.target.value;
      if (!t) delete S.values.boundaries[p.name];
      else S.values.boundaries[p.name] = { type: t };
      buildBoundaries(); changed();
    });
    div.addEventListener("mouseenter", () => { S.hl = p.name; draw(); });
    div.addEventListener("mouseleave", () => { S.hl = null; draw(); });
    bs.appendChild(div);
  }
}

function missingBoundaries() {
  if (!S.user) return null;
  return S.user.patches.filter((p) => !((S.values.boundaries || {})[p.name] || {}).type).map((p) => p.name);
}

function fieldEl(sec, e) {
  const div = document.createElement("div");
  const given = (S.values[sec] || {})[e.key];
  const req = e.default === "required";
  div.className = "field" + (req ? " required" : "") + (given !== undefined ? " changed" : "");
  const id = `f_${sec}_${e.key}`;
  let input;
  if (e.choices) {
    input = `<select id="${id}">` + e.choices.map((c) =>
      `<option ${fmtDefault(given ?? e.default) === c ? "selected" : ""}>${c}</option>`).join("") + "</select>";
  } else {
    const ph = req ? "required" : (e.default === null ? "off" : fmtDefault(e.default));
    const val = given === undefined || given === null ? "" : fmtDefault(given);
    input = `<input id="${id}" value="${val.replace(/"/g, "&quot;")}" placeholder="${ph}" spellcheck="false">`;
  }
  div.innerHTML = `<div class="row"><label for="${id}" title="${e.key}">${e.key}</label>${input}` +
    `<span class="unit">${e.unit || ""}</span></div><div class="desc">${e.text}` +
    (req ? "" : ` <i>Default: ${e.default === null ? "off" : fmtDefault(e.default)}.</i>`) + "</div>";
  const inp = div.querySelector("input,select");
  inp.addEventListener("input", () => {
    S.values[sec] = S.values[sec] || {};
    const v = inp.value.trim();
    if (v === "" || (e.choices && v === fmtDefault(e.default) && given === undefined)) delete S.values[sec][e.key];
    else S.values[sec][e.key] = v;
    div.classList.toggle("changed", S.values[sec][e.key] !== undefined);
    changed();
  });
  return div;
}

function config() {
  const c = { case: S.values.case || "celCase", template: S.template };
  for (const sec of Object.keys(S.init.schema[S.template])) c[sec] = { ...(S.values[sec] || {}) };
  if (S.template === "userMesh") c.boundaries = JSON.parse(JSON.stringify(S.values.boundaries || {}));
  return c;
}

let tGeom = null, tMesh = null;
function changed() {
  clearTimeout(tGeom); clearTimeout(tMesh);
  $("mesher").value = mesherNow();
  S.other = null; S.otherReq++;
  tGeom = setTimeout(updateGeometry, 150);
  tMesh = setTimeout(updateMesh, 1200);
}

function mesherNow() { return (S.values.mesh || {}).mesher || "blockMesh"; }

async function updateGeometry() {
  const id = ++S.geomReq;
  const g = await api("/api/geometry", { config: config() });
  if (id !== S.geomReq) return;
  if (g.error) { showError(g.error); return; }
  S.geom = g;
  $("summary").innerHTML = `<b>${g.summary}</b>` + (g.mass != null ? `; material ${Math.round(g.mass).toLocaleString()} kg` : "");
  if (!S.view) fit();
  draw();
}

async function updateMesh() {
  if (!S.init.openfoam) {
    $("checks").innerHTML = `<div class="error">OpenFOAM is not loaded (blockMesh not found): ` +
      `the mesh preview and checks need it. Source OpenFOAM's bashrc and start celSetup --gui again.</div>`;
    return;
  }
  if (S.template === "userMesh") {
    const miss = missingBoundaries();
    if (miss === null || miss.length) {
      S.meshReq++;
      $("checks").innerHTML = `<div class="error">${miss === null ? "Load the mesh first (Mesh: file, then Load mesh)."
        : "Give a condition for every patch: " + miss.join(", ")}. The checks run once every patch has one.</div>`;
      return;
    }
  }
  const id = ++S.meshReq;
  $("busy").hidden = false;
  const m = await api("/api/mesh", { config: config() });
  if (id !== S.meshReq) return;
  $("busy").hidden = true;
  if (m.error) { showError(m.error); return; }
  S.mesh = m;
  showChecks();
  draw();
  if ($("compare").checked) updateCompare();
}

// G.2: the same case meshed with the other mesher, drawn on top
async function updateCompare() {
  const id = ++S.otherReq;
  const other = mesherNow() === "gmsh" ? "blockMesh" : "gmsh";
  const info = (S.init.meshers || {})[other];
  if (info && !info.ok) { toast(`Cannot compare: ${info.note}`); $("compare").checked = false; return; }
  $("busy").textContent = `building the ${NAMES[other]} mesh…`; $("busy").hidden = false;
  const m = await api("/api/compare", { config: config() });
  $("busy").hidden = true; $("busy").textContent = "building mesh…";
  if (id !== S.otherReq) return;
  if (m.error) { toast(`${NAMES[other]}: ${m.error}`, 9000); return; }
  S.other = m;
  showChecks();
  draw();
}

function showError(msg) {
  $("busy").hidden = true;
  $("checks").innerHTML = `<div class="error">${msg}</div>`;
  for (const el of document.querySelectorAll(".field")) el.classList.remove("bad");
  const m = /^(\w+): '(\w+)'/.exec(msg) || /^(\w+): unknown entr\w+ (\w+)/.exec(msg);
  if (m) { const el = $(`f_${m[1]}_${m[2]}`); if (el) el.closest(".field").classList.add("bad"); }
}

function showChecks() {
  const m = S.mesh;
  for (const el of document.querySelectorAll(".field")) el.classList.remove("bad");
  const h4 = m.hours["4"] ?? m.hours[4];
  $("summary").innerHTML = (S.geom ? `<b>${S.geom.summary}</b>; ` : "") +
    `${m.nCells.toLocaleString()} cells, time step ${m.dt.toExponential(3)} s, about ` +
    `${h4 >= 1 ? h4.toFixed(1) + " h" : Math.round(h4 * 60) + " min"} on 4 cores. ${m.fill}. ` +
    `Verdict: <span class="badge ${m.verdict}" style="display:inline-block;width:46px">${m.verdict}</span>`;
  $("checks").innerHTML = compareTable() + m.checks.map((c) =>
    `<div class="check"><span class="badge ${c.status}">${c.status}</span><span><b>${c.name}</b>: ${c.text}</span></div>`).join("");
}

function compareTable() {
  const a = S.mesh, b = S.other;
  if (!a || !b || !$("compare").checked) return "";
  const h = (m) => { const v = m.hours["4"] ?? m.hours[4]; return v >= 1 ? v.toFixed(1) + " h" : Math.round(v * 60) + " min"; };
  const rows = [
    ["cells", (m) => m.nCells.toLocaleString()],
    ["time step (s)", (m) => m.dt.toExponential(4)],
    ["non-orthogonality max (deg)", (m) => m.nonOrth.toFixed(1)],
    ["run time on 4 cores", h],
    ["verdict", (m) => `<span class="badge ${m.verdict}">${m.verdict}</span>`],
  ];
  return `<table class="cmp"><tr><th></th><th>${NAMES[a.mesher]} (colours)</th>` +
    `<th>${NAMES[b.mesher]} (<span class="otherline">lines</span>)</th></tr>` +
    rows.map(([k, f]) => `<tr><td>${k}</td><td>${f(a)}</td><td>${f(b)}</td></tr>`).join("") + "</table>";
}

// ---------- drawing -------------------------------------------------------
const SPECTRAL = [[158, 1, 66], [213, 62, 79], [244, 109, 67], [253, 174, 97], [254, 224, 139],
  [230, 245, 152], [171, 221, 164], [102, 194, 165], [50, 136, 189], [94, 79, 162]];
const SAND = [[250, 246, 236], [240, 214, 150], [214, 160, 70], [150, 95, 30]];
function ramp(stops, t) {
  t = Math.min(1, Math.max(0, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(t)), f = t - i;
  const a = stops[i], b = stops[i + 1];
  return `rgb(${a.map((x, k) => Math.round(x + f * (b[k] - x))).join(",")})`;
}

function bounds() {
  if (S.user) { const [[x0, y0], [x1, y1]] = S.user.box; return [x0, y0, x1, y1]; }
  const pts = [];
  if (S.geom) { for (const s of S.geom.solids) pts.push(...s); pts.push(...S.geom.fill); }
  if (!pts.length) return [0, 0, 1, 1];
  let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  for (const [x, y] of pts) { x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
  return [x0, y0, x1, y1];
}

function fit() {
  const c = $("canvas"), [x0, y0, x1, y1] = bounds();
  const w = c.clientWidth, h = c.clientHeight, pad = 30;
  const s = Math.min((w - 2 * pad) / Math.max(x1 - x0, 1e-9), (h - 2 * pad) / Math.max(y1 - y0, 1e-9));
  S.view = { s, ox: pad - x0 * s + (w - 2 * pad - (x1 - x0) * s) / 2, oy: h - pad + y0 * s - (h - 2 * pad - (y1 - y0) * s) / 2 };
  draw();
}

function draw() {
  const c = $("canvas"), dpr = window.devicePixelRatio || 1;
  const w = c.clientWidth, h = c.clientHeight;
  if (c.width !== w * dpr || c.height !== h * dpr) { c.width = w * dpr; c.height = h * dpr; }
  const g = c.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);
  if (!S.view) return;
  const { s, ox, oy } = S.view;
  const X = (x) => ox + x * s, Y = (y) => oy - y * s;
  const dark = matchMedia("(prefers-color-scheme: dark)").matches;
  const mesh = S.mesh;
  if (mesh) {
    const vals = S.colour === "size" ? mesh.size : mesh.alpha;
    let lo = 0, hi = 1;
    if (S.colour === "size") {
      const sorted = [...mesh.size].sort((a, b) => a - b);
      lo = sorted[0]; hi = sorted[Math.floor(0.98 * (sorted.length - 1))];
    }
    const edges = $("edges").checked && s * (mesh.size[0] / 1000) > 3;
    g.lineWidth = 0.5;
    g.strokeStyle = dark ? "rgba(0,0,0,.45)" : "rgba(60,70,80,.45)";
    for (let i = 0; i < mesh.polys.length; i++) {
      const q = mesh.polys[i];
      g.beginPath();
      g.moveTo(X(q[0][0]), Y(q[0][1]));
      for (let k = 1; k < q.length; k++) g.lineTo(X(q[k][0]), Y(q[k][1]));
      g.closePath();
      const t = (vals[i] - lo) / Math.max(hi - lo, 1e-12);
      g.fillStyle = S.colour === "size" ? ramp(SPECTRAL, t) : ramp(SAND, t);
      g.fill();
      if (edges) g.stroke();
    }
    legend(lo, hi);
    if (S.other && $("compare").checked) {
      // only when the cells are big enough on screen to tell the lines apart
      const visible = s * (S.other.size[0] / 1000) > 4;
      if (visible) {
        g.strokeStyle = dark ? "#ff7bd5" : "#c2187a"; g.lineWidth = 0.9;
        for (const q of S.other.polys) { poly(g, q, X, Y); g.stroke(); }
      }
      $("legend").innerHTML += `<div><span class="otherline">━</span> ${NAMES[S.other.mesher]} cell edges` +
        (visible ? "" : " (zoom in to see them)") + "</div>";
    }
    const L = mesh.limit;
    const lx = L.r ?? L.centre[0], ly = L.y ?? L.centre[1];
    g.strokeStyle = dark ? "#fff" : "#000"; g.lineWidth = 2;
    g.beginPath(); g.arc(X(lx), Y(ly), 9, 0, 2 * Math.PI); g.stroke();
  } else if (S.user) {
    g.lineWidth = 0.5; g.strokeStyle = dark ? "rgba(200,210,220,.5)" : "rgba(60,70,80,.45)";
    g.fillStyle = dark ? "#262c33" : "#eef1f4";
    for (const q of S.user.polys) { poly(g, q, X, Y); g.fill(); if (s * 0.004 > 3) g.stroke(); }
    if (S.geom && S.geom.fill) {
      g.fillStyle = dark ? "rgba(214,160,70,.45)" : "rgba(214,160,70,.55)";
      poly(g, S.geom.fill, X, Y); g.fill();
    }
    $("legend").innerHTML = "<div>your mesh (checks run when every patch has a condition)</div>";
  } else if (S.geom && S.geom.fill) {
    g.fillStyle = dark ? "rgba(214,160,70,.45)" : "rgba(214,160,70,.55)";
    poly(g, S.geom.fill, X, Y); g.fill();
    $("legend").innerHTML = "<div>geometry (mesh building…)</div>";
  }
  if (S.user) drawPatches(g, X, Y, dark);
  if (S.geom) {
    g.strokeStyle = dark ? "#ef6f61" : "#c0392b"; g.lineWidth = 1.5;
    for (const sp of S.geom.solids) { poly(g, sp, X, Y); g.stroke(); }
  }
}

function drawPatches(g, X, Y, dark) {
  const segs = {};
  if (S.mesh && S.mesh.patches) Object.assign(segs, S.mesh.patches);
  else for (const p of S.user.patches) segs[p.name] = p.segments;
  const used = new Set();
  for (const [name, list] of Object.entries(segs)) {
    const t = ((S.values.boundaries || {})[name] || {}).type;
    if (t === "wedge" || t === "empty") continue;
    const hl = S.hl === name;
    let col = BCOL[t] || BCOL.none;
    if (t === "wall" && dark) col = "#d0d6dd";
    used.add(t || "none");
    g.strokeStyle = hl ? (dark ? "#ffff66" : "#ffbf00") : col;
    g.lineWidth = hl ? 6 : 3;
    g.setLineDash(t === "symmetry" ? [6, 4] : []);
    g.beginPath();
    for (const sg of list) { g.moveTo(X(sg[0][0]), Y(sg[0][1])); for (const q of sg.slice(1)) g.lineTo(X(q[0]), Y(q[1])); }
    g.stroke();
  }
  g.setLineDash([]);
  $("legend").innerHTML += [...used].map((t) =>
    `<div><span class="swatch" style="background:${BCOL[t]}"></span>${t === "none" ? "no condition yet" : t}</div>`).join("");
}

function poly(g, pts, X, Y) {
  g.beginPath();
  g.moveTo(X(pts[0][0]), Y(pts[0][1]));
  for (const [x, y] of pts.slice(1)) g.lineTo(X(x), Y(y));
  g.closePath();
}

function legend(lo, hi) {
  const stops = S.colour === "size" ? SPECTRAL : SAND;
  const grad = stops.map((c, i) => `rgb(${c.join(",")}) ${(100 * i / (stops.length - 1)).toFixed(0)}%`).join(",");
  const title = S.colour === "size" ? "cell size, smaller side (mm)" : "material fraction at the start";
  const f = (v) => S.colour === "size" ? v.toFixed(1) : v.toFixed(1);
  $("legend").innerHTML = `<div>${title}</div><div class="bar" style="background:linear-gradient(90deg,${grad})"></div>` +
    `<div class="ends"><span>${f(lo)}</span><span>${f(hi)}${S.colour === "size" ? "+" : ""}</span></div>` +
    `<div style="margin-top:4px">○ cell that sets the time step</div>`;
}

function setupCanvas() {
  const c = $("canvas");
  let drag = null;
  c.addEventListener("wheel", (ev) => {
    ev.preventDefault();
    if (!S.view) return;
    const r = c.getBoundingClientRect(), mx = ev.clientX - r.left, my = ev.clientY - r.top;
    const k = Math.exp(-ev.deltaY * 0.0015);
    S.view.ox = mx - (mx - S.view.ox) * k; S.view.oy = my - (my - S.view.oy) * k; S.view.s *= k;
    S.moved = true;
    draw();
  }, { passive: false });
  c.addEventListener("mousedown", (ev) => { drag = [ev.clientX, ev.clientY]; c.style.cursor = "grabbing"; });
  window.addEventListener("mouseup", () => { drag = null; c.style.cursor = "grab"; });
  window.addEventListener("mousemove", (ev) => {
    if (drag && S.view) {
      S.view.ox += ev.clientX - drag[0]; S.view.oy += ev.clientY - drag[1];
      S.moved = true;
      drag = [ev.clientX, ev.clientY]; draw();
    }
    hover(ev);
  });
  // refit while the user has not zoomed or panned; the view keeps its size
  // changes (checks appearing below it) from cutting the drawing
  new ResizeObserver(() => (S.moved ? draw() : fit())).observe(c);
}

function hover(ev) {
  const c = $("canvas"), r = c.getBoundingClientRect(), hv = $("hover");
  if (!S.view || !S.mesh || ev.clientX < r.left || ev.clientX > r.right || ev.clientY < r.top || ev.clientY > r.bottom) {
    hv.hidden = true; return;
  }
  const x = (ev.clientX - r.left - S.view.ox) / S.view.s, y = (S.view.oy - (ev.clientY - r.top)) / S.view.s;
  const m = S.mesh;
  for (let i = 0; i < m.polys.length; i++) {
    if (inside(x, y, m.polys[i])) {
      hv.hidden = false;
      hv.textContent = `${m.wedge ? "r" : "x"} ${(x * 1000).toFixed(1)} mm, y ${(y * 1000).toFixed(1)} mm: ` +
        `cell ${m.size[i].toFixed(2)} mm, material ${m.alpha[i].toFixed(2)}`;
      return;
    }
  }
  hv.hidden = true;
}

function inside(x, y, q) {
  let c = false;
  for (let i = 0, j = q.length - 1; i < q.length; j = i++) {
    if ((q[i][1] > y) !== (q[j][1] > y) && x < q[i][0] + (y - q[i][1]) * (q[j][0] - q[i][0]) / (q[j][1] - q[i][1])) c = !c;
  }
  return c;
}

// ---------- save and run ---------------------------------------------------
function setupButtons() {
  $("saveBtn").onclick = () => {
    $("yamlPath").value = $("yamlPath").value || S.init.yamlPath;
    $("caseDir").value = $("caseDir").value || S.init.caseDir;
    $("saveDlg").showModal();
  };
  $("saveGo").onclick = async (ev) => {
    ev.preventDefault();
    $("saveDlg").close();
    toast("Saving and meshing…");
    const r = await api("/api/save", { config: config(), yamlPath: $("yamlPath").value, caseDir: $("caseDir").value });
    toast(r.error ? "Not saved: " + r.error : `Saved ${r.yaml}\nCase ${r.case}\n${r.mesh}`, 9000);
    if (!r.error) S.runCase = r.case;
  };
  $("runBtn").onclick = async () => {
    const dir = S.runCase || $("caseDir").value || S.init.caseDir;
    if (!confirm(`Start ./Allrun in\n${dir}\n?\nIt runs in the background; closing this page does not stop it.`)) return;
    const r = await api("/api/run", { caseDir: dir });
    if (r.error) { toast(r.error); return; }
    S.runCase = dir;
    poll();
  };
  for (const el of document.querySelectorAll("input[name=colour]")) {
    el.onchange = () => { S.colour = el.value; draw(); };
  }
  $("edges").onchange = draw;
  $("fitBtn").onclick = () => { S.moved = false; fit(); };
  $("limBtn").onclick = () => {
    if (!S.mesh) return;
    const L = S.mesh.limit, c = $("canvas"), [x0, y0, x1, y1] = bounds();
    const lx = L.r ?? L.centre[0], ly = L.y ?? L.centre[1];
    const w = 0.12 * Math.max(x1 - x0, y1 - y0), s = Math.min(c.clientWidth, c.clientHeight) / w;
    // the cell a little left of and below the middle, so its surroundings show
    S.view = { s, ox: c.clientWidth * 0.35 - lx * s, oy: c.clientHeight * 0.7 + ly * s };
    S.moved = true; draw();
  };
  $("mesher").onchange = () => {
    const f = $("f_mesh_mesher");
    if (f) { f.value = $("mesher").value; f.dispatchEvent(new Event("input")); }
  };
  $("compare").onchange = () => {
    if ($("compare").checked) updateCompare(); else { S.other = null; S.otherReq++; showChecks(); draw(); }
  };
  $("template").onchange = () => {
    // keep the material and run entries; the template's own sections start afresh
    S.template = $("template").value;
    S.values = { case: S.values.case, material: S.values.material, run: S.values.run };
    S.user = null; S.mesh = null; S.other = null; S.geom = null; S.view = null;
    buildForm(); meshToolbar(); changed(); setTimeout(fit, 300);
  };
}

async function poll() {
  const r = await api("/api/status?case=" + encodeURIComponent(S.runCase));
  const frac = r.time && r.endTime ? r.time / r.endTime : 0;
  toast(`Run ${r.running ? "running" : "stopped"}: t = ${r.time ?? 0} of ${r.endTime ?? "?"} s ` +
    `(${(100 * frac).toFixed(1)} %)\n${S.runCase}`, 5000);
  if (r.running) setTimeout(poll, 4000);
}

// ---------- start ------------------------------------------------------------
async function main() {
  S.init = await api("/api/init");
  const cfg = S.init.config;
  S.template = cfg.template || Object.keys(S.init.schema)[0];
  $("template").innerHTML = Object.keys(S.init.schema).map((t) => `<option>${t}</option>`).join("");
  $("template").value = S.template;
  S.values = { case: cfg.case };
  for (const sec of Object.keys(S.init.schema[S.template])) {
    S.values[sec] = {};
    for (const [k, v] of Object.entries(cfg[sec] || {})) S.values[sec][k] = fmtDefault(v);
  }
  if (cfg.boundaries) {
    S.values.boundaries = {};
    for (const [p, d] of Object.entries(cfg.boundaries)) {
      S.values.boundaries[p] = {};
      for (const [k, v] of Object.entries(typeof d === "object" ? d : { type: d })) S.values.boundaries[p][k] = fmtDefault(v);
    }
  }
  buildForm();
  for (const o of $("mesher").options) {
    const info = (S.init.meshers || {})[o.value];
    if (info && !info.ok) { o.disabled = true; o.textContent += " (not found)"; o.title = info.note; }
  }
  $("mesher").value = mesherNow();
  meshToolbar();
  setupCanvas();
  setupButtons();
  if (S.template === "userMesh" && (S.values.mesh || {}).file) { await loadMesh(); return; }
  await updateGeometry();
  fit();
  updateMesh();
}
main();
