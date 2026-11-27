"""Initially fluidised dam break in air (WP F.8, decisions 0125, 0133).

Shared by tests/gas/collapse_checks.py (done-when test,
tests/manual/test_fluidised_collapse.sh).

Targets: Roche, Montserrat, Nino and Tamburrino (2008), J. Geophys. Res.
113, B12203, doi 10.1029/2008JB005664 (HAL hal-00381402; the F.0 table cites
it as 2010). Only the abstract could be read: glass beads 60 to 90
micrometres, columns initially fluidised and expanded 2.5 to 4.5 %, height
over length 0.5 to 3; the front moves at a constant U = sqrt(2 g h0) from
t/sqrt(h0/g) = 1.5 to about 4, the flow stops at about 6.5 with run-out
x/h0 = 5.5 to 6. We read x as measured from the gate.

Setup (F.8c, decision 0133) from the open-access restatement of the same
experiments, Aravena, Chupin, Dubois and Roche (2021), Bull. Volcanol. 83,
doi 10.1007/s00445-021-01495-1: reservoir 0.2 m long, channel 3 m, smooth
horizontal floor (their model: basal friction 15 deg), beads 75 um, rho_s
2500, bulk density 1450 (gas fraction 0.42), the gate opens in under 0.1 s,
the initial pore pressure 0.9 of the weight (measured, Montserrat et al.
2012), and the air supply at the reservoir base kept on after release. Their
model needs an effective diffusivity about 9 times below the theoretical one
to reach the measured run-out; with the theoretical one (as here) it reaches
less than 60 %. They report u_max/sqrt(g h0) 1.0 to 1.15 and a run-out of
about 3.8 h0 for a = 2 on the horizontal (reported, not checked).

Case: a column X0 = 0.2 m long and H0 = 0.4 m high (a = 2) against a smooth
left wall, Mohr-Coulomb with mu(I) (phi 24, mu2 tan 33, I0 0.3, c 0, psi 0, E 1 MPa;
flowImpedance, fullBlend; F.8d) with a compaction cap for the expanded
packing (CAP_LAMBDA below) and overfillStress, floor friction 19 deg, open
top. The column starts at rest with p = 0.9 rho_b g (H0 - y) and the rest of
the weight on the grains (sigma'_yy = -0.1 rho_b g (H0 - y), sigma'_xx = K0
times that). The reservoir floor (x < X0, patch floorRes, made by topoSet and
createPatch) blows air at the superficial velocity that carries 0.9 of the
weight through the packing (Darcy, Q_SUPPLY); the channel floor is sealed.
The gate opens at once (no gate). 10 mm cells, NY = 45 rows, NX = 300 columns
(3 m, the channel). Every gas run uses ergunWenYu. Route B (default, F.8e):
model poreGas with the pore pressure carried by the grains and diffusing at
KAPPA_FIT; kappa = 0 gives the compressible gas of F.8c (model gasFlow,
surfaceCarry yes; decision 0133).

Dry (gas off): same column, no gasProperties, no floorRes, geostatic K0 stress.
"""

import math
import os
import re
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402
import fluid_cases as fl  # noqa: E402

G_ACC = 9.81
D = 7.5e-5
RHOP = 2500.0
E0 = 0.42
RHO_B = (1.0 - E0) * RHOP
X0, H0 = 0.2, 0.4
CELL = 0.01
NX, NY = 300, 45
N = NX * NY
E_MOD = 1.0e6
# Grains (F.8d, decision 0133): mu(I) of glass beads with the static angle
# measured by Gueugneau et al. (2017) (24 deg, slope of a pile at rest),
# mu2 tan 33 deg and I0 0.3 (Jop et al. 2006 as quoted by Chupin and Dubois
# 2023), floor 19 deg (Gueugneau: beads glued on a plate sliding on the
# channel base). With flowImpedance and fullBlend the dry column runs out
# 2.35 h0 against Roche's pressure-balanced 1.05 m (2.6 h0, Breard et al.
# 2019); plain Mohr-Coulomb with these angles stops at 0.9 h0
PHI = 24.0
MU2 = math.tan(math.radians(33.0))
I0 = 0.3
FLOOR = 22.0
# WP F.8k (decision 0135, Sia 2026-11-27): the dam break runs with
# yieldShearJump no and slowFlowJump yes (SLOW_FLOW); with the face damping
# of slow flow gone the dry column ran 3.27 h0 on 5 mm cells with
# Gueugneau's 19 deg floor (+26 %), 2.73 h0 (+5 %) with 22 deg, which Sia
# chose. The v1 dry reference (make(cap=False), check B) keeps 19 deg and
# the v1 options
FLOOR_V1 = 19.0
SLOW_FLOW = "yieldShearJump no;\nslowFlowJump yes;"
DRY_RUNOUT = 1.05 / 0.4                  # x/h0, pressure-balanced (Breard et al. 2019)
P_FRAC = 1.0                          # initial pore pressure / weight (Roche 2008: fully supported)
T_SCALE = math.sqrt(H0 / G_ACC)
WRITE = 0.02
T_END = WRITE * math.ceil(9.0 * T_SCALE / WRITE - 1e-9)   # 1.82 s; Roche: stops at 6.5
MU = gr.air_viscosity(pc.T)
RG = gr.air_density(pc.P0, pc.T)
# Roche et al. (2008) abstract
U_FRONT = math.sqrt(2.0 * G_ACC * H0)
RUNOUT = 5.75                         # x/h0 5.5 to 6
PHASE = (1.5, 4.0)                    # constant-velocity phase, t/sqrt(h0/g)
# Aravena et al. (2021), reported only
ARAVENA_U = (1.0, 1.15)               # u_max / sqrt(g h0)
ARAVENA_RUNOUT = 3.8                  # x/h0, a = 2, horizontal
# Route B (F.8e, decision 0133): pore pressure carried with the grains and
# diffusing (gasProperties diffusivity, model poreGas), as Gueugneau et al.
# (2017) and Aravena et al. (2021). KAPPA_FIT is fitted to Roche's run-out
# (Aravena fit 0.015 m^2/s); kappa_theory() is k p0 / (mu e).
KAPPA_FIT = 0.22                      # F.8h refit on 5 mm cells (theory k p0/(mu e) 0.11)
# F.8g: route B as published, the pore pressure only lowers the friction
FRICTION_ONLY = True
# Expanded (fluidised) packing: soft until compacted (decision 0133): the
# compaction cap (WP K.1, decision 0113) normally consolidated from 1 Pa, with
# lambda* = 0.0044 so that loading to the column's base stress (about 3 kPa)
# compacts it by the 3.5 % it was expanded (0.035 / ln 3000)
CAP_LAMBDA, CAP_P0 = 0.0044, 1.0
MC = (f"type mohrCoulomb;\nE {E_MOD!r};\nnu 0.3;\nphi {PHI!r};\nc 0;\npsi 0;\n"
      f"mu2 {MU2!r};\nI0 {I0!r};\nd {D!r};\nrhoS {RHOP!r};\nflowImpedance yes;\nfullBlend yes;\n"
      f"capLambda {CAP_LAMBDA!r};\ncapP0 {CAP_P0!r};\noverfillStress yes;")


def diffusion_time():
    """H0^2/c with the rigid-bed diffusivity c = k p0/(mu e) [s]."""
    k = gr.permeability(D, E0)
    return H0 * H0 / (k * pc.P0 / (MU * E0))


def kappa_theory():
    """Pore pressure diffusivity of the rigid packing, k p0 / (mu e) [m^2/s]."""
    return gr.permeability(D, E0) * pc.P0 / (MU * E0)


def supply():
    """Superficial air velocity at the reservoir floor that carries P_FRAC of
    the weight through the packing (Darcy) [m/s]."""
    return gr.permeability(D, E0) / MU * P_FRAC * RHO_B * G_ACC


def describe():
    return (f"Dam break: column {X0} x {H0} m, {NX} x {NY} cells of {CELL * 1e3:g} mm; glass beads "
            f"d {D * 1e6:g} um, rho_p {RHOP:g}, e {E0}, rho_b {RHO_B:g} kg/m^3; H0^2/c "
            f"{diffusion_time():.3f} s; sqrt(h0/g) {T_SCALE:.4f} s; Roche: front {U_FRONT:.3f} m/s, "
            f"run-out {RUNOUT} h0 from the gate; floor {FLOOR:g} deg, initial p {P_FRAC} of the weight, "
            f"reservoir air {supply() * 1e3:.2f} mm/s")


def make(path, *, gas=True, decomp=(1, 1, 1), end=T_END, extra_gas="", cap=True, kappa=KAPPA_FIT,
         friction_only=None, slow=None):
    """kappa > 0: route B (model poreGas, diffusivity kappa); kappa = 0: the
    compressible gas of F.8c (model gasFlow, surfaceCarry yes).
    friction_only (route B only, F.8g): the pore pressure only lowers the
    friction (gasProperties frictionOnly); sigma starts as the total stress,
    isotropic -rho_b g (H0 - y). slow (default: with the cap) adds SLOW_FLOW
    (WP F.8k); cap=False is the v1 material with the 19 deg floor."""
    if slow is None:
        slow = cap
    floor = FLOOR if cap else FLOOR_V1
    if friction_only is None:
        friction_only = FRICTION_ONLY and kappa > 0
    if friction_only:
        extra_gas = "frictionOnly yes;\n" + extra_gas
    if os.path.isdir(path):
        shutil.rmtree(path)
    pc.make_case(path, size=(NX * CELL, NY * CELL), cells=(NX, NY), e=E0, d=D, end=end,
                 write_interval=WRITE, p_patches={"top": "type fixedValue; value uniform 0;"},
                 decomp=decomp, E=E_MOD, coupled=True)
    g = f"{path}/constant/gasProperties"
    pc.write(g, open(g).read().replace(f"rhoP    {pc.RHOP!r};", f"rhoP    {RHOP!r};"))
    inflow = f"inflow {{ floorRes {supply()!r}; }}\n" if gas else ""
    route = f"diffusivity {kappa!r};\n" if kappa > 0 else "surfaceCarry yes;\n"
    fl._common(path, rho=RHO_B, E=E_MOD, extra_gas=route + inflow + extra_gas,
               material=(MC if cap else MC.split("capLambda")[0].rstrip("\n"))
               + ("\n" + SLOW_FLOW if slow else ""),
               walls=f"    left {{ friction 0; }}\n    right {{ friction 0; }}\n"
                     f"    bottom {{ friction {floor!r}; }}"
                     + (f"\n    floorRes {{ friction {floor!r}; }}" if gas else ""))
    if gas and kappa > 0:
        pc.write(g, open(g).read().replace("model   gasFlow;", "model   poreGas;"))
    k0 = 0.3 / 0.7
    alpha, rho, sig, p = [], [], [], []
    for j in range(NY):
        y = (j + 0.5) * CELL
        for i in range(NX):
            bed = (i + 0.5) * CELL < X0 and y < H0
            alpha.append(1.0 if bed else 0.0)
            rho.append(RHO_B if bed else 0.0)
            w = RHO_B * G_ACC * (H0 - y) if bed else 0.0
            ws = (1.0 - P_FRAC) * w if gas else w
            if gas and friction_only:
                sig.append((-w, -w))
            else:
                sig.append((-k0 * ws, -ws))
            p.append(P_FRAC * w if gas else 0.0)
    fl._write_fields(path, alpha, rho, sig, N)
    txt = open(f"{path}/0/p").read()
    lst = f"nonuniform List<scalar> {N}\n(\n" + "\n".join(repr(v) for v in p) + "\n)"
    txt = txt.replace("internalField uniform 0;", f"internalField {lst};")
    if not gas:
        pc.write(f"{path}/0/p", txt)
        os.remove(g)
        return
    # the reservoir floor blows air (patch floorRes, made in run())
    pc.write(f"{path}/0/p", txt.replace(
        "    frontAndBack { type empty; }",
        "    floorRes { type fixedGradient; gradient uniform 0; }\n    frontAndBack { type empty; }"))
    pc.write(f"{path}/system/topoSetDict", pc.HEADER.format(cls="dictionary", obj="topoSetDict") + f"""
actions
(
    {{ name floorRes; type faceSet; action new; source patchToFace; patch bottom; }}
    {{ name floorRes; type faceSet; action subset; source boxToFace; box (-1 -1 -1) ({X0!r} 0.001 1); }}
);
""")
    pc.write(f"{path}/system/createPatchDict", pc.HEADER.format(cls="dictionary", obj="createPatchDict") + """
pointSync false;
patches ( { name floorRes; patchInfo { type patch; } constructFrom set; set floorRes; } );
""")


def _latest(path):
    ts = pc.times(path)
    return ts[-1] if ts else None


def run(path, nProcs=1):
    """Mesh, decompose and run, or resume from the last written time if the
    case already ran part way (checkpoint every WRITE s). Returns the wall
    time of the celFoam call(s) [s], summed over resumes (log.wallTime)."""
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    cd = f"{path}/system/controlDict"
    end = float(re.search(r"endTime\s+(\S+);", open(cd).read()).group(1))
    wt = f"{path}/log.wallTime"
    done = float(open(wt).read()) if os.path.exists(wt) else 0.0
    if nProcs > 1:
        pdir = f"{path}/processor0"
        last = _latest(pdir) if os.path.isdir(pdir) else None
    else:
        last = _latest(path) if os.path.isdir(f"{path}/constant/polyMesh") else None
    if last is not None and float(last) >= end - 1e-9:
        if nProcs > 1 and _latest(path) != last:
            pc.sh(path, ["reconstructPar"], "reconstructPar")
        return done
    if last is None:
        pc.sh(path, ["blockMesh"], "blockMesh")
        if os.path.exists(f"{path}/system/createPatchDict"):
            pc.sh(path, ["topoSet"], "topoSet")
            pc.sh(path, ["createPatch", "-overwrite"], "createPatch")
        if nProcs > 1:
            pc.sh(path, ["decomposePar", "-force"], "decomposePar")
        log = "celFoam"
    else:
        txt = open(cd).read().replace("startFrom       startTime;", "startFrom       latestTime;")
        open(cd, "w").write(txt)
        k = len([f for f in os.listdir(path) if f.startswith("log.celFoam")])
        log = f"celFoam.{k}"
    cmd = ["celFoam"] if nProcs == 1 else \
        ["mpirun", "-np", str(nProcs), "--oversubscribe", "celFoam", "-parallel"]
    t0 = time.perf_counter()
    pc.sh(path, cmd, log, env)
    done += time.perf_counter() - t0
    open(wt, "w").write(repr(done))
    if nProcs > 1:
        pc.sh(path, ["reconstructPar"], "reconstructPar")
    return done


def logs(path):
    """celFoam logs in run order (the first run, then each resume)."""
    fs = [f for f in os.listdir(path) if re.fullmatch(r"log\.celFoam(\.\d+)?", f)]
    return [f"{path}/{f}" for f in sorted(fs, key=lambda f: int(f.split(".")[2]) if f.count(".") == 2 else 0)]


def columns(path, t):
    """Height of material in each column (sum alpha dy) and its mass in
    cells' worth (sum rho / rho_b)."""
    al = pc.field(path, t, "alpha", N).reshape(NY, NX)
    rho = pc.field(path, t, "rho", N).reshape(NY, NX)
    return al.sum(axis=0) * CELL, rho.sum(axis=0) / RHO_B


def front(path, t, thresh=0.5):
    """Front: right edge of the last column holding at least thresh of a
    cell's mass [m] (as the dry collapse test, decision 0025)."""
    _, m = columns(path, t)
    f = np.where(m >= thresh)[0]
    return (f.max() + 1) * CELL if f.size else 0.0


def deposit_history(path):
    """(t, front) at each write with the front at the end of the continuous
    deposit from the gate (first column with less than half a cell of
    material): grains sprayed ahead of the flow do not count (F.8e) [m]."""
    rows = []
    for t in pc.times(path):
        hc, _ = columns(path, t)
        k = int(round(X0 / CELL))
        while k < NX and hc[k] >= 0.5 * CELL:
            k += 1
        rows.append((float(t), k * CELL))
    return np.array(rows)


def history(path, thresh=0.5):
    """Rows (t, front x) at every write, with t = 0."""
    rows = [(0.0, X0)]
    for t in pc.times(path):
        rows.append((float(t), front(path, t, thresh)))
    return np.array(rows)


def front_speed(h):
    """Least-squares slope of the front over Roche's constant-velocity phase
    t/sqrt(h0/g) in PHASE [m/s]."""
    tau = h[:, 0] / T_SCALE
    sel = (tau >= PHASE[0]) & (tau <= PHASE[1])
    return float(np.polyfit(h[sel, 0], h[sel, 1], 1)[0])


def stop_time(h, tol=CELL):
    """First write after which the front moves less than one cell [s]."""
    xf = h[-1, 1]
    for t, x in h:
        if xf - x < tol:
            return float(t)
    return float(h[-1, 0])


def gas_mass_drift(path):
    """Largest |gas mass - net inflow - its start| over the first gas mass,
    for each run and resume of the case apart (netIn restarts at 0 on a
    resume, decision 0129)."""
    worst = 0.0
    for lg in logs(path):
        rows = re.findall(r"Pore gas: \d+ solves \(\d+ Newton, \d+ iterations\), change \S+ Pa, "
                          r"mass (\S+) kg, net in (\S+) kg", open(lg).read())
        if not rows:
            continue
        a = np.array(rows, dtype=float)
        worst = max(worst, float(np.abs((a[:, 0] - a[:, 1]) - (a[0, 0] - a[0, 1])).max() / a[0, 0]))
    return worst


def steps(path):
    return sum(len(re.findall(r"^Time = \S+  deltaT = ", open(lg).read(), re.M)) for lg in logs(path))


if __name__ == "__main__":
    print(f"glass beads d {D * 1e6:.0f} um, e {E0}, rho_b {RHO_B:.0f} kg/m^3; "
          f"k {gr.permeability(D, E0):.3e} m^2; H0^2/c {diffusion_time():.3f} s; "
          f"sqrt(h0/g) {T_SCALE:.4f} s; Roche: U {U_FRONT:.3f} m/s, run-out {RUNOUT * H0:.3f} m "
          f"from the gate, stop at {6.5 * T_SCALE:.3f} s; U_mf Wen-Yu "
          f"{gr.wen_yu_umf(D, RHOP, RG, MU):.4f} m/s")
