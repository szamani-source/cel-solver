"""Fluidisation cases (WP F.7, decisions 0125, 0132).

Shared by tests/gas/fluidisation_checks.py (done-when test) and the gas
validation pack test 2, tests/validation/gas/fluidOnset. Every case uses the
ergunWenYu drag law (Sia, 2026-11-24), model gasFlow (gas in the pores and
the void) and plain Mohr-Coulomb grains (phi 30, c 0, psi 0, no extra rules:
the v1 T0 rules part serial from parallel, decision 0130).

Bed (make_bed): quartz sand (d 0.2 mm, rho_p 2650) packed at gas fraction
E_MF = 0.38 (Wen and Yu's correlation stands for e_mf near 0.38: its
constants are 1/(phi_s e^3) = 14 and (1-e)/(phi_s^2 e^3) = 11), HB = 5 cm
deep in a column NX cells wide, NB bed rows and NV void rows above (2.5 mm
cells), smooth walls, the top open (p fixedValue 0). Air is blown in through
the floor (fixedGradient p, gasProperties inflow) at a superficial velocity
ramped from 0 to U_END over T_RAMP. The bed starts at rest in its geostatic
stress. The bed lifts when its pressure drop reaches its weight per area,
rho g HB, at the minimum fluidisation velocity U_mf.

Slab (make_slab): a slab of fine sand (d 50 um) suspended at gas fraction
E_SLAB = 0.9, HS thick, starts at rest in the middle of a closed column
(sealed floor, open top) with void above and below. It falls; the air
under it must pass through it, so it settles at the hindered settling speed,
Richardson and Zaki's U_t e^n.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402
import venting_cases as vc  # noqa: E402

MU = gr.air_viscosity(pc.T)
RG = gr.air_density(pc.P0, pc.T)
G_ACC = 9.81
MC = "type mohrCoulomb;\nE {E!r};\nnu 0.3;\nphi 30;\nc 0;\npsi 0;"

# Bed
D_BED = 2.0e-4
E_MF = 0.38
RHO_BED = (1 - E_MF) * pc.RHOP
HB, NB, NV, NX = 0.05, 20, 20, 4
DX = HB / NB
E_BED = 1.0e7
W_BED = RHO_BED * G_ACC * HB              # weight per area [Pa]
U_WY = gr.wen_yu_umf(D_BED, pc.RHOP, RG, MU)
U_MODEL = gr.model_umf(D_BED, pc.RHOP, RG, MU, E_MF)
U_END = 1.4 * U_WY
T_RAMP = 0.28


def _write_fields(path, alpha, rho, sig, n):
    def lst(cls, vals):
        return f"nonuniform List<{cls}> {n}\n(\n" + "\n".join(vals) + "\n)"

    common = 'boundaryField { ".*" { type zeroGradient; } frontAndBack { type empty; } }'
    pc.write(f"{path}/0/alpha", pc.HEADER.format(cls="volScalarField", obj="alpha")
             + f"dimensions [0 0 0 0 0 0 0];\ninternalField {lst('scalar', [repr(a) for a in alpha])};\n{common}\n")
    pc.write(f"{path}/0/rho", pc.HEADER.format(cls="volScalarField", obj="rho")
             + f"dimensions [1 -3 0 0 0 0 0];\ninternalField {lst('scalar', [repr(r) for r in rho])};\n{common}\n")
    rows = [f"({a!r} 0 0 {b!r} 0 {a!r})" for a, b in sig]
    pc.write(f"{path}/0/sigma", pc.HEADER.format(cls="volSymmTensorField", obj="sigma")
             + f"dimensions [1 -1 -2 0 0 0 0];\ninternalField {lst('symmTensor', rows)};\n{common}\n")


def _common(path, *, rho, E, walls, extra_gas, material=None):
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {rho!r};
freeSurface yes;
{material or MC.format(E=E)}
walls
{{
{walls}
}}
""")
    pc.write(f"{path}/constant/g", pc.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G_ACC!r} 0);\n")
    g = f"{path}/constant/gasProperties"
    txt = open(g).read().replace("model   poreGas;", "model   gasFlow;")
    pc.write(g, txt + extra_gas)
    with open(f"{path}/system/controlDict", "a") as f:
        f.write("damping 0;\n")


def make_bed(path, *, end=T_RAMP, write_interval=T_RAMP / 56, decomp=(1, 1, 1), gas=True,
             u_end=U_END, t_ramp=T_RAMP, extra_gas="", nx=NX, material=None):
    ny = NB + NV
    pc.make_case(path, size=(nx * DX, ny * DX), cells=(nx, ny), e=E_MF, d=D_BED, end=end,
                 write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;",
                            "bottom": "type fixedGradient; gradient uniform 0;"},
                 inflow={"bottom": 0.0}, decomp=decomp, E=E_BED, coupled=True)
    g = f"{path}/constant/gasProperties"
    pc.write(g, open(g).read().replace("    bottom 0.0;", f"    bottom table ((0 0) ({t_ramp!r} {u_end!r}));"))
    _common(path, rho=RHO_BED, E=E_BED, extra_gas=extra_gas, material=material,
            walls="    left { friction 0; }\n    right { friction 0; }\n    bottom { friction 0; }")
    k0 = 0.3 / 0.7
    alpha, rho, sig = [], [], []
    for j in range(ny):
        y = (j + 0.5) * DX
        bed = j < NB
        for _ in range(nx):
            alpha.append(1.0 if bed else 0.0)
            rho.append(RHO_BED if bed else 0.0)
            s = -RHO_BED * G_ACC * (HB - y) if bed else 0.0
            sig.append((k0 * s, s))
    _write_fields(path, alpha, rho, sig, nx * ny)
    if not gas:
        os.remove(g)


def u_in(t, u_end=U_END, t_ramp=T_RAMP):
    return u_end * min(t / t_ramp, 1.0)


def bed_state(path, t, nx=NX):
    """Pressure drop over the bed [Pa] (the floor's face pressure), the
    bed's mean vertical velocity [m/s] and its centre of mass height [m]."""
    ny = NB + NV
    n = nx * ny
    p = pc.field(path, t, "p", n)
    pb = float(np.mean(p[:nx] + pc.patch_gradient(path, t, "p", "bottom") * DX / 2))
    U = pc.field(path, t, "U", n)
    rho = pc.field(path, t, "rho", n)
    y = (np.repeat(np.arange(ny), nx) + 0.5) * DX
    m = rho.sum()
    return pb, float((rho * U[:, 1]).sum() / m), float((rho * y).sum() / m)


def onset(path, nx=NX):
    """U_mf from the run: where the pressure drop before lift (Ergun-like,
    fitted u + u^2 through the origin on points below 0.85 W) meets the
    bed weight W; and where the bed starts to rise (its velocity after lift,
    fitted linear in U, extrapolated to 0). Returns (U_dp, U_lift, rows)."""
    rows = []
    for t in pc.times(path):
        pb, v, yc = bed_state(path, t, nx)
        rows.append((float(t), u_in(float(t)), pb, v, yc))
    a = np.array(rows)
    # Before lift: up to the first write at 0.85 W (after lift the
    # pressure drop of a wide bed swings as it breaks up)
    k = int(np.argmax(a[:, 2] >= 0.85 * W_BED)) if (a[:, 2] >= 0.85 * W_BED).any() else len(a)
    lo = a[:k]
    lo = lo[lo[:, 1] > 0]
    A = np.vstack([lo[:, 1], lo[:, 1] ** 2]).T
    c1, c2 = np.linalg.lstsq(A, lo[:, 2], rcond=None)[0]
    u_dp = float((-c1 + np.sqrt(c1 * c1 + 4 * c2 * W_BED)) / (2 * c2)) if c2 != 0 else W_BED / c1
    # Rise speed just after lift (up to 1.15 U_mf)
    up = a[(a[:, 3] > 2e-3) & (a[:, 1] <= 1.15 * u_dp)]
    u_lift = float("nan")
    if len(up) >= 3:
        s, i = np.polyfit(up[:, 1], up[:, 3], 1)
        u_lift = -i / s
    return float(u_dp), u_lift, a


def compare_runs(a, b, n, t_end, fields=("p", "Ug", "U")):
    """Two runs of the same case (processor counts): the largest difference
    of each field up to t_end over its largest |value| in that window; the
    same per write time (each over that time's largest); and the largest
    per-time difference after t_end (report)."""
    import gas_flow_checks as gf
    before, per, after = {}, {}, {}
    for f in fields:
        A, B, d_after = [], [], 0.0
        for t in pc.times(a):
            x, y = gf.vfield(a, t, f, n), gf.vfield(b, t, f, n)
            if float(t) <= t_end:
                A.append(x)
                B.append(y)
            else:
                d_after = max(d_after, gf.rel(x, y))
        scale = max(float(np.abs(x).max()) for x in A)
        before[f] = max(float(np.abs(x - y).max()) for x, y in zip(A, B)) / scale
        per[f] = max(gf.rel(x, y) for x, y in zip(A, B))
        after[f] = d_after
    return before, per, after


# Slab
D_SLAB = 5.0e-5
E_SLAB = 0.9
RHO_SLAB = (1 - E_SLAB) * pc.RHOP
DY_S = 1.0e-3
NS_BELOW, NS_SLAB, NS_ABOVE = 30, 20, 20
NX_S = 2
E_SLAB_MOD = 1.0e6
T_SLAB = 0.08
V_MODEL = gr.model_settling(D_SLAB, pc.RHOP, RG, MU, E_SLAB)
V_RZ = gr.richardson_zaki(D_SLAB, pc.RHOP, RG, MU, E_SLAB)


def make_slab(path, *, end=T_SLAB, write_interval=T_SLAB / 16, decomp=(1, 1, 1), extra_gas=""):
    ny = NS_BELOW + NS_SLAB + NS_ABOVE
    pc.make_case(path, size=(NX_S * DY_S, ny * DY_S), cells=(NX_S, ny), e=E_SLAB, d=D_SLAB,
                 end=end, write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;"},
                 decomp=decomp, E=E_SLAB_MOD, coupled=True)
    _common(path, rho=RHO_SLAB, E=E_SLAB_MOD, extra_gas=extra_gas,
            walls="    left { friction 0; }\n    right { friction 0; }\n    bottom { friction 0; }")
    alpha, rho, sig = [], [], []
    for j in range(ny):
        slab = NS_BELOW <= j < NS_BELOW + NS_SLAB
        for _ in range(NX_S):
            alpha.append(1.0 if slab else 0.0)
            rho.append(RHO_SLAB if slab else 0.0)
            sig.append((0.0, 0.0))
    _write_fields(path, alpha, rho, sig, NX_S * ny)


def slab_state(path, t):
    """The slab's mean vertical velocity [m/s] and centre of mass [m]."""
    ny = NS_BELOW + NS_SLAB + NS_ABOVE
    n = NX_S * ny
    U = pc.field(path, t, "U", n)
    rho = pc.field(path, t, "rho", n)
    y = (np.repeat(np.arange(ny), NX_S) + 0.5) * DY_S
    m = rho.sum()
    return float((rho * U[:, 1]).sum() / m), float((rho * y).sum() / m)


# Relaxation of the pore gas momentum (frozen bed, coarse grains)
D_RLX, E_RLX = 5.0e-3, 0.4
L_RLX, N_RLX = 0.02, 20
DP_RLX = 2.0e-3                     # [Pa]
E_RLX_MOD = 1.0e4                   # sets the step (1.4e-4 s), solid frozen


def make_relax(path, *, end, write_interval, decomp=(1, 1, 1), extra_gas=""):
    """Column filled with a frozen bed (coupled no), the gas at rest, the
    floor at DP_RLX over the top from t = 0: the flux relaxes to Ergun's
    over tau = e rho_g / K."""
    pc.make_case(path, size=(L_RLX / N_RLX, L_RLX), cells=(1, N_RLX), e=E_RLX, d=D_RLX,
                 end=end, write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;",
                            "bottom": f"type fixedValue; value uniform {DP_RLX!r};"},
                 decomp=decomp, E=E_RLX_MOD, coupled=False)
    rho = (1 - E_RLX) * pc.RHOP
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {rho!r};
freeSurface yes;
type linearElastic;
E {E_RLX_MOD!r};
nu 0.3;
""")
    g = f"{path}/constant/gasProperties"
    pc.write(g, open(g).read().replace("model   poreGas;", "model   gasFlow;") + extra_gas)


def relax_reference(times, law=None):
    """Superficial flux q(t) of e rho dq/dt = -e^2 G - K(q/e) q (RK4),
    G = DP/L, and its steady value."""
    law = law or gr.drag_ergun_wen_yu
    e, G = E_RLX, DP_RLX / L_RLX

    def f(q):
        return (-e * e * (-G) - law(e, abs(q) / e, D_RLX, MU, RG) * q) / (e * RG)

    out, q, t, h = [], 0.0, 0.0, 1e-6
    for T in times:
        while t < T - 1e-15:
            hh = min(h, T - t)
            k1 = f(q)
            k2 = f(q + 0.5 * hh * k1)
            k3 = f(q + 0.5 * hh * k2)
            k4 = f(q + hh * k3)
            q += hh * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            t += hh
        out.append(q)
    q_inf = e * gr.model_slip(e, G, D_RLX, MU, RG, law)
    tau = e * RG / law(e, q_inf / e, D_RLX, MU, RG)
    return np.array(out), q_inf, tau


if __name__ == "__main__":
    print(f"bed: rho {RHO_BED:.1f} kg/m^3, weight {W_BED:.2f} Pa, U_mf Wen-Yu {U_WY:.5f} m/s, "
          f"model (Ergun e {E_MF}) {U_MODEL:.5f} m/s ({100 * (U_MODEL / U_WY - 1):+.2f} %), "
          f"U_end {U_END:.5f}")
    print(f"slab: rho {RHO_SLAB:.1f}, settling model {V_MODEL:.5f} m/s, Richardson-Zaki {V_RZ:.5f} m/s "
          f"({100 * (V_MODEL / V_RZ - 1):+.2f} %)")
