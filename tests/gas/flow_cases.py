"""Gas flow in the void (WP F.6, decisions 0125, 0131).

Shared by tests/gas/gas_flow_checks.py (done-when test). Every case uses the
ergunWenYu drag law (Sia, 2026-11-24; no grains here, so it does not act)
and a freeSurface solid with no material: the domain is void, the gas
moves by its own momentum (model gasFlow). The solid's material only sets
the time step (P-wave speed of the reference density over the cell size).

Channel (make_channel): plane channel H x L, NY cells across, NX along, no
slip walls (top, bottom), the inlet (left) at DP over p0 and the outlet
(right) at p0, Ug zeroGradient there (fully developed). Starts at rest;
the steady flow is Poiseuille's: u(y) = DP/(2 mu L) y (H - y), mean
DP H^2/(12 mu L). The slowest mode decays as exp(-pi^2 nu t/H^2), so END
(about 1.2 H^2/nu) leaves 1e-5 of it.

Box (make_box): closed box, non-orthogonal cells (top edge shifted), the gas
at p0 and at rest: it must stay at rest.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

MU = gr.air_viscosity(pc.T)
RHO_G = pc.P0 / (287.05 * pc.T)
NU_G = MU / RHO_G

H, L = 2.0e-3, 2.0e-3
NY, NX = 20, 10
DP = 1.0e-3                       # [Pa], u_max 1.4 cm/s, Re 1.9
END = 1.2 * H * H / NU_G
U_MEAN = DP * H * H / (12 * MU * L)
U_MAX = 1.5 * U_MEAN

# Solid: no material in the domain; E sets the time step only
E_SOLID, RHO_SOLID = 1.0e4, 1590.0


def u_exact(y):
    return DP / (2 * MU * L) * y * (H - y)


def _void(path, model, extra_gas):
    common = 'boundaryField { ".*" { type zeroGradient; } frontAndBack { type empty; } }'
    for f, dim in (("alpha", "[0 0 0 0 0 0 0]"), ("rho", "[1 -3 0 0 0 0 0]")):
        pc.write(f"{path}/0/{f}", pc.HEADER.format(cls="volScalarField", obj=f)
                 + f"dimensions {dim};\ninternalField uniform 0;\n{common}\n")
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO_SOLID!r};
freeSurface yes;
type linearElastic;
E {E_SOLID!r};
nu 0.3;
""")
    txt = open(f"{path}/constant/gasProperties").read()
    txt = txt.replace("model   poreGas;", f"model   {model};")
    pc.write(f"{path}/constant/gasProperties", txt + extra_gas)


def make_channel(path, *, end=END, write_interval=END / 4, decomp=(1, 1, 1), gas=True,
                 model="gasFlow", extra_gas="pTolerance 1e-15;\n", nx=NX, ny=NY):
    pc.make_case(path, size=(L, H), cells=(nx, ny), e=0.4, d=1e-3, end=end,
                 write_interval=write_interval,
                 p_patches={"left": f"type fixedValue; value uniform {DP!r};",
                            "right": "type fixedValue; value uniform 0;"},
                 decomp=decomp, E=E_SOLID, coupled=True)
    _void(path, model, extra_gas)
    if not gas:
        os.remove(f"{path}/constant/gasProperties")


def make_box(path, *, end, write_interval, decomp=(1, 1, 1), skew=0.3e-3):
    pc.make_case(path, size=(L, H), cells=(NX, NY), e=0.4, d=1e-3, end=end,
                 write_interval=write_interval, p_patches={}, skew=skew,
                 decomp=decomp, E=E_SOLID, coupled=True)
    _void(path, "gasFlow", "pTolerance 1e-15;\n")


if __name__ == "__main__":
    print(f"mu {MU:.6g} Pa s, nu {NU_G:.6g} m^2/s, H^2/nu {H * H / NU_G:.4g} s, "
          f"end {END:.4g} s, u_max {U_MAX:.6g} m/s, Re {U_MEAN * H / NU_G:.3g}")
