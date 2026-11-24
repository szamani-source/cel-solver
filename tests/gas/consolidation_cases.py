"""Gas-filled column loaded at the top (WP F.3, decisions 0125, 0128).

Shared by tests/gas/consolidation_checks.py (done-when test) and the gas
validation pack test 3, tests/validation/gas/terzaghiColumn.

A column of fine quartz powder (grains d, density 2650 kg/m^3, gas
fraction e) on a sealed, smooth base between smooth sealed side walls, so
it strains in one dimension (oedometer); the top is drained (p = 0 gauge)
and carries a load q from t = 0. Linear elastic skeleton (E, nu), so its
constrained modulus M = E (1 - nu) / ((1 + nu)(1 - 2 nu)). No gravity.

Linearised about p0 the pore pressure obeys Terzaghi's equation with
c = k / (mu (e/p0 + 1/M)) (decision 0125 item 4), k the Ergun
permeability; the sudden load starts it at the undrained pressure
p_i = q (1/M) / (e/p0 + 1/M), uniform. Degree of consolidation
U = 1 - mean(p)/p_i against the series of tools/gas_reference.py.

The load also starts elastic waves (period 4 H / c_u, about 0.05 s here)
that the quasi-static theory leaves out; the velocity damping D of
eulerianSolid takes them out. It adds a drag rho D u that delays the
consolidation by about rho D c / M (0.8 % at D = 50 1/s). Drag law
ergunWenYu (Sia, 2026-11-24).
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

E_GAS = 0.40            # gas fraction
D_GRAIN = 1.0e-4        # 0.1 mm grains
E_MOD, NU = 1.0e6, 0.3  # soft powder skeleton
H = 1.0                 # column height [m]
Q = 2000.0              # load [Pa]
DAMPING = 50.0          # [1/s]

MU = gr.air_viscosity(pc.T)
M_MOD = E_MOD * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
C_V = gr.diffusivity_deformable(D_GRAIN, E_GAS, MU, M_MOD, pc.P0)
P_I = Q * (1 / M_MOD) / (E_GAS / pc.P0 + 1 / M_MOD)
RHO = (1 - E_GAS) * pc.RHOP


def t_of(T):
    return T * H * H / C_V


def make(path, *, ny, end, write_interval, decomp=(1, 1, 1), gas=True,
         load=Q, at_rest=False, gravity=0.0, damping=DAMPING, extra="",
         nx=1):
    """Column of nx x ny square cells. at_rest: start in equilibrium with
    the load (and gravity), the effective stress set to the final state
    (nx = 1 only)."""
    dx = H / ny
    pc.make_case(path, size=(nx * dx, H), cells=(nx, ny), e=E_GAS, d=D_GRAIN,
                 end=end, write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;"},
                 decomp=decomp, E=E_MOD, coupled=True, extra_gas=extra)
    with open(f"{path}/system/controlDict", "a") as f:
        f.write(f"damping {damping!r};\n")
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO!r};
type linearElastic;
E {E_MOD!r};
nu {NU!r};
walls
{{
    left   {{ friction 0; }}
    right  {{ friction 0; }}
    bottom {{ friction 0; }}
}}
loads
{{
    top {{ pressure {load!r}; startTime 0; }}
}}
""")
    pc.write(f"{path}/constant/g", pc.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-gravity!r} 0);\n")
    if at_rest:
        # Final (drained) state: vertical effective stress -(q + rho g
        # depth), horizontal K0 times it (oedometer), at the cell centres
        k0 = NU / (1 - NU)
        y = (np.arange(ny) + 0.5) * dx
        syy = [float(s) for s in -(load + RHO * gravity * (H - y))]
        rows = "\n".join(f"({k0 * s!r} 0 0 {s!r} 0 {k0 * s!r})" for s in syy)
        pc.write(f"{path}/0/sigma", pc.HEADER.format(cls="volSymmTensorField", obj="sigma") + f"""
dimensions [1 -1 -2 0 0 0 0];
internalField nonuniform List<symmTensor> {ny}
(
{rows}
);
boundaryField {{ ".*" {{ type zeroGradient; }} frontAndBack {{ type empty; }} }}
""")
    if not gas:
        os.remove(f"{path}/constant/gasProperties")


def degree(path, t):
    """U = 1 - mean(p)/p_i at write time t (uniform cells)."""
    p = pc.field(path, t, "p")
    return 1.0 - float(np.mean(p)) / P_I


def theory(T):
    return gr.terzaghi_degree(T)


def profile_theory(ny, T):
    """p/p_i at the cell centres (z/H from the drained top)."""
    dx = H / ny
    y = (np.arange(ny) + 0.5) * dx
    return np.array([gr.terzaghi_pressure((H - yy) / H, T) for yy in y])


def deltaT(path):
    """Time steps of the run (from the log)."""
    import re
    txt = open(f"{path}/log.celFoam").read()
    return np.array([float(x) for x in re.findall(r"^Time = \S+\s+deltaT = (\S+)", txt, re.M)])


if __name__ == "__main__":
    print(f"M = {M_MOD:.6g} Pa, c = {C_V:.6g} m^2/s, p_i = {P_I:.6g} Pa, "
          f"t(T=1) = {t_of(1):.4g} s, rho D c / M = {RHO * DAMPING * C_V / M_MOD:.3g}")
