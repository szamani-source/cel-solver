"""Gas with a free surface, void and an outlet (WP F.4, decisions 0125, 0129).

Shared by tests/gas/venting_checks.py (done-when test).

Venting column (make_column): fine quartz powder (d 0.1 mm, gas fraction
0.40, linear elastic E 100 MPa, nu 0.3) fills the lower NB cells of a column
one cell wide, NV void cells above it, smooth sealed walls and base, the top
patch at p0 (fixedValue 0 gauge). Gravity on, the bed starts at rest in its
geostatic effective stress (K0 sideways). The gas in the bed starts at
p_i = 100 Pa over p0, the void at p0, and the bed at rest under it: its
effective stress is the geostatic one less p_i (the uplift of the pressure
jump at the surface), K0 sideways. The gas vents through the free surface
(the face between the bed and the void), which holds p0. The total stress is
constant, so the pore pressure obeys Terzaghi's equation with
c = k / (mu (e/p0 + 1/M)) and the degree of venting
U = 1 - mean(p)/p_i follows the series of tools/gas_reference.py, as the
consolidation of WP F.3 does. p_i is below the weight of the top cell
(rho g dy = 370 Pa), so the gas does not lift the surface.

Outlet box (make_outlet): a box of the same powder with void above and an
outlet in the middle of the floor (materialProperties outlets) that opens
at OPEN_TIME. Mohr-Coulomb (phi 30, c 0, plain, no extra rules), coarser
grains (d 1 mm) so the gas lets the material flow. The gas starts at p0;
the outlet is fixedValue 0 in 0/p, sealed until it opens, then the gas
leaves through it and with the falling material.

Drag law ergunWenYu (Sia, 2026-11-24).
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import gas_reference as gr  # noqa: E402
import pore_cases as pc  # noqa: E402

E_GAS = 0.40
D_GRAIN = 1.0e-4
# Stiff, so the surface settles by about 1e-5 of the top cell. The soft
# bed of WP F.3 (E 1 MPa, make_column(E=1e6)) settles by 1e-3 of it; its
# surface cell then creeps (v1 partly filled cell, also without gas, known
# issue in docs/STATUS.md), which the venting check reports too
E_MOD, NU = 1.0e8, 0.3
RHO = (1 - E_GAS) * pc.RHOP
G_ACC = 9.81
DAMPING = 50.0

# Column: NB bed cells, NV void cells; with 28 cells, 4 processors cut at
# 7, 14 and 21 cells: through the bed and exactly at the free surface
NB, NV = 21, 7
HB = 0.5
DY = HB / NB
P_I = 100.0

MU = gr.air_viscosity(pc.T)
M_MOD = E_MOD * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
C_V = gr.diffusivity_deformable(D_GRAIN, E_GAS, MU, M_MOD, pc.P0)


def c_of(E):
    M = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    return gr.diffusivity_deformable(D_GRAIN, E_GAS, MU, M, pc.P0)


def t_of(T, E=E_MOD):
    return T * HB * HB / c_of(E)


def _fields(path, alpha, rho, p, sigma, n):
    def lst(cls, vals):
        return f"nonuniform List<{cls}> {n}\n(\n" + "\n".join(vals) + "\n)"

    common = 'boundaryField { ".*" { type zeroGradient; } frontAndBack { type empty; } }'
    pc.write(f"{path}/0/alpha", pc.HEADER.format(cls="volScalarField", obj="alpha")
             + f"dimensions [0 0 0 0 0 0 0];\ninternalField {lst('scalar', [repr(a) for a in alpha])};\n{common}\n")
    pc.write(f"{path}/0/rho", pc.HEADER.format(cls="volScalarField", obj="rho")
             + f"dimensions [1 -3 0 0 0 0 0];\ninternalField {lst('scalar', [repr(r) for r in rho])};\n{common}\n")
    if sigma is not None:
        rows = [f"({a!r} 0 0 {b!r} 0 {a!r})" for a, b in sigma]
        pc.write(f"{path}/0/sigma", pc.HEADER.format(cls="volSymmTensorField", obj="sigma")
                 + f"dimensions [1 -1 -2 0 0 0 0];\ninternalField {lst('symmTensor', rows)};\n{common}\n")
    txt = open(f"{path}/0/p").read()
    txt = txt.replace("internalField uniform 0;",
                      f"internalField {lst('scalar', [repr(v) for v in p])};")
    pc.write(f"{path}/0/p", txt)


def make_column(path, *, end, write_interval, decomp=(1, 1, 1), gas=True, p_i=P_I,
                nx=1, E=E_MOD):
    """Venting column (see the module doc), nx x (NB + NV) square cells."""
    ny = NB + NV
    pc.make_case(path, size=(nx * DY, ny * DY), cells=(nx, ny), e=E_GAS, d=D_GRAIN,
                 end=end, write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;"},
                 decomp=decomp, E=E, coupled=True)
    with open(f"{path}/system/controlDict", "a") as f:
        f.write(f"damping {DAMPING!r};\n")
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO!r};
freeSurface yes;
type linearElastic;
E {E!r};
nu {NU!r};
walls
{{
    left   {{ friction 0; }}
    right  {{ friction 0; }}
    bottom {{ friction 0; }}
}}
""")
    pc.write(f"{path}/constant/g", pc.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G_ACC!r} 0);\n")
    k0 = NU / (1 - NU)
    alpha, rho, p, sig = [], [], [], []
    for j in range(ny):
        y = (j + 0.5) * DY
        bed = j < NB
        for _ in range(nx):
            alpha.append(1.0 if bed else 0.0)
            rho.append(RHO if bed else 0.0)
            p.append(p_i if bed else 0.0)
            # At rest with the gas: the effective stress carries the
            # weight less the uplift p_i of the pressure jump at the surface
            s = -RHO * G_ACC * (HB - y) + p_i if bed else 0.0
            sig.append((k0 * s, s))
    _fields(path, alpha, rho, p, sig, nx * ny)
    if not gas:
        os.remove(f"{path}/constant/gasProperties")


def degree(path, t, nx=1):
    """U = 1 - mean(p in the bed)/p_i at write time t."""
    p = pc.field(path, t, "p", nx * (NB + NV))
    return 1.0 - float(np.mean(p[: nx * NB])) / P_I


# Outlet box
OUT_CELL = 0.01
OUT_NX, OUT_NY, OUT_FILL = 20, 24, 18   # cells; material in the lower 18 rows
OUT_N = 4                               # outlet cells, centred in the floor
OUT_D = 1.0e-3
OUT_RHO = 1500.0
OPEN_TIME = 0.02
MATERIALS = {
    "mohrCoulomb": "type mohrCoulomb;\nE 1e6;\nnu 0.3;\nphi 30;\nc 0;\npsi 0;",
    "linearElastic": "type linearElastic;\nE 1e6;\nnu 0.3;",
}


def make_outlet(path, *, end, write_interval, decomp=(1, 1, 1), gas=True,
                material="mohrCoulomb", refine=1, extra_gas=""):
    """Outlet box (see the module doc); refine divides the cells. Run
    with set_fields=True."""
    L = OUT_CELL / refine
    nx, ny = OUT_NX * refine, OUT_NY * refine
    n_out, fill = OUT_N * refine, OUT_FILL * refine
    ns = (nx - n_out) // 2
    xs = [0.0, ns * L, (ns + n_out) * L, nx * L]
    Hd = ny * L
    verts = ([f"({x!r} 0 0)" for x in xs] + [f"({x!r} {Hd!r} 0)" for x in xs]
             + [f"({x!r} 0 {L!r})" for x in xs] + [f"({x!r} {Hd!r} {L!r})" for x in xs])
    ncell = [ns, n_out, ns]
    blocks = "\n".join(f"    hex ({b} {b + 1} {5 + b} {4 + b} {8 + b} {9 + b} {13 + b} {12 + b}) "
                       f"({ncell[b]} {ny} 1) simpleGrading (1 1 1)" for b in range(3))

    def bot(b):
        return f"({b} {b + 1} {9 + b} {8 + b})"

    def top(b):
        return f"({4 + b} {12 + b} {13 + b} {5 + b})"

    pc.write(f"{path}/system/blockMeshDict", pc.HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices ({' '.join(verts)});
blocks
(
{blocks}
);
boundary
(
    left {{ type patch; faces ((0 8 12 4)); }}
    right {{ type patch; faces ((3 7 15 11)); }}
    bottom {{ type patch; faces ({bot(0)} {bot(2)}); }}
    outlet {{ type patch; faces ({bot(1)}); }}
    top {{ type patch; faces ({top(0)} {top(1)} {top(2)}); }}
    frontAndBack {{ type empty; faces (
        {' '.join(f'({b} {4 + b} {5 + b} {b + 1})' for b in range(3))}
        {' '.join(f'({8 + b} {9 + b} {13 + b} {12 + b})' for b in range(3))}); }}
);
""")
    # The rest of the case from pore_cases, then the fields and materials
    os.makedirs(f"{path}/system", exist_ok=True)
    keep = open(f"{path}/system/blockMeshDict").read()
    pc.make_case(path, size=(nx * L, Hd), cells=(nx, ny), e=1 - OUT_RHO / pc.RHOP, d=OUT_D,
                 end=end, write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;",
                            "outlet": "type fixedValue; value uniform 0;"},
                 decomp=decomp, coupled=True, extra_gas=extra_gas)
    pc.write(f"{path}/system/blockMeshDict", keep)
    with open(f"{path}/system/controlDict", "a") as f:
        f.write("damping 0;\n")
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {OUT_RHO!r};
freeSurface yes;
{MATERIALS[material]}
walls
{{
    left   {{ friction 20; }}
    right  {{ friction 20; }}
    bottom {{ friction 20; }}
}}
outlets
{{
    outlet {{ friction 20; openTime {OPEN_TIME!r}; }}
}}
""")
    pc.write(f"{path}/constant/g", pc.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G_ACC!r} 0);\n")
    # Fill by setFields (three blocks: the cells are not in rows); run
    # with pc.run(path, n, set_fields=True)
    common = 'boundaryField { ".*" { type zeroGradient; } frontAndBack { type empty; } }'
    for f, dim in (("alpha", "[0 0 0 0 0 0 0]"), ("rho", "[1 -3 0 0 0 0 0]")):
        pc.write(f"{path}/0/{f}", pc.HEADER.format(cls="volScalarField", obj=f)
                 + f"dimensions {dim};\ninternalField uniform 0;\n{common}\n")
    pc.write(f"{path}/system/setFieldsDict", pc.HEADER.format(cls="dictionary", obj="setFieldsDict") + f"""
defaultFieldValues ();
regions ( boxToCell {{ box (-1 -1 -1) (1 {fill * L!r} 1);
    fieldValues (volScalarFieldValue alpha 1 volScalarFieldValue rho {OUT_RHO!r}); }} );
""")
    if not gas:
        os.remove(f"{path}/constant/gasProperties")


if __name__ == "__main__":
    print(f"column: M {M_MOD:.6g} Pa, c {C_V:.6g} m^2/s, T = 1 at {t_of(1):.4g} s, "
          f"top-cell weight {RHO * G_ACC * DY:.1f} Pa, p_i {P_I} Pa")
