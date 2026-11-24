"""Air-impeded hopper discharge (WP F.5, decision 0130).

Shared by tests/gas/air_impeded_checks.py (done-when test) and the
validation pack test tests/validation/gas/fineHopper.

Plane-flow (slot) hopper, full width, on a body-fitted mesh (Sia,
2026-11-19: walls aligned with the mesh, never Cartesian cut-cell). Two
blockMesh blocks: the hopper is one trapezoid whose column lines fan from the
outlet edge to the top of the sloping walls (rows horizontal, graded fine
towards the outlet), the bin above has the same columns, vertical walls.

    bin half width R 0.1 m, outlet half width RO 0.01 m (20 mm slot),
    wall half angle BETA 30 deg from vertical, hopper height HT 0.156 m,
    bin 0.2 m, fill to 0.15 m above the transition.

Patches: outlet (floor, solid outlet that opens at OPEN_TIME, gas p0),
hopperWall and binWall (walls, friction PHI_W, no gas through them), top
(p0, void above the fill). The grains are quartz (rho_p 2650 kg/m^3)
packed at gas fraction 0.40 (bulk 1590 kg/m^3), Mohr-Coulomb phi 30, c 0,
psi 1 (Sia, 2026-11-19: low dilation), E 5 MPa. The grain diameter d is the
only change between cases: it sets the permeability (Ergun) and so how hard
the air resists the dilating flow. Gravity on from rest, velocity damping
while the outlet is closed (off when it opens, as in v1), fullBlend void.

Drag law ergunWenYu (Sia, 2026-11-24).
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import pore_cases as pc  # noqa: E402

R, RO, BETA = 0.10, 0.03, 30.0
HT = (R - RO) / math.tan(math.radians(BETA))
HC = 0.20
FILL = HT + 0.15
E_GAS = 0.40
RHO = (1 - E_GAS) * pc.RHOP
PHI_W = 20.0
G_ACC = 9.81
OPEN_TIME = 0.15
DAMPING = 50.0
MATERIAL = "type mohrCoulomb;\nE 1e6;\nnu 0.3;\nphi 30;\nc 0;\npsi 1;\nsigmaT 0;\nfullBlend yes;"
# v1's T0 rule set for flowing grains (decisions 0102, 0105), as in the Couto
# hopper; mu(I) with the grain diameter d of the case
T0_RULES = ("mu2 0.8412;\nI0 0.279;\nd {d!r};\nrhoS {rhoS!r};\nshearCap yes;\n"
            "spreadLayers yes;\nflowImpedance yes;\nnormalCap yes;\n")


NCOL, OUTLET_RATIO = 16, 0.25


def depth(ncol=NCOL):
    """Mesh thickness (one cell) [m]: rates per metre divide by it."""
    return 2 * R / ncol


def ncells(ncol=NCOL, outlet_ratio=OUTLET_RATIO):
    nh, nc = mesh_rows(ncol, outlet_ratio)
    return ncol * (nh + nc)


def outlet_row(ncol=NCOL, outlet_ratio=OUTLET_RATIO):
    """Rows of the hopper block and the outlet row's height [m]."""
    nh, _ = mesh_rows(ncol, outlet_ratio)
    r = (1 / outlet_ratio) ** (1 / (nh - 1))
    return nh, HT * (r - 1) / (r ** nh - 1)


def mesh_rows(ncol, outlet_ratio):
    """Rows of the hopper block graded from the outlet (height about
    outlet_ratio x the hopper's top row) and of the bin (square-ish)."""
    top = 2 * R / ncol                      # top cell width
    q = 1 / outlet_ratio

    def top_row(m):
        r = q ** (1 / max(m - 1, 1))
        return HT * (r - 1) / (r ** m - 1) * r ** (m - 1)
    nh = next(m for m in range(2, 400) if top_row(m) <= top * 1.0001)
    nc = int(round(HC / top))
    return nh, nc


def write_mesh(path, ncol, outlet_ratio):
    nh, nc = mesh_rows(ncol, outlet_ratio)
    dz = 2 * R / ncol
    v2 = [(-RO, 0), (RO, 0), (R, HT), (-R, HT), (R, HT + HC), (-R, HT + HC)]
    v = [(x, y, 0.0) for x, y in v2] + [(x, y, dz) for x, y in v2]
    vs = "\n".join(f"    ({x!r} {y!r} {z!r})" for x, y, z in v)
    pc.write(f"{path}/system/blockMeshDict", pc.HEADER.format(cls="dictionary", obj="blockMeshDict") + f"""
scale 1;
vertices
(
{vs}
);
blocks
(
    hex (0 1 2 3 6 7 8 9) ({ncol} {nh} 1) simpleGrading (1 {1 / outlet_ratio!r} 1)
    hex (3 2 4 5 9 8 10 11) ({ncol} {nc} 1) simpleGrading (1 1 1)
);
boundary
(
    outlet {{ type patch; faces ((0 6 7 1)); }}
    hopperWall {{ type patch; faces ((1 7 8 2) (0 3 9 6)); }}
    binWall {{ type patch; faces ((2 8 10 4) (3 5 11 9)); }}
    top {{ type patch; faces ((5 4 10 11)); }}
    frontAndBack {{ type empty; faces ((0 1 2 3) (3 2 4 5) (6 9 8 7) (9 11 10 8)); }}
);
""")
    return ncol * (nh + nc), nh, nc


def make_hopper(path, *, d, end, write_interval, gas=True, decomp=(1, 1, 1), ncol=NCOL,
                outlet_ratio=OUTLET_RATIO, extra_gas="", open_time=OPEN_TIME, rules=True):
    """Hopper case with grains of diameter d [m]; gas False drops
    constant/gasProperties (dry run, same mesh). Run with set_fields=True.
    Returns the cell count."""
    pc.make_case(path, size=(2 * R, HT + HC), cells=(ncol, 10), e=E_GAS, d=d,
                 end=end, write_interval=write_interval,
                 p_patches={"top": "type fixedValue; value uniform 0;",
                            "outlet": "type fixedValue; value uniform 0;"},
                 decomp=decomp, coupled=True, extra_gas=extra_gas)
    ncell, _, _ = write_mesh(path, ncol, outlet_ratio)
    with open(f"{path}/system/controlDict", "a") as f:
        f.write(f"damping {DAMPING!r};\n")
    pc.write(f"{path}/constant/materialProperties",
             pc.HEADER.format(cls="dictionary", obj="materialProperties") + f"""
rho {RHO!r};
freeSurface yes;
{MATERIAL}
{T0_RULES.format(d=d, rhoS=pc.RHOP) if rules else ""}walls
{{
    hopperWall {{ friction {PHI_W!r}; }}
    binWall    {{ friction {PHI_W!r}; }}
}}
outlets
{{
    outlet {{ friction {PHI_W!r}; openTime {open_time!r}; }}
}}
""")
    pc.write(f"{path}/constant/g", pc.HEADER.format(cls="uniformDimensionedVectorField", obj="g")
             + f"dimensions [0 1 -2 0 0 0 0];\nvalue (0 {-G_ACC!r} 0);\n")
    common = 'boundaryField { ".*" { type zeroGradient; } frontAndBack { type empty; } }'
    for f, dim in (("alpha", "[0 0 0 0 0 0 0]"), ("rho", "[1 -3 0 0 0 0 0]")):
        pc.write(f"{path}/0/{f}", pc.HEADER.format(cls="volScalarField", obj=f)
                 + f"dimensions {dim};\ninternalField uniform 0;\n{common}\n")
    pc.write(f"{path}/system/setFieldsDict", pc.HEADER.format(cls="dictionary", obj="setFieldsDict") + f"""
defaultFieldValues ();
regions ( boxToCell {{ box (-1 -1 -1) (1 {FILL!r} 1);
    fieldValues (volScalarFieldValue alpha 1 volScalarFieldValue rho {RHO!r}); }} );
""")
    if not gas:
        os.remove(f"{path}/constant/gasProperties")
    return ncell


def beverloo_slot(d, C=1.03, k=1.5):
    """Plane-flow Beverloo rate per metre of slot [kg/(m s)]:
    W = C rho_b sqrt(g) (B - k d)^1.5 (B the slot width). Reference scale
    only; v1's own dry rate is the control."""
    return C * RHO * math.sqrt(G_ACC) * (2 * RO - k * d) ** 1.5


if __name__ == "__main__":
    nh, nc = mesh_rows(NCOL, OUTLET_RATIO)
    print(f"hopper height {HT:.4f} m, rows {nh} + {nc}, fill to {FILL:.4f} m, "
          f"bulk {RHO} kg/m^3")


def draw_mesh(png, ncol=NCOL, outlet_ratio=OUTLET_RATIO, title=""):
    """Picture of the mesh (cell lines, fill level, outlet) for Sia."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nh, nc = mesh_rows(ncol, outlet_ratio)
    r = (1 / outlet_ratio) ** (1 / (nh - 1))
    ys = [0.0]
    for i in range(nh):
        ys.append(ys[-1] + HT * (r - 1) / (r ** nh - 1) * r ** i)
    ys += [HT + HC * (j + 1) / nc for j in range(nc)]

    def half(y):
        return RO + (R - RO) * min(y, HT) / HT
    fig, axs = plt.subplots(1, 2, figsize=(10, 6.5), gridspec_kw=dict(width_ratios=[1.2, 1]))
    for k, ax in enumerate(axs):
        ax.fill([-RO, RO, R, R, -R, -R], [0, 0, HT, FILL, FILL, HT], color="#e8d5a8", zorder=0)
        for y in ys:
            ax.plot([-half(y), half(y)], [y, y], color="#444", lw=0.5)
        for i in range(ncol + 1):
            s = -1 + 2 * i / ncol
            ax.plot([s * RO, s * R, s * R], [0, HT, HT + HC], color="#444", lw=0.5)
        ax.plot([-RO, RO], [0, 0], color="#d62728", lw=3, label="outlet (opens at 0.15 s)")
        ax.set_aspect("equal")
        if k == 1:
            ax.set_xlim(-0.06, 0.06)
            ax.set_ylim(-0.005, 0.07)
            ax.set_title(f"outlet close up: {ncol} cells across {2 * RO * 1000:g} mm, "
                         f"{2 * RO / ncol * 1000:.2f} x {ys[1] * 1000:.1f} mm", fontsize=9)
        else:
            ax.set_title(title or f"{ncol * (nh + nc)} cells, sand colour = fill", fontsize=9)
            ax.legend(loc="upper right", fontsize=8)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
    fig.suptitle(f"Plane-flow hopper for WP F.5: slot {2 * RO * 1000:g} mm, walls "
                 f"{BETA:g} deg from vertical, bin {2 * R:g} m wide, fill {FILL:.3f} m", fontsize=10)
    fig.tight_layout()
    fig.savefig(png, dpi=150)
