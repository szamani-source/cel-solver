"""Initial fields of the fluidisedDamBreak tutorial: 0/alpha, rho, p, sigma.

    python3 makeColumn.py [cell in mm] [dry]

The column 0.2 m wide and 0.4 m high stands against the left wall. Its pore
pressure carries the whole weight (Roche et al. 2008: fully fluidised),
p = rho_b g (0.4 - y); the stress starts as the total stress, isotropic
-rho_b g (0.4 - y) (with frictionOnly the pore pressure only lowers the
friction). dry: no pore pressure, stress at rest with K0 = nu/(1 - nu).
Uniform mesh of square cells, numbered row by row as blockMesh does.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
G, RHO_B, X0, H0, LX, LY = 9.81, 1450.0000000000002, 0.2, 0.4, 3.0, 0.45
HEAD = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def write(case, name, cls, dims, values, boundary):
    with open(os.path.join(case, "0", name), "w") as f:
        f.write(HEAD.format(cls=cls, obj=name))
        typ = "symmTensor" if cls == "volSymmTensorField" else "scalar"
        f.write(f"dimensions {dims};\ninternalField nonuniform List<{typ}> {len(values)}\n(\n")
        f.write("\n".join(values))
        f.write(f"\n);\nboundaryField\n{{\n{boundary}\n}}\n")


def main(cell, dry, case=HERE):
    nx, ny = int(round(LX / cell)), int(round(LY / cell))
    k0 = 0.3 / 0.7
    alpha, rho, p, sig = [], [], [], []
    for j in range(ny):
        y = (j + 0.5) * cell
        for i in range(nx):
            bed = (i + 0.5) * cell < X0 and y < H0
            w = RHO_B * G * (H0 - y) if bed else 0.0
            alpha.append("1" if bed else "0")
            rho.append(repr(RHO_B) if bed else "0")
            p.append("0" if dry else repr(w))
            sx, sy = (-k0 * w, -w) if dry else (-w, -w)
            sig.append(f"({sx!r} 0 0 {sy!r} 0 {sx!r})")
    os.makedirs(os.path.join(case, "0"), exist_ok=True)
    zg = '    ".*" { type zeroGradient; }\n    frontAndBack { type empty; }'
    write(case, "alpha", "volScalarField", "[0 0 0 0 0 0 0]", alpha, zg)
    write(case, "rho", "volScalarField", "[1 -3 0 0 0 0 0]", rho, zg)
    write(case, "sigma", "volSymmTensorField", "[1 -1 -2 0 0 0 0]", sig, zg)
    if not dry:
        write(case, "p", "volScalarField", "[1 -1 -2 0 0 0 0]", p,
              '    ".*" { type zeroGradient; }\n'
              '    top { type fixedValue; value uniform 0; }\n'
              '    floorRes { type fixedGradient; gradient uniform 0; }   // flux set by gasProperties inflow\n'
              '    frontAndBack { type empty; }')


if __name__ == "__main__":
    cell = float(sys.argv[1]) / 1000 if len(sys.argv) > 1 else 0.01
    main(cell, "dry" in sys.argv[2:], os.getcwd())
