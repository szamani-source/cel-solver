"""Writes 0/Ug and 0/p of heatedAirChannel: the developed laminar
(Poiseuille) profile at a mean speed of 0.53 m/s and the linear pressure
drop 12 mu u_m L/H^2 that drives it. Edit H, L, NX, NY if you change
system/blockMeshDict. Called by Allrun after blockMesh.
"""
H, L = 2.0e-3, 60.0e-3      # channel height and length [m]
NX, NY = 200, 16            # cells (as system/blockMeshDict)
U_M = 0.53                  # mean speed [m/s]
T = 293.15
MU = 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)   # Sutherland
DP = 12 * MU * U_M * L / H ** 2

u, p = [], []
for j in range(NY):
    y = (j + 0.5) * H / NY
    for i in range(NX):
        x = (i + 0.5) * L / NX
        u.append(DP / (2 * MU * L) * y * (H - y))
        p.append(DP * (1 - x / L))


def replace(name, new):
    txt = open(f"0/{name}").read()
    a = txt.index("internalField")
    b = txt.index("boundaryField")
    open(f"0/{name}", "w").write(txt[:a] + new + "\n" + txt[b:])


replace("Ug", f"internalField nonuniform List<vector> {len(u)} ("
        + " ".join(f"({v!r} 0 0)" for v in u) + ");")
replace("p", f"internalField nonuniform List<scalar> {len(p)} ("
        + " ".join(repr(v) for v in p) + ");")
print(f"Poiseuille profile, mean {U_M} m/s, pressure drop {DP:.6g} Pa")
