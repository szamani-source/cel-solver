"""WP 4.3c done-when checks: STL cut-cell walls (decision 0047) with the free
surface (void, decision 0031); the combination is decision 0051. Run by
tests/test_stl_void.sh after celVoidCheck and celSiloCheck are built:

    python3 stl_void_checks.py WORKDIR

Material T0 (cohesionless Mohr-Coulomb, rho 1800, E 10 MPa, nu 0.3, phi 35),
10 mm cells.

1. Lake at rest in an STL tank: STL side walls cutting cells at a quarter
   (open 0.75) and three quarters (open 0.25, a small cell mixed with its
   neighbour), mesh-wall floor (an STL floor makes corner cells with one
   oblique wall face, not in discrete equilibrium: decision 0051), 0.1 m wide,
   filled to 0.455 m (top cell half full), K0 stress of the weight above,
   full gravity, 1000 steps: speeds below 1e-10 m/s, mass to 1e-12.
2. Block falling onto a 10 deg STL floor: a 0.1 x 0.1 m block of T0 falls
   about 0.07 m onto a 3-cell-thick floor slab turned 10 deg against the mesh
   (wall friction 20 deg), open mesh below the slab. Material below the slab
   (through it) <= 1e-6 of the total; mass to 1e-10; the block stops (it
   slumps into a heap, as a cohesionless block must): at 0.6 s the kinetic
   energy <= 1e-3 of its peak and the centroid moved <= 1 mm along the
   slope in the last 0.1 s.
3. Tank with a void top: the WP 4.3a tank turned 10 deg against the mesh
   (tests/stlwalls/stl_wall_checks.py, case 3) without its lid, filled to
   the tank top with void above (freeSurface yes), damped to 0.5 s: mean
   stress along the tank axis within 5 % (the WP 4.3a tilted tolerance) of
   the same tank with mesh walls and freeSurface at every depth >= W; mass
   to 1e-12. (The brief compared with the lid version, without void; the
   free-surface tank itself differs from that by up to 17 % at the bottom,
   with mesh walls as with STL walls: WP 4.3c log, decision 0051.)

Prints PASS/FAIL per check, exits non-zero on any failure.
"""

import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "silo"))
sys.path.insert(0, os.path.join(HERE, "..", "void"))
sys.path.insert(0, os.path.join(HERE, "..", "lagrange"))
import cases as cs  # noqa: E402
import silo_checks as sc  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_wall_checks as sw  # noqa: E402

CELL = 0.01
TILT = 10.0
TANK_T = 0.5        # tank cases: damped settling to this time
TANK_TOL = 0.05     # as the WP 4.3a tilted tank


def add_bodies(path, bodies):
    """bodies: name -> (lo, hi, rot, friction): STL boxes in
    constant/rigidBodies."""
    body = ""
    for name, (lo, hi, rot, fr) in bodies.items():
        cs.write(f"{path}/constant/triSurface/{name}.stl", sw.box_stl(lo, hi, rot))
        body += f"    {name}\n    {{\n        surface {name}.stl;\n        friction {fr!r};\n    }}\n"
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def history(path):
    return np.loadtxt(f"{path}/history.dat")


# history columns of celVoidCheck
STEP, T, MASS, PX, PY, CX, CY, SMAX, SEMPTY, UMAX = range(10)
KE = 12
FXB, FYB, MBELOW = 15, 16, 17


def check_lake(work):
    print("[1] Lake at rest in an STL tank: side walls cut cells at 1/4 and 3/4, mesh")
    print("    floor, filled to 0.455 m (top cell half full), 1000 steps")
    # mesh x in [-0.0225, 0.1275] (15 cells): left wall face at x = 0 cuts
    # the cell [-0.0025, 0.0075] at a quarter; right wall face at x = 0.1
    # cuts [0.0975, 0.1075] at a quarter (open 0.25)
    p = vc.make_case(f"{work}/lake", 15, 50, vc.MC,
                     {"fill": [-1.0, -1.0, 1.0, 0.455], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "nSteps": 1000})
    move_mesh(p, -0.0225, 0.0)
    add_bodies(p, {
        "leftWall": ((-0.2, -0.1, -0.01), (0.0, 1.0, 0.02), 0.0, 0.0),
        "rightWall": ((0.1, -0.1, -0.01), (0.3, 1.0, 0.02), 0.0, 0.0)})
    t = vc.run(p)
    h = history(p)
    vmax = h[1:, UMAX].max()
    dm = abs(h[-1, MASS] / h[0, MASS] - 1)
    ok = vmax <= 1e-10 and dm <= 1e-12
    return sc.report("Lake at rest in an STL tank", ok, [
        f"largest speed {vmax:.1e} m/s (<= 1e-10), mass change {dm:.0e} (<= 1e-12), "
        f"{int(h[-1, STEP])} steps, {t:.0f} s"])


def move_mesh(path, x0, y0):
    """Translate the box mesh of a voidCheck case by (x0, y0)."""
    fn = f"{path}/system/blockMeshDict"
    lines = open(fn).read()
    i = lines.index("vertices")
    j = lines.index(");", i)
    head, verts, tail = lines[:i], lines[i:j], lines[j:]
    def mv(m):
        x, y, z = (float(v) for v in m.group(1).split())
        return f"({x + x0!r} {y + y0!r} {z!r})"
    verts = re.sub(r"\(([-0-9.e]+ [-0-9.e]+ [-0-9.e]+)\)", mv, verts)
    open(fn, "w").write(head + verts + tail)


def check_block(work):
    print("\n[2] T0 block falling onto a floor slab turned 10 deg against the mesh")
    th = math.radians(TILT)
    thick = 0.03
    # floor slab: top face through the origin at 10 deg, 3 cells thick
    p = vc.make_case(f"{work}/block", 40, 40, vc.MC,
                     {"fill": [0.15, 0.13, 0.25, 0.23], "U0": [0, 0, 0],
                      "tEnd": 0.6, "maxCo": 0.4, "writeTimes": [0.5]})
    move_mesh(p, 0.0, -0.1)
    add_bodies(p, {"floor": ((-0.5, -thick, -0.01), (1.0, 0.0, 0.02), TILT, 20.0)})
    t = vc.run(p)
    h = history(p)
    mass = h[0, MASS]
    dm = abs(h[-1, MASS] / mass - 1)
    leak = h[:, MBELOW].max() / mass
    c1 = np.loadtxt(f"{p}/cells.dat")
    c0 = np.loadtxt(f"{p}/cells_0.5.dat")

    def centroid(c):
        m = c[:, 3] * c[:, 10]
        return np.array([(m * c[:, 0]).sum(), (m * c[:, 1]).sum()]) / m.sum()
    # below the slab: signed distance below its bottom face, any material
    d = -math.sin(th) * c1[:, 0] + math.cos(th) * c1[:, 1]
    below = (c1[:, 3] * c1[:, 10] * (d < -thick - 0.5 * CELL)).sum() * CELL ** 2 / mass
    leak = max(leak, below)
    move = abs(np.dot(centroid(c1) - centroid(c0), [math.cos(th), math.sin(th)]))
    ke = h[-1, KE] / h[:, KE].max()
    ok = leak <= 1e-6 and dm <= 1e-10 and ke <= 1e-3 and move <= 1e-3
    return sc.report("Block on a 10 deg STL floor", ok, [
        f"material below the slab {leak:.1e} of the total (<= 1e-6), mass change "
        f"{dm:.0e} (<= 1e-10)",
        f"at {h[-1, T]:.2f} s: kinetic energy {ke:.1e} of its peak (<= 1e-3), centroid moved "
        f"{1e3 * move:.2f} mm along the slope since 0.5 s (<= 1 mm)",
        f"force of the floor on the material ({h[-1, FXB]:.1f}, {h[-1, FYB]:.1f}) N, "
        f"weight {mass * vc.G_ACC:.1f} N; {int(h[-1, STEP])} steps, {t:.0f} s"])


def main():
    work = sys.argv[1]
    which = sys.argv[2:] or ["lake", "block", "tank"]
    t0 = time.perf_counter()
    ok = True
    if "lake" in which:
        ok &= check_lake(work)
    if "block" in which:
        ok &= check_block(work)
    if "tank" in which:
        ok &= check_tank(work)
    print(f"\nWP 4.3c checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


def check_tank(work):
    print("\n[3] The 10 deg tank of WP 4.3a without its lid, void top, against the same")
    print("    tank with mesh walls and a void top")
    gm = sc.G_ACC
    W, H, T = sw.W, sw.H, 0.3
    th = math.radians(TILT)
    c, s_ = math.cos(th), math.sin(th)
    settle = dict(sw.SETTLE, tMax=TANK_T)
    # reference: mesh walls, freeSurface, filled to the top
    ref = sc.make_case(f"{work}/ref", W, H, True, settle, sc.T_RAMP, sc.ALPHA)
    open(f"{ref}/constant/materialProperties", "a").write("freeSurface yes;\n")
    open(f"{ref}/system/siloCheckDict", "a").write(f"fillBelow (0 1 0 {H!r});\n")
    sec0 = sc.run(ref)
    corners = np.array([[0, 0], [W, 0], [W, H], [0, H]]) @ np.array([[c, s_], [-s_, c]])
    lo, hi = corners.min(axis=0) - 2.3 * CELL, corners.max(axis=0) + 2.3 * CELL
    nx, ny = int(math.ceil((hi[0] - lo[0]) / CELL)), int(math.ceil((hi[1] - lo[1]) / CELL))
    z0 = (-0.01, 0.02)
    walls = {"leftWall": ((-T, -T, z0[0]), (0.0, H + T, z0[1]), TILT, sw.DELTA_W),
             "rightWall": ((W, -T, z0[0]), (W + T, H + T, z0[1]), TILT, sw.DELTA_W),
             "bottomWall": ((0.0, -T, z0[0]), (W, 0.0, z0[1]), TILT, sw.DELTA_W)}
    sw.SETTLE, keep = settle, sw.SETTLE
    top = sw.make_stl_case(f"{work}/void", float(lo[0]), nx, float(lo[1]), ny + 3,
                           walls, [], (gm * s_, -gm * c))
    sw.SETTLE = keep
    open(f"{top}/constant/materialProperties", "a").write("freeSurface yes;\n")
    open(f"{top}/system/siloCheckDict", "a").write(f"fillBelow ({-s_!r} {c!r} 0 {H!r});\n")
    sec1 = sc.run(top)
    zr, svr, _ = sw.profile(ref, 0.0)
    z, sv, _ = sw.profile(top, TILT)
    s = sc.summary(top)
    deep = z >= W - 1e-9
    e = np.abs(sv[deep] / svr[deep] - 1.0)
    det = [f"damped {int(s['steps'])} steps to {s['t']:.2f} s, then {int(s['contSteps'])} "
           f"undamped: stress change {100 * s['contMaxStressChange'] / (sc.RHO * gm * H):.2f} "
           f"% of rho g H, KE / strain energy {s['contKEratio']:.1e}",
           f"mean stress along the axis against the mesh-wall tank, depth >= W: "
           f"largest difference {100 * e.max():.2f} % (<= {100 * TANK_TOL:g} %)"]
    for frac in (0.5, 0.75, 1.0):
        i = int(np.argmin(np.abs(z - frac * H)))
        det.append(f"  depth {z[i]:.3f} m: {sv[i]:8.1f} Pa, mesh walls {svr[i]:8.1f} Pa")
    det.append(f"mass change {s['massChange']:.1e} (<= 1e-12); wall clock mesh walls "
               f"{sec0:.0f} s, STL {sec1:.0f} s")
    ok = e.max() <= TANK_TOL and abs(s["massChange"]) <= 1e-12
    return sc.report("Tank turned 10 deg with a void top", ok, det)


if __name__ == "__main__":
    main()
