"""WP 4.3d done-when check: free rigid bodies (decision 0058), whose motion
follows from the force of the material on them and gravity. Run by
tests/test_free_body.sh after celVoidCheck is built:

    python3 free_body_checks.py WORKDIR

[1] Free flight: a block thrown in the void above a thin layer, v0 (0.5, 1)
    m/s, omega 5 rad/s about z. Its centre follows the parabola and its
    angular velocity stays constant (no force, no moment) to 1e-9.
[2] Drop: a steel block 0.095 m wide and 0.1 m tall (rho 7850, so 0.746 kg
    in the 10 mm slab the mesh holds; bearing pressure 7.7 kPa) falls from
    7.5 mm onto a 0.15 m T0 layer (as WP 4.3b, decision 0053), wall
    friction 20 deg, and comes to rest. No inertia given, so it does not
    rotate (with rotation it tilts on impact, decision 0058). Checks:
    - mass kept to 1e-10, nothing lost in hand-overs (no leakage);
    - momentum of material and block together against gravity and the
      mesh walls: largest difference <= 1 % of the block's weight impulse;
    - at rest: over the last 0.05 s the block's speed below 1 % of its
      impact speed, and the mean force of the material on it within 2 %
      of its weight;
    - largest material density within 2 % of rho0.

Prints PASS/FAIL, exits non-zero on failure.
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("silo", "void", "lagrange", "stlwalls"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import silo_checks as sc  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_wall_checks as sw  # noqa: E402
import stl_push_checks as sp  # noqa: E402

CELL = 0.01
DEPTH = 0.15
X0, X1 = 0.1025, 0.1975
Y0, HEIGHT = 0.1575, 0.1
RHO_B = 7850.0
CX, CY, VX, VY, WZ = 23, 24, 25, 26, 27


def add_free_body(path, name, lo, hi, mass, inertia, v0, w0, fr):
    cs.write(f"{path}/constant/triSurface/{name}.stl", sw.box_stl(lo, hi))
    J = inertia
    inert = f"        inertia ({J!r} 0 0 {J!r} 0 {J!r});\n" if J > 0 else ""
    body = (f"    {name}\n    {{\n        surface {name}.stl;\n        friction {fr!r};\n"
            f"        mass {mass!r};\n{inert}"
            f"        velocity ({v0[0]!r} {v0[1]!r} 0);\n        omega (0 0 {w0!r});\n    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def block_mass():
    w, h = X1 - X0, HEIGHT
    m = RHO_B * w * h * CELL
    return m, m * (w * w + h * h) / 12


def check_flight(work):
    print("[1] Free flight: block thrown in the void, v0 (0.5, 1) m/s, omega 5 rad/s")
    tEnd = 0.1
    p = vc.make_case(f"{work}/flight", 30, 30, vc.MC + sp.T0_RULES,
                     {"fill": [-1.0, -1.0, 1.0, 0.03], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "tEnd": tEnd, "maxCo": 0.4})
    m, J = block_mass()
    v0, w0 = (0.5, 1.0), 5.0
    lo, hi = (0.08, 0.1, -0.01), (0.12, 0.14, 0.02)
    add_free_body(p, "block", lo, hi, m, J, v0, w0, 0.0)
    t = vc.run(p)
    h = np.loadtxt(f"{p}/history.dat")
    c0 = h[0, [CX, CY]]     # STL points are single precision
    tt = h[:, sp.T]
    # the scheme (central difference, staggered velocity): each step moves
    # the centre with the velocity before it, then adds g dt
    dt = np.diff(np.r_[0.0, tt])
    vy = v0[1] - vc.G_ACC * np.cumsum(dt)
    cy = c0[1] + np.cumsum(dt * np.r_[v0[1], vy[:-1]])
    ex = np.abs(h[:, CX] - (c0[0] + v0[0] * tt)).max()
    ey = max(np.abs(h[:, CY] - cy).max(), np.abs(h[:, VY] - vy).max())
    par = np.abs(h[:, CY] - (c0[1] + v0[1] * tt - 0.5 * vc.G_ACC * tt ** 2)).max()
    ew = np.abs(h[:, WZ] - w0).max()
    ok = max(ex, ey) <= 1e-9 and ew <= 1e-9
    return sc.report("Free flight", ok, [
        f"centre and velocity against the central-difference recurrence: "
        f"{max(ex, ey):.1e} (<= 1e-9); against the exact parabola {par:.1e} m "
        f"(first order in dt)",
        f"omega change {ew:.1e} rad/s (<= 1e-9)",
        f"{int(h[-1, sp.STEP])} steps, {t:.0f} s"])


def check_drop(work):
    print("[2] Steel block dropped from 7.5 mm onto a T0 layer, comes to rest")
    tEnd = 0.4
    p = vc.make_case(f"{work}/drop", 30, 30, vc.MC + sp.T0_RULES,
                     {"fill": [-1.0, -1.0, 1.0, DEPTH], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "tEnd": tEnd, "maxCo": 0.4})
    m, J = block_mass()
    # no inertia: the block does not rotate (decision 0058: with rotation
    # the impact kicks it to 0.5 rad/s and its corner squeezes a cut cell)
    add_free_body(p, "block", (X0, Y0, -0.01), (X1, Y0 + HEIGHT, 0.02), m, 0.0,
                  (0.0, 0.0), 0.0, 20.0)
    t = vc.run(p)
    h = np.loadtxt(f"{p}/history.dat")
    mass = h[0, sp.MASS]
    dm = abs(h[-1, sp.MASS] / mass - 1)
    lost = h[:, sp.LOST].max() / mass
    rho = h[:, sp.RHOMAX].max()
    W = m * vc.G_ACC

    # momentum of material and block: d(P + m v) = sum dt (F_walls + (M + m) g)
    dt = np.diff(h[:, sp.T])
    pyb = m * h[:, VY]
    pxb = m * h[:, VX]
    fx = h[1:, sp.FXW]
    fy = h[1:, sp.FYW] - (h[1:, sp.MASS] + m) * vc.G_ACC
    ex = (h[1:, sp.PX] + pxb[1:] - h[0, sp.PX] - pxb[0]) - np.cumsum(dt * fx)
    ey = (h[1:, sp.PY] + pyb[1:] - h[0, sp.PY] - pyb[0]) - np.cumsum(dt * fy)
    bal = max(np.abs(ex).max(), np.abs(ey).max()) / (W * tEnd)

    vImp = math.sqrt(2 * vc.G_ACC * (Y0 - DEPTH))
    late = h[:, sp.T] >= tEnd - 0.05
    speed = np.hypot(h[late, VX], h[late, VY]).max() / vImp
    fb = -h[late, sp.FYB].mean() / W     # force of the material on the block
    sink = DEPTH - (h[-1, CY] - HEIGHT / 2)
    ok = (dm <= 1e-10 and lost == 0 and bal <= 0.01 and speed <= 0.01
          and abs(fb - 1) <= 0.02 and rho <= 1.02)
    return sc.report("Block drop", ok, [
        f"mass change {dm:.0e} (<= 1e-10), lost in hand-overs {lost:.0e} (0)",
        f"momentum of material and block: largest difference {100 * bal:.3f} % of the "
        f"weight impulse (<= 1 %)",
        f"last 0.05 s: speed {100 * speed:.2f} % of the impact speed {vImp:.3f} m/s "
        f"(<= 1 %), force / weight {fb:.4f} (1 +- 0.02)",
        f"sinks {1000 * sink:.1f} mm into the layer; omega at the end {h[-1, WZ]:.1e} rad/s",
        f"largest material density {rho:.4f} rho0 (<= 1.02)",
        f"{int(h[-1, sp.STEP])} steps, {t:.0f} s"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    ok = check_flight(work)
    ok = check_drop(work) and ok
    print(f"\nWP 4.3d checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
