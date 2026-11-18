"""WP 4.3b done-when check: a moving STL body pushes the material through its
cut cells (decision 0042) and reports the force on it; the moving cut cells
are decision 0053. Run by tests/test_stl_push.sh after celVoidCheck is built:

    python3 stl_push_checks.py WORKDIR

Footing push: a 0.3 m wide box (mesh walls, no wall friction) holds a 0.15 m
layer of T0 (cohesionless Mohr-Coulomb, rho 1800, E 10 MPa, nu 0.3, phi 35)
with the T0 rule set of WP 3.6e (mu(I), flowImpedance, normalCap, shearCap,
spreadLayers) at its K0 stress, void above, 10 mm cells. (Plain
Mohr-Coulomb without flowImpedance has no nearly-full snap: the footing's
corner then squeezes a nearly full surface cell to 1.39 rho0, decision 0053.) A rigid footing 0.095 m wide and
0.1 m tall (an STL box, faces off the mesh lines: sides cut cells at a quarter,
bottom at three quarters; wall friction 20 deg) starts 7.5 mm above the
surface, moves down at 0.25 m/s for 0.2 s (4.25 cm into the layer), then up
at 0.25 m/s for 0.2 s, back to its start. Checks:

- mass kept to 1e-10 and no mass lost in hand-overs from covered cells
  (no leakage: material cannot enter the body, so leakage would show as
  lost mass);
- while pushing (0.1 s and 0.2 s): no void under the footing (cells below its
  bottom within its width, and the cut cells there, alpha >= 0.99);
- material density rho/alpha within 2 % of rho0 in every cell at every step
  (no material squeezed into the swept volume);
- momentum of the material against the forces on it (bodies, mesh walls,
  gravity), summed over the run: largest difference <= 1 % of the impulse
  of the footing force;
- after pulling out: over the last 0.05 s the force on the footing below
  1 % of its peak (no tension, the footing separates).

Prints PASS/FAIL, exits non-zero on failure.
"""

import os
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
SPEED = 0.25        # m/s, down then up
T_PUSH = 0.2        # s
DEPTH = 0.15        # layer
X0, X1 = 0.1025, 0.1975
Y0, HEIGHT = 0.1575, 0.1
FRICTION = 20.0
T0_RULES = ("mu2 0.961813;\nI0 0.279;\nd 0.001;\nrhoS 2650.0;\nflowImpedance yes;\n"
            "normalCap yes;\nshearCap yes;\nspreadLayers yes;\n")

# history columns of celVoidCheck
STEP, T, MASS, PX, PY = range(5)
FXB, FYB = 15, 16
FXW, FYW, LOST, RHOMAX = 18, 19, 20, 21


def add_moving_body(path, name, lo, hi, velocity, fr):
    cs.write(f"{path}/constant/triSurface/{name}.stl", sw.box_stl(lo, hi))
    body = (f"    {name}\n    {{\n        surface {name}.stl;\n        friction {fr!r};\n"
            f"        velocity {velocity};\n    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def check_push(work):
    print("[1] Footing pushed 4.25 cm into a T0 layer at 0.25 m/s, then pulled out")
    tEnd = 2 * T_PUSH
    p = vc.make_case(f"{work}/push", 30, 30, vc.MC + T0_RULES,
                     {"fill": [-1.0, -1.0, 1.0, DEPTH], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "tEnd": tEnd, "maxCo": 0.4,
                      "writeTimes": [0.1, T_PUSH]})
    v = SPEED
    table = (f"table ((0 (0 {-v!r} 0)) ({T_PUSH!r} (0 {-v!r} 0)) "
             f"({T_PUSH + 1e-4!r} (0 {v!r} 0)) ({tEnd + 1!r} (0 {v!r} 0)))")
    add_moving_body(p, "footing", (X0, Y0, -0.01), (X1, Y0 + HEIGHT, 0.02), table, FRICTION)
    t = vc.run(p)
    h = np.loadtxt(f"{p}/history.dat")
    mass = h[0, MASS]
    dm = abs(h[-1, MASS] / mass - 1)
    lost = h[:, LOST].max() / mass
    rho = h[:, RHOMAX].max()

    # no void under the footing while pushing
    gap = 0.0
    for tw in (0.1, T_PUSH):
        c = np.loadtxt(f"{p}/cells_{tw:g}.dat")
        yb = Y0 - SPEED * tw
        under = ((c[:, 0] > X0 + 0.5 * CELL) & (c[:, 0] < X1 - 0.5 * CELL)
                 & (c[:, 1] < yb) & (c[:, 1] > yb - 2 * CELL) & (c[:, 10] > 0))
        gap = max(gap, 1 - c[under, 2].min())

    # momentum balance: dP = sum dt (F_body + F_walls + m g)
    dt = np.diff(h[:, T])
    fx = h[1:, FXB] + h[1:, FXW]
    fy = h[1:, FYB] + h[1:, FYW] - h[1:, MASS] * vc.G_ACC
    ex = (h[1:, PX] - h[0, PX]) - np.cumsum(dt * fx)
    ey = (h[1:, PY] - h[0, PY]) - np.cumsum(dt * fy)
    imp = np.sum(dt * np.abs(h[1:, FYB]))
    bal = max(np.abs(ex).max(), np.abs(ey).max()) / imp

    # separation after pulling out
    peak = np.abs(h[:, FYB]).max()
    late = h[:, T] >= tEnd - 0.05
    tail = np.abs(h[late, FYB]).max() / peak

    ok = (dm <= 1e-10 and lost == 0 and gap <= 0.01 and rho <= 1.02
          and bal <= 0.01 and tail <= 0.01)
    ipk = int(np.argmax(np.abs(h[:, FYB])))
    return sc.report("Footing push and pull-out", ok, [
        f"mass change {dm:.0e} (<= 1e-10), mass lost in hand-overs {lost:.0e} (0)",
        f"void under the footing at 0.1 and 0.2 s: {gap:.1e} (<= 0.01)",
        f"largest material density {rho:.4f} rho0 (<= 1.02)",
        f"momentum against the forces: largest difference {100 * bal:.3f} % of the "
        f"footing impulse (<= 1 %)",
        f"footing force: peak {peak:.0f} N at {h[ipk, T]:.3f} s "
        f"(weight of the layer {mass * vc.G_ACC:.0f} N), last 0.05 s {100 * tail:.2f} % "
        f"of the peak (<= 1 %)",
        f"{int(h[-1, STEP])} steps, {t:.0f} s"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    ok = check_push(work)
    print(f"\nWP 4.3b checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
