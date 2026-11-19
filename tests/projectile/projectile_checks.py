"""WP 5.4 done-when check: a projectile into clay (target problem T3), the
projectile free (no drive) and coupled both ways. Run by
tests/test_projectile.sh after celVoidCheck is built:

    python3 projectile_checks.py WORKDIR

A flat-nosed slug 0.095 m wide and 0.1 m tall (plane strain) starts 7.5 mm
above a 0.15 m layer of Tresca clay (cu 10 kPa, E 3 MPa, nu 0.45 (undrained), rho 1800,
gravity) moving down at 1 m/s, smooth (friction 0). Its mass is that of a
tungsten slug (rho 19300) in the 10 mm slab of the mesh: 1.83 kg. Twice:

- rigid: a free rigid body (decision 0058) of that mass;
- deformable: a steel hex8 body (decision 0072, 0073; E 200 GPa, 10 x 10 x
  1 elements, 0.03 m deep, its density set so its mass is the same).

Reference: the energy of the slug at impact plus the work of gravity over
the penetration d equals the work of the smooth strip footing's bearing
resistance, Prandtl's (2 + pi) cu B with Skempton's depth factor
(1 + 0.2 d / B):

    0.5 m v_i^2 + m g d = (2 + pi) cu B s (d + 0.1 d^2 / B),

v_i the speed at impact, s the slab depth (0.01 m). Inertia of the clay and
rate effects are left out. Checks (decision 0074):

- deformable within 5 % of rigid (the reference run: the rigid machinery
  is checked against Prandtl in WP 4.4);
- final penetration of each within 30 % of the energy estimate (a sanity
  band: on 2026-11-18 both read 26 % short, decision 0074);
- both at rest (speed below 2 % of the impact speed over the last 0.05 s);
- mass kept to 1e-10, nothing lost in hand-overs, density within 3 % (the
  clay's elastic compression under the impact, bulk modulus 10 MPa).

Prints PASS/FAIL, exits non-zero on failure.
"""

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("silo", "void", "lagrange", "stlwalls", "freebody"):
    sys.path.insert(0, os.path.join(HERE, "..", d))
import cases as cs  # noqa: E402
import silo_checks as sc  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_push_checks as sp  # noqa: E402
import free_body_checks as fb  # noqa: E402

CELL = 0.01
DEPTH = 0.15                    # clay layer
X0, X1 = 0.1025, 0.1975         # slug, B = 0.095 m
Y0, HEIGHT = 0.1575, 0.1
B = X1 - X0
V0 = 1.0                        # m/s, down
RHO_SLUG = 19300.0              # tungsten
MASS = RHO_SLUG * B * HEIGHT * CELL
DEPTH_B = 0.03                  # deformable body depth in z (-0.01 to 0.02)
CU, E_CLAY, NU_CLAY = 10e3, 3e6, 0.45
T_END = 0.3
CX, CY, VX, VY = 23, 24, 25, 26
TOL, TOL_RIGID, RHO_TOL = 0.30, 0.05, 1.03


def reference():
    vi2 = V0 ** 2 + 2 * vc.G_ACC * (Y0 - DEPTH)
    R = (2 + math.pi) * CU * B * CELL
    # 0.1 R / B d^2 + (R - m g) d - 0.5 m vi^2 = 0
    a, b, c = 0.1 * R / B, R - MASS * vc.G_ACC, -0.5 * MASS * vi2
    return (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a), math.sqrt(vi2)


def add_deformable(path):
    rho = MASS / (B * HEIGHT * DEPTH_B)
    body = (f"    slug\n    {{\n        deformable\n        {{\n"
            f"            min ({X0!r} {Y0!r} -0.01);\n"
            f"            max ({X1!r} {Y0 + HEIGHT!r} 0.02);\n"
            f"            cells (10 10 1);\n            density {rho!r};\n"
            f"            material {{ type linearElastic; E 2e11; nu 0.3; }}\n"
            f"            planeStrain yes;\n            damping 0.2;\n        }}\n"
            f"        velocity constant (0 {-V0!r} 0);\n        friction 0;\n    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def run_case(work, kind):
    mat = (f"type mohrCoulomb;\nE {E_CLAY!r};\nnu {NU_CLAY!r};\nphi 0;\nc {CU!r};\npsi 0;\n"
           "fullBlend yes;\n")
    p = vc.make_case(f"{work}/{kind}", 30, 30, mat,
                     {"fill": [-1.0, -1.0, 1.0, DEPTH], "U0": [0, 0, 0],
                      "initK0": NU_CLAY / (1 - NU_CLAY), "tEnd": T_END, "maxCo": 0.4})
    if kind == "rigid":
        fb.add_free_body(p, "slug", (X0, Y0, -0.01), (X1, Y0 + HEIGHT, 0.02), MASS, 0.0,
                         (0.0, -V0), 0.0, 0.0)
    else:
        add_deformable(p)
    sec = vc.run(p)
    return np.loadtxt(f"{p}/history.dat"), sec


def summary(name, h, s, ref, vImp):
    mass = h[0, sp.MASS]
    dm = abs(h[-1, sp.MASS] / mass - 1)
    lost = h[:, sp.LOST].max() / mass
    rho = h[:, sp.RHOMAX].max()
    d = -(h[-1, CY] - h[0, CY]) - (Y0 - DEPTH)
    late = h[:, sp.T] >= T_END - 0.05
    speed = np.hypot(h[late, VX], h[late, VY]).max() / vImp
    err = d / ref - 1
    ok = dm <= 1e-10 and lost == 0 and rho <= RHO_TOL and speed <= 0.02 and abs(err) <= TOL
    line = (f"{name}: penetration {1000 * d:.2f} mm ({100 * err:+.1f} % of the reference, "
            f"{100 * TOL:.0f} %), last 0.05 s speed {100 * speed:.2f} % of impact (2 %), "
            f"mass change {dm:.0e}, lost {lost:.0e}, density {rho:.4f} rho0 ({RHO_TOL}), {s:.0f} s")
    return ok, d, line


def check(work):
    ref, vImp = reference()
    print(f"[1] Slug ({MASS:.3f} kg, B {B} m) into Tresca clay at {V0} m/s: rigid and steel hex8")
    hr, sr = run_case(work, "rigid")
    hd, sd = run_case(work, "steel")
    okr, dr, lr = summary("rigid", hr, sr, ref, vImp)
    okd, dd, ld = summary("steel", hd, sd, ref, vImp)
    e = dd / dr - 1
    ok = okr and okd and abs(e) <= TOL_RIGID
    for name, h in (("rigid", hr), ("steel", hd)):
        np.savetxt(f"{work}/{name}_penetration.dat",
                   np.c_[h[:, sp.T], -(h[:, CY] - h[0, CY]) - (Y0 - DEPTH), h[:, VY]],
                   header="t penetration vy")
    return sc.report("Projectile into clay", ok, [
        f"reference (Prandtl with depth factor, energy): {1000 * ref:.2f} mm, impact speed "
        f"{vImp:.3f} m/s", lr, ld,
        f"steel against rigid: {100 * e:+.2f} % ({100 * TOL_RIGID:.0f} %)"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    ok = check(work)
    print(f"\nWP 5.4 checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
