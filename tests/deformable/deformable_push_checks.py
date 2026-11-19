"""WP 5.3a done-when check: two-way coupling of a deformable body (hex8,
WP 5.1) with the material (decision 0072). Run by
tests/test_deformable_push.sh after celVoidCheck is built:

    python3 deformable_push_checks.py WORKDIR [--soft]

The footing push of WP 4.3b (decision 0053: 0.095 m wide, 0.1 m tall,
pushed 4.25 cm into a T0 layer at 0.25 m/s and pulled out) twice: rigid
(STL) and as a steel body of 10 x 10 x 1 hex8 elements (one element deep,
plane strain, damping 0.2) driven at its top nodes. Checks:

- both runs: mass kept to 1e-10, nothing lost in hand-overs, material
  density within 2 % of rho0, finite forces;
- the steel footing's force follows the rigid one: mean over the push
  (0.06 to 0.19 s) within 5 % (single values swing with the cut cells for
  both, from about 20 to 64 N);
- after pulling out, the force below 1 % of its peak (both).

--soft (WP 5.3b, tests/test_soft_push.sh, decision 0073): a soft body
(E 20 MPa, nu 0.3, rho 1100, hourglass 0.3), drive ramped over 0.02 s:
shortening against the elastic column in plane strain,
dh = F h (1 - nu^2) / (E A), A its width times its depth (0.03 m), within
20 % over 0.15 to 0.2 s, and the material density within 2 %. On
2026-11-18: ratio 0.90, density 1.016 rho0 (before decision 0073: 1.44 and
1.62 rho0).

Prints PASS/FAIL, exits non-zero on failure.
"""

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
import stl_push_checks as sp  # noqa: E402

DEPTH_B = 0.03           # body depth in z (-0.01 to 0.02)
HGT = 27                 # history column: body height (deformable)
WINDOW = (0.06, 0.19)


def table(tEnd, ramp=0.0):
    v = sp.SPEED
    if ramp > 0:
        return (f"table ((0 (0 0 0)) ({ramp!r} (0 {-v!r} 0)) ({sp.T_PUSH!r} (0 {-v!r} 0)) "
                f"({sp.T_PUSH + ramp!r} (0 {v!r} 0)) ({tEnd + 1!r} (0 {v!r} 0)))")
    return (f"table ((0 (0 {-v!r} 0)) ({sp.T_PUSH!r} (0 {-v!r} 0)) "
            f"({sp.T_PUSH + 1e-4!r} (0 {v!r} 0)) ({tEnd + 1!r} (0 {v!r} 0)))")


def add_deformable(path, E, nu, rho, tbl, hourglass=None):
    lo = (sp.X0, sp.Y0, -0.01)
    hi = (sp.X1, sp.Y0 + sp.HEIGHT, 0.02)
    top = sp.Y0 + sp.HEIGHT
    body = (f"    footing\n    {{\n        deformable\n        {{\n"
            f"            min ({lo[0]!r} {lo[1]!r} {lo[2]!r});\n"
            f"            max ({hi[0]!r} {hi[1]!r} {hi[2]!r});\n"
            f"            cells (10 10 1);\n            density {rho!r};\n"
            f"            material {{ type linearElastic; E {E!r}; nu {nu!r}; }}\n"
            f"            planeStrain yes;\n            damping 0.2;\n"
            + (f"            hourglass {hourglass!r};\n" if hourglass else "") +
            f"            drive ({sp.X0 - 0.01!r} {top - 1e-6!r} -1) "
            f"({sp.X1 + 0.01!r} {top + 1e-6!r} 1);\n        }}\n"
            f"        velocity {tbl};\n        friction {sp.FRICTION!r};\n"
            f"    }}\n")
    cs.write(f"{path}/constant/rigidBodies",
             cs.HEADER.format(cls="dictionary", obj="rigidBodies")
             + "bodies\n{\n" + body + "}\n")


def run_case(work, name, kind, ramp=0.0, hourglass=None):
    tEnd = 2 * sp.T_PUSH
    p = vc.make_case(f"{work}/{name}", 30, 30, vc.MC + sp.T0_RULES,
                     {"fill": [-1.0, -1.0, 1.0, sp.DEPTH], "U0": [0, 0, 0],
                      "initK0": vc.NU / (1 - vc.NU), "tEnd": tEnd, "maxCo": 0.4})
    if kind == "rigid":
        sp.add_moving_body(p, "footing", (sp.X0, sp.Y0, -0.01),
                           (sp.X1, sp.Y0 + sp.HEIGHT, 0.02), table(tEnd), sp.FRICTION)
    else:
        E, nu, rho = kind
        add_deformable(p, E, nu, rho, table(tEnd, ramp), hourglass)
    sec = vc.run(p)
    return np.loadtxt(f"{p}/history.dat"), sec


def basic(name, h, s):
    mass = h[0, sp.MASS]
    dm = abs(h[-1, sp.MASS] / mass - 1)
    lost = h[:, sp.LOST].max() / mass
    rho = h[:, sp.RHOMAX].max()
    fin = bool(np.all(np.isfinite(h[:, sp.FYB])))
    good = dm <= 1e-10 and lost == 0 and rho <= 1.02 and fin
    return good, (f"{name}: mass change {dm:.0e}, lost {lost:.0e}, density {rho:.4f} rho0 "
                  f"(1.02), peak force {np.abs(h[:, sp.FYB]).max():.1f} N, {s:.0f} s")


def mean_force(h):
    w = (h[:, sp.T] >= WINDOW[0]) & (h[:, sp.T] <= WINDOW[1])
    return h[w, sp.FYB].mean()


def check_steel(work):
    print("[1] Footing push of WP 4.3b: rigid and steel (hex8, two-way)")
    tEnd = 2 * sp.T_PUSH
    hr, sr = run_case(work, "rigid", "rigid")
    hs, ss = run_case(work, "steel", (2e11, 0.3, 7800.0))
    ok = True
    det = []
    for name, h, s in (("rigid", hr, sr), ("steel", hs, ss)):
        good, line = basic(name, h, s)
        ok = ok and good
        det.append(line)
    e = mean_force(hs) / mean_force(hr) - 1
    ok = ok and abs(e) <= 0.05
    det.append(f"mean force over {WINDOW[0]} to {WINDOW[1]} s: steel {-mean_force(hs):.2f} N, "
               f"rigid {-mean_force(hr):.2f} N ({100 * e:+.2f} %, 5 %)")
    for name, h in (("rigid", hr), ("steel", hs)):
        late = h[:, sp.T] >= tEnd - 0.05
        tail = np.abs(h[late, sp.FYB]).max() / np.abs(h[:, sp.FYB]).max()
        ok = ok and tail <= 0.01
        det.append(f"{name} after pulling out: {100 * tail:.2f} % of the peak (1 %)")
    return sc.report("Steel footing (deformable, two-way) against the rigid one", ok, det)


def check_soft(work):
    print("[2] Soft footing (hex8, E 20 MPa), drive ramped over 0.02 s")
    E, nu = 2e7, 0.3
    h, s = run_case(work, "soft", (E, nu, 1100.0), ramp=0.02, hourglass=0.3)
    good, line = basic("soft", h, s)
    w = (h[:, sp.T] >= 0.15) & (h[:, sp.T] <= 0.2)
    F = -h[w, sp.FYB]
    dh = sp.HEIGHT - h[w, HGT]
    A = (sp.X1 - sp.X0) * DEPTH_B
    ref = F * sp.HEIGHT * (1 - nu ** 2) / (E * A)
    ratio = dh.mean() / ref.mean()
    ok = good and abs(ratio - 1) <= 0.2
    return sc.report("Soft footing shortening", ok, [
        line,
        f"shortening {1000 * dh.mean():.3f} mm against the elastic column "
        f"{1000 * ref.mean():.3f} mm (ratio {ratio:.3f}, 20 %)"])


def main():
    work = sys.argv[1]
    t0 = time.perf_counter()
    if "--soft" in sys.argv[2:]:
        ok = check_soft(work)
        tag = "5.3b"
    else:
        ok = check_steel(work)
        tag = "5.3a"
    print(f"\nWP {tag} checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
