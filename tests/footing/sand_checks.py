"""WP 4.4c done-when check: smooth strip footing on sand, bearing capacity
factor N_gamma (q = 0.5 gamma B N_gamma) against the exact value of
Martin (2005, "Exact bearing capacity calculations using the method of
characteristics", Proc. IACMAG, Turin), as tabulated in the Optum G2/GX
verification "Bearing capacity" (strip footing on sand): phi 30 deg, smooth
footing, N_gamma = 7.653 (rough 14.754). Run by tests/manual/test_sand_footing.sh
after celVoidCheck is built:

    python3 sand_checks.py WORKDIR [--fine] [--psi DEG] [--friction DEG]
        [--noblend] [--maxco CO] [--x0 X] [--gap G]

--fine runs 5 mm cells; --psi sets the dilation angle (default 30, associated);
--friction the footing's base friction (0, smooth); --noblend drops fullBlend;
--maxco, --x0 (left edge) and --gap (base above the surface) as named
(decision 0091); --width, --depth (layer, m), --speed (m/s), --E (Pa) and
--smax (settlement / B) as named, --novoidflow drops voidPlasticFlow, and
--stable checks the A.4 done-when instead of N_gamma (decision 0095).
The reading depends on the push speed times the wave impedance rho c
(numerical viscosity, decision 0095): slower pushes or a softer E read lower.
Not yet passing (decision 0057): see tests/manual/test_sand_footing.sh.

Sand: cohesionless Mohr-Coulomb, phi 30, psi 30 (associated, as the exact
solution assumes), rho 1800 (gamma 17.7 kN/m3), E 10 MPa (loose sand),
nu 0.3, at its K0 stress, 0.6 m wide and 0.3 m deep, frictionless mesh
walls, void above, 10 mm cells, fullBlend yes (decision 0056). Smooth rigid
footing B = 0.1 m starting 2.5 mm above the surface, speed ramped to
0.01 m/s over 0.05 s. The footing sinks as it loads, adding a surcharge
gamma s on its sides; the reference at settlement s is
q_ref = 0.5 gamma B N_gamma + gamma s N_q (N_q = 18.40, Prandtl-Reissner).
Check: the peak of q / q_ref over the run within 15 % of 1; mass kept.

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
import silo_checks as sc  # noqa: E402
import void_checks as vc  # noqa: E402
import stl_push_checks as sp  # noqa: E402

CELL = 0.01
B = 0.1
PHI = 30.0
E_SAND = 10e6
NU = 0.3
RHO = 1800.0
NX, NY = 60, 34
DEPTH = 0.3
X0 = 0.25 + 0.0025
GAP = 0.0025
SPEED, T_RAMP = 0.01, 0.05
S_MAX = 0.06 * B
NGAMMA = 7.653          # Martin (2005), phi 30, smooth
NQ = math.exp(math.pi * math.tan(math.radians(PHI))) * math.tan(math.radians(45 + PHI / 2)) ** 2
TOL = 0.15
PSI = PHI
NOBLEND = False
MAXCO = 0.4
FRICTION = 0.0         # footing base friction angle [deg]; 0 = smooth
VOIDFLOW = True        # voidPlasticFlow (decision 0095); --novoidflow drops it
STABLE = False         # --stable: A.4 done-when (runs through, bounded), not N_gamma


def travel(t):
    t = np.asarray(t)
    return np.where(t < T_RAMP, 0.5 * SPEED * t ** 2 / T_RAMP,
                    0.5 * SPEED * T_RAMP + SPEED * (t - T_RAMP))


def check_sand(work):
    print(f"[1] Smooth strip footing B = 0.1 m ({B / CELL:.0f} cells) on sand, phi 30, psi {PSI:g},")
    print("    gamma 17.7 kN/m3, E 10 MPa, pushed to a settlement of 0.06 B")
    tEnd = T_RAMP + (GAP + S_MAX) / SPEED - 0.5 * T_RAMP
    mat = (f"type mohrCoulomb;\nE {E_SAND!r};\nnu {NU!r};\nphi {PHI!r};\nc 0;\n"
           f"psi {PSI!r};\n" + ("" if NOBLEND else "fullBlend yes;\n")
           + ("voidPlasticFlow yes;\n" if VOIDFLOW else ""))
    k0 = 1 - math.sin(math.radians(PHI))
    p = vc.make_case(f"{work}/sand", NX, NY, mat,
                     {"fill": [-1.0, -1.0, 10.0, DEPTH], "U0": [0, 0, 0],
                      "initK0": k0, "tEnd": tEnd, "maxCo": MAXCO})
    v = SPEED
    table = f"table ((0 (0 0 0)) ({T_RAMP!r} (0 {-v!r} 0)) ({tEnd + 1!r} (0 {-v!r} 0)))"
    y0 = DEPTH + GAP
    sp.add_moving_body(p, "footing", (X0, y0, -0.01), (X0 + B, y0 + 0.1, 0.02), table, FRICTION)
    sec = vc.run(p)
    h = np.loadtxt(f"{p}/history.dat")
    gamma = RHO * vc.G_ACC
    s = travel(h[:, sp.T]) - GAP
    q = -h[:, sp.FYB] / (B * CELL)
    qref = 0.5 * gamma * B * NGAMMA + gamma * np.maximum(s, 0) * NQ
    r = q / qref
    on = s > 0
    ipk = int(np.argmax(np.where(on, r, -1)))
    mass = h[0, sp.MASS]
    dm = abs(h[-1, sp.MASS] / mass - 1)
    lost = h[:, sp.LOST].max() / mass
    err = r[ipk] - 1
    rhoMax = h[:, sp.RHOMAX].max()
    if STABLE:
        # A.4 (decision 0095): the associated push runs to 0.06 B with the
        # load bounded (no pressure build-up in the surface cells), mass
        # kept, no packing. N_gamma itself is the manual check.
        done = h[-1, sp.T] >= tEnd - 1e-6 and s[-1] >= S_MAX - 1e-6
        late = r[s > 0.5 * S_MAX]
        ok = (done and dm <= 1e-10 and lost == 0 and rhoMax <= 1.01
              and np.all(np.isfinite(q)) and late.max() < 3.0
              and late.max() / late.min() < 1.3)
    else:
        ok = abs(err) <= TOL and dm <= 1e-10 and lost == 0
    det = [f"peak q / q_ref {r[ipk]:.3f} at s = {s[ipk] / B:.3f} B (limit +-{100 * TOL:.0f} %); "
           f"q {q[ipk]:.0f} Pa, 0.5 gamma B N_gamma = {0.5 * gamma * B * NGAMMA:.0f} Pa, "
           f"N_q {NQ:.2f}"]
    for sb in (0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15):
        if sb > S_MAX / B + 1e-9:
            break
        i = int(np.argmin(np.abs(s / B - sb)))
        det.append(f"  s = {s[i] / B:.3f} B: q {q[i]:7.0f} Pa, q / q_ref {r[i]:.3f}, "
                   f"implied N_gamma {(q[i] - gamma * max(s[i], 0) * NQ) / (0.5 * gamma * B):.2f}")
    if STABLE:
        det.append(f"stable push (A.4): reached s = {s[-1] / B:.3f} B of {S_MAX / B:.2f} B; "
                   f"q / q_ref over the second half {late.min():.2f} to {late.max():.2f} "
                   f"(below 3, spread below 1.3); density {rhoMax:.4f} rho0 (<= 1.01)")
    det.append(f"mass change {dm:.0e} (<= 1e-10), lost in hand-overs {lost:.0e}; "
               f"{int(h[-1, sp.STEP])} steps, {sec:.0f} s")
    np.savetxt(f"{p}/q.dat", np.c_[s / B, q, r], header="s/B q[Pa] q/qref")
    name = ("Associated sand footing runs through (A.4)" if STABLE
            else "Strip footing on sand (N_gamma)")
    return sc.report(name, ok, det)


def main():
    global VOIDFLOW, STABLE, CELL, NX, NY, PSI, FRICTION, NOBLEND, MAXCO, X0, GAP, DEPTH, SPEED, E_SAND, S_MAX
    work = sys.argv[1]
    if "--fine" in sys.argv[2:]:
        vc.CELL = CELL = 0.005
        NX, NY = 120, 68
    if "--psi" in sys.argv[2:]:
        PSI = float(sys.argv[sys.argv.index("--psi") + 1])
    if "--friction" in sys.argv[2:]:
        FRICTION = float(sys.argv[sys.argv.index("--friction") + 1])
    NOBLEND = "--noblend" in sys.argv[2:]
    VOIDFLOW = "--novoidflow" not in sys.argv[2:]
    STABLE = "--stable" in sys.argv[2:]
    if "--x0" in sys.argv[2:]:
        X0 = float(sys.argv[sys.argv.index("--x0") + 1])
    if "--gap" in sys.argv[2:]:
        GAP = float(sys.argv[sys.argv.index("--gap") + 1])
    if "--width" in sys.argv[2:]:
        # layer width [m]; the footing stays centred (decision 0095)
        w = float(sys.argv[sys.argv.index("--width") + 1])
        X0 += 0.5 * (w - NX * CELL)
        NX = int(round(w / CELL))
    if "--depth" in sys.argv[2:]:
        d = float(sys.argv[sys.argv.index("--depth") + 1])
        NY += int(round((d - DEPTH) / CELL))
        DEPTH = d
    if "--smax" in sys.argv[2:]:
        S_MAX = float(sys.argv[sys.argv.index("--smax") + 1]) * B
    if "--speed" in sys.argv[2:]:
        SPEED = float(sys.argv[sys.argv.index("--speed") + 1])
    if "--E" in sys.argv[2:]:
        E_SAND = float(sys.argv[sys.argv.index("--E") + 1])
    if "--maxco" in sys.argv[2:]:
        MAXCO = float(sys.argv[sys.argv.index("--maxco") + 1])
    t0 = time.perf_counter()
    ok = check_sand(work)
    print(f"\n{'A.4' if STABLE else 'WP 4.4c'} checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t0:.0f} s)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
