"""WP 3.6a checks of the mu(I) friction (Python reference, decision 0037).

1. Steady simple shear at fixed pressure: all three normal stresses held at -P, the
   shear strain xy applied at a constant rate gamma_dot for 300 steps. The friction
   mobilised at the end, read from the stress alone (Drucker-Prager compression fit:
   sqrt(J2)/P -> phi_m; Mohr-Coulomb: sin phi_m = (s1 - s3) / -(s1 + s3)), must equal
   mu(I) at I = gamma_dot d / sqrt(P / rho_s) to a relative 1e-6, for I = 1e-4 to 1.
   Cases: placeholder sand (mu_s = tan 35 deg, plain Jop law) with Drucker-Prager
   (psi = 0 and 10 deg) and Mohr-Coulomb (psi = 0), and Jop's glass beads with the
   Barker-Gray branch (alpha = 1.9, I_1^N = 0.019, so I = 1e-4 and 1e-3 are on it).
2. mu2 = tan(phi): the reference paths of the C++ kernel test (uniaxial strain,
   triaxial compression and extension, simple shear, 1000 random increments) give the
   rate-independent model bit for bit.

Usage: python3 check_mu_i.py   (exit status 1 on failure)
"""

import math
import sys

import numpy as np

import driver
import models
import mu_i

E, NU = 10e6, 0.3
P = 10e3


def measured_mu(model_name, sigma):
    if model_name == "dp":
        eta = math.sqrt(models.j2_invariant(sigma)) / -models.mean_stress(sigma)
        s = 3.0 * math.sqrt(3.0) * eta / (6.0 + math.sqrt(3.0) * eta)
    else:
        w = np.linalg.eigvalsh(models.stress_to_tensor(sigma))
        s = (w[-1] - w[0]) / -(w[-1] + w[0])
    return math.tan(math.asin(s))


def steady_shear(model, params, I):
    gdot = I * math.sqrt(P / params["rho_s"]) / params["d"]
    dgam = 1e-4
    p = dict(params, dt=dgam / gdot)
    s0 = -P * models.I_V
    d = np.zeros(6)
    d[3] = dgam
    run = driver.run_path(model, p, s0, 300, d, [False, False, False, True, True, True], s0)
    sig = run["sigma"]
    return sig[-1], np.max(np.abs(sig[-1] - sig[-2])) / P


def check_steady():
    sand = dict(mu_i.sand_placeholder(), E=E, nu=NU, c=0.0)
    glass = dict(phi=20.9,
                 mu2=math.tan(math.radians(32.76)), I0=0.279, d=0.53e-3, rho_s=2450.0,
                 muI_alpha=mu_i.ALPHA_BARKER_GRAY, E=E, nu=NU, c=0.0)
    cases = (
        ("sand DP psi 0", "dp", mu_i.mu_i_drucker_prager, dict(sand, psi=0.0)),
        ("sand DP psi 10", "dp", mu_i.mu_i_drucker_prager, dict(sand, psi=10.0)),
        ("sand MC psi 0", "mc", mu_i.mu_i_mohr_coulomb, dict(sand, psi=0.0)),
        ("glass BG DP psi 0", "dp", mu_i.mu_i_drucker_prager, dict(glass, psi=0.0)),
        ("glass BG MC psi 0", "mc", mu_i.mu_i_mohr_coulomb, dict(glass, psi=0.0)),
    )
    worst = 0.0
    ok = True
    for name, kind, model, p in cases:
        law = mu_i.friction_law(p)
        extra = "" if law.alpha == 0.0 else "  (I_1^N %.4g)" % law.I1N
        print("%s%s" % (name, extra))
        for I in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
            sig, change = steady_shear(model, p, I)
            mu_m = measured_mu(kind, sig)
            mu_t = law.mu(I)
            err = abs(mu_m - mu_t) / mu_t
            worst = max(worst, err)
            good = err <= 1e-6 and change <= 1e-9
            ok = ok and good
            print("  I %.0e  mu(I) %.9f  measured %.9f  rel err %.1e  last change %.1e  %s"
                  % (I, mu_t, mu_m, err, change, "ok" if good else "FAIL"))
    print("Steady simple shear: worst relative error %.2e (limit 1e-6)" % worst)
    return ok


def check_bit_for_bit():
    rng = np.random.default_rng(2022)
    steps = [10 ** rng.uniform(-6, -2.5) * rng.normal(size=6) for _ in range(1000)]
    pc = 100e3
    sets = (
        ("dpCohesive", models.drucker_prager, mu_i.mu_i_drucker_prager,
         dict(E=20e6, nu=0.3, phi=30.0, c=10e3, psi=5.0, fit="compression")),
        ("dpSand", models.drucker_prager, mu_i.mu_i_drucker_prager,
         dict(E=20e6, nu=0.3, phi=35.0, c=0.0, psi=0.0, fit="compression")),
        ("mcCohesive", models.mohr_coulomb, mu_i.mu_i_mohr_coulomb,
         dict(E=20e6, nu=0.3, phi=30.0, c=10e3, psi=5.0)),
        ("mcSand", models.mohr_coulomb, mu_i.mu_i_mohr_coulomb,
         dict(E=20e6, nu=0.3, phi=35.0, c=0.0, psi=0.0)),
        ("mcTension", models.mohr_coulomb, mu_i.mu_i_mohr_coulomb,
         dict(E=20e6, nu=0.3, phi=30.0, c=10e3, psi=5.0, sigma_t=5e3)),
    )
    ok = True
    for tag, base, rd, p in sets:
        q = dict(p, mu2=math.tan(math.radians(p["phi"])), I0=0.279, d=1e-3,
                 rho_s=2650.0, dt=1e-4)
        paths = {
            "triaxialCompression": driver.triaxial(base, p, pc, -0.02, 200),
            "triaxialExtension": driver.triaxial(base, p, pc, 0.02, 200, extension=True),
            "simpleShear": driver.simple_shear(base, p, pc, 0.05, 200),
        }
        incs = {k: np.diff(r["eps"], axis=0) for k, r in paths.items()}
        d = np.zeros(6)
        incs["uniaxialStrain"] = ([d + [0, 0, -2e-4, 0, 0, 0]] * 100
                                  + [d + [0, 0, 2e-4, 0, 0, 0]] * 150)
        incs["random1000"] = steps
        n_same = 0
        for name, dl in incs.items():
            s1 = s2 = -pc * models.I_V
            h1 = h2 = models.init_history()
            same = True
            for de in dl:
                s1, h1 = base(de, s1, h1, p)
                s2, h2 = rd(de, s2, h2, q)
                same = same and np.array_equal(s1, s2) and np.array_equal(
                    h1["eps_p"], h2["eps_p"]) and h1["eps_p_eq"] == h2["eps_p_eq"] \
                    and h1["apex"] == h2["apex"]
            n_same += same
            ok = ok and same
        print("  %-11s mu2 = tan(phi): %d of %d paths identical" % (tag, n_same, len(incs)))
    print("mu2 = mu_s bit for bit: %s" % ("ok" if ok else "FAIL"))
    return ok


if __name__ == "__main__":
    a = check_steady()
    b = check_bit_for_bit()
    print("PASS" if a and b else "FAIL")
    sys.exit(0 if a and b else 1)
