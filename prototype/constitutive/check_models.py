"""WP 1.4 verification checks (run by tests/test_proto_constitutive.sh).

Each check prints one line with the computed and expected values and exits non-zero
on the first failure. Units: Pa, stress positive in tension.
"""

import math
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__file__))

import driver  # noqa: E402
import models  # noqa: E402

failures = []


def check(name, ok, detail):
    print(("PASS  " if ok else "FAIL  ") + name + ": " + detail)
    if not ok:
        failures.append(name)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-30)


# Reference sand-like material: E = 20 MPa, nu = 0.3, phi = 30 deg, c = 10 kPa, psi = 10 deg
E, nu, phi, c, psi = 20e6, 0.3, 30.0, 10e3, 10.0
pc = 100e3  # confining pressure
s, co = math.sin(math.radians(phi)), math.cos(math.radians(phi))

mc = dict(E=E, nu=nu, phi=phi, c=c, psi=psi)
dpc = dict(mc, fit="compression")
dpe = dict(mc, fit="extension")
dpp = dict(mc, fit="plane_strain")

# 1. Linear elastic, uniaxial strain against closed form
K, G = models.bulk_shear(dict(E=E, nu=nu))
r = driver.uniaxial_strain(models.linear_elastic, dict(E=E, nu=nu), -1e-3, 10)
sz, sx = r["sigma"][-1, 2], r["sigma"][-1, 0]
check("elastic uniaxial strain", rel(sz, -(K + 4 * G / 3) * 1e-3) < 1e-12
      and rel(sx, -(K - 2 * G / 3) * 1e-3) < 1e-12,
      "sigma_zz=%.6g sigma_xx=%.6g" % (sz, sx))

# 2. J2: uniaxial stress plateau at sigma_y (perfect) and linear hardening slope
j2 = dict(E=E, nu=nu, sigma_y=200e3, H=0.0)
r = driver.triaxial(models.von_mises, j2, 0.0, -0.02, 200, extension=False)
check("J2 uniaxial yield (perfect plasticity)", rel(-r["sigma"][-1, 2], 200e3) < 1e-9,
      "sigma_zz=%.6g expected -2e5" % r["sigma"][-1, 2])
H = 2e6
r = driver.triaxial(models.von_mises, dict(j2, H=H), 0.0, -0.02, 200, extension=False)
epq = r["history"][-1]["eps_p_eq"]
check("J2 linear hardening", rel(-r["sigma"][-1, 2], 200e3 + H * epq) < 1e-9,
      "sigma_zz=%.6g expected %.6g" % (r["sigma"][-1, 2], -(200e3 + H * epq)))

# 3. Mohr-Coulomb triaxial compression and extension strengths
sa_c = (2 * c * co + pc * (1 + s)) / (1 - s)        # axial compressive stress at failure
sa_e = (pc * (1 - s) - 2 * c * co) / (1 + s)        # axial compressive stress at failure
r = driver.triaxial(models.mohr_coulomb, mc, pc, -0.05, 500, extension=False)
check("MC triaxial compression strength", rel(-r["sigma"][-1, 2], sa_c) < 1e-3,
      "%.6g Pa, closed form %.6g Pa, error %.2e" % (-r["sigma"][-1, 2], sa_c, rel(-r["sigma"][-1, 2], sa_c)))
r = driver.triaxial(models.mohr_coulomb, mc, pc, 0.05, 500, extension=True)
check("MC triaxial extension strength", rel(-r["sigma"][-1, 2], sa_e) < 1e-3,
      "%.6g Pa, closed form %.6g Pa, error %.2e" % (-r["sigma"][-1, 2], sa_e, rel(-r["sigma"][-1, 2], sa_e)))

# Cohesionless MC compression, and pure tension of cohesionless material drops to zero stress
r = driver.triaxial(models.mohr_coulomb, dict(mc, c=0.0), pc, -0.05, 500, extension=False)
target = pc * (1 + s) / (1 - s)
check("MC cohesionless compression (Kp * pc)", rel(-r["sigma"][-1, 2], target) < 1e-3,
      "%.6g Pa, closed form %.6g Pa" % (-r["sigma"][-1, 2], target))
sig, _ = models.mohr_coulomb(np.array([1e-2, 1e-2, 1e-2, 0, 0, 0]), -pc * models.I_V,
                             models.init_history(), dict(mc, c=0.0))
check("MC cohesionless apex return under tension", np.max(np.abs(sig)) < 1e-6,
      "max |sigma| = %.3g Pa" % np.max(np.abs(sig)))

# 4. Drucker-Prager fitted to the compression meridian gives the MC compression strength;
#    extension fit gives the MC extension strength
r = driver.triaxial(models.drucker_prager, dpc, pc, -0.05, 500, extension=False)
check("DP (compression fit) triaxial compression strength", rel(-r["sigma"][-1, 2], sa_c) < 1e-3,
      "%.6g Pa, MC closed form %.6g Pa" % (-r["sigma"][-1, 2], sa_c))
r = driver.triaxial(models.drucker_prager, dpe, pc, 0.05, 500, extension=True)
check("DP (extension fit) triaxial extension strength", rel(-r["sigma"][-1, 2], sa_e) < 1e-3,
      "%.6g Pa, MC closed form %.6g Pa" % (-r["sigma"][-1, 2], sa_e))

# 5. 10,000 random strain increments: returned stress on or inside the yield surface
rng = np.random.default_rng(1234)


def random_increment():
    mag = 10 ** rng.uniform(-6, -2)
    kind = rng.integers(5)
    if kind == 0:   # volumetric stretch: drives the stress to the apex
        return mag * models.I_V * rng.uniform(0.5, 2)
    if kind == 1:   # triaxial compression-like: two equal principal strains (corner)
        a, b = rng.uniform(-1, 1, 2)
        return mag * np.array([a, a, b, 0, 0, 0])
    if kind == 2:   # triaxial extension-like
        a, b = rng.uniform(-1, 1, 2)
        return mag * np.array([a, b, b, 0, 0, 0])
    return mag * rng.normal(size=6)


cases = [
    ("MC", models.mohr_coulomb, mc, models.mc_yield),
    ("MC c=0", models.mohr_coulomb, dict(mc, c=0.0), models.mc_yield),
    ("MC associated", models.mohr_coulomb, dict(mc, psi=phi), models.mc_yield),
    ("DP compression", models.drucker_prager, dpc, models.dp_yield),
    ("DP extension c=0", models.drucker_prager, dict(dpe, c=0.0), models.dp_yield),
    ("DP plane strain", models.drucker_prager, dpp, models.dp_yield),
]
n_per = 10000
for name, model, p, yf in cases:
    worst, n_apex, n_plastic = -np.inf, 0, 0
    sig, hist = -pc * models.I_V, models.init_history()
    for i in range(n_per):
        if i % 50 == 0:  # restart from a random compressive state
            sig = -pc * rng.uniform(0.01, 3) * models.I_V
        sig_new, hist_new = model(random_increment(), sig, hist, p)
        f = yf(sig_new, p)
        scale = max(np.max(np.abs(sig_new)), p["c"], 1.0)
        worst = max(worst, f / scale)
        if np.any(hist_new["eps_p"] != hist["eps_p"]):
            n_plastic += 1
        if hist_new.get("apex", False):
            n_apex += 1
        sig, hist = sig_new, hist_new
    check("%s: 10000 random increments on/inside yield surface" % name, worst < 1e-9,
          "max f/scale = %.2e, plastic steps %d, apex returns %d" % (worst, n_plastic, n_apex))
    if name in ("MC c=0", "DP extension c=0"):
        check("%s: apex cases exercised" % name, n_apex > 100, "%d apex returns" % n_apex)

# Corner returns in MC are exercised by the triaxial paths above (two equal principal stresses).

# 6. Dilatancy in simple shear at constant normal stress (plane strain, steady plastic flow)
sn = 100e3
r = driver.simple_shear(models.mohr_coulomb, mc, sn, 0.2, 2000)
ev, gam = r["eps"][:, 1], r["eps"][:, 3]
k = len(gam) // 2   # steady part of the curve
slope = (ev[-1] - ev[k]) / (gam[-1] - gam[k])
check("MC simple shear dilatancy d(eps_v)/d(gamma) = tan(psi)",
      rel(slope, math.tan(math.radians(psi))) < 1e-3,
      "%.6f, expected %.6f" % (slope, math.tan(math.radians(psi))))

r = driver.simple_shear(models.drucker_prager, dpp, sn, 0.2, 2000)
ev, gam = r["eps"][:, 1], r["eps"][:, 3]
slope = (ev[-1] - ev[k]) / (gam[-1] - gam[k])
check("DP (plane-strain fit) simple shear dilatancy = tan(psi)",
      rel(slope, math.tan(math.radians(psi))) < 1e-3,
      "%.6f, expected %.6f" % (slope, math.tan(math.radians(psi))))

r = driver.simple_shear(models.mohr_coulomb, dict(mc, psi=0.0), sn, 0.2, 2000)
ev = r["eps"][:, 1]
check("MC simple shear with psi = 0 has no plastic volume change",
      abs(ev[-1] - ev[k]) < 1e-9, "delta eps_v = %.2e" % (ev[-1] - ev[k]))

print()
if failures:
    print("%d check(s) failed: %s" % (len(failures), ", ".join(failures)))
    sys.exit(1)
print("All checks passed.")
