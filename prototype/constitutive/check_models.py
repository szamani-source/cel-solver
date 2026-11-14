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

# 7. Mohr-Coulomb tension cut-off sigma_t (WP 2.5c)
ccot = c * co / s                                   # apex c cot(phi), 17.3 kPa
mc_t_strength = 2 * c * co / (1 + s)                # MC uniaxial tensile strength, 11.5 kPa
for frac in (0.0, 0.5, 0.9):
    st = frac * ccot
    p = dict(mc, sigma_t=st)
    # Uniaxial tension (radial stress zero): strength min(sigma_t, 2 c cos(phi)/(1 + sin(phi)))
    r = driver.triaxial(models.mohr_coulomb, p, 0.0, 0.01, 200, extension=True)
    target = min(st, mc_t_strength)
    check("MC sigma_t = %.2f c cot(phi): uniaxial tensile strength" % frac,
          abs(r["sigma"][-1, 2] - target) < 1e-9 * c,
          "%.6g Pa, expected %.6g Pa" % (r["sigma"][-1, 2], target))
    # Hydrostatic tension from 100 kPa compression is capped at sigma_t
    sig, hist = -pc * models.I_V, models.init_history()
    for _ in range(100):
        sig, hist = models.mohr_coulomb(1e-4 * models.I_V, sig, hist, p)
    check("MC sigma_t = %.2f c cot(phi): hydrostatic tension capped at sigma_t" % frac,
          np.max(np.abs(sig - st * models.I_V)) < 1e-9 * c and hist["apex"],
          "sigma = %s Pa, expected %.6g Pa, apex %s" % (np.array2string(sig[:3], precision=6),
                                                       st, hist["apex"]))

# Unloading from the cut-off is elastic
p = dict(mc, sigma_t=0.5 * ccot)
r = driver.triaxial(models.mohr_coulomb, p, 0.0, 0.01, 50, extension=True)
sig, hist = r["sigma"][-1], r["history"][-1]
back = np.array([0.0, 0.0, -1e-5, 0, 0, 0])
sig2, hist2 = models.mohr_coulomb(back, sig, hist, p)
K, G = models.bulk_shear(p)
check("MC tension cut-off: unloading is elastic",
      np.array_equal(hist2["eps_p"], hist["eps_p"])
      and np.allclose(sig2, sig + models.elastic_matrix(K, G) @ back, rtol=0, atol=1e-9 * c),
      "plastic strain unchanged, sigma_zz %.6g -> %.6g Pa" % (sig[2], sig2[2]))

# Compression paths never reach the cut-off and are unchanged bit for bit
same = True
for fn in (lambda q: driver.triaxial(models.mohr_coulomb, q, pc, -0.05, 200),
           lambda q: driver.triaxial(models.mohr_coulomb, q, pc, 0.05, 200, extension=True),
           lambda q: driver.simple_shear(models.mohr_coulomb, q, pc, 0.2, 200)):
    a, b = fn(mc), fn(dict(mc, sigma_t=0.0))
    same = same and np.array_equal(a["sigma"], b["sigma"]) and np.array_equal(a["eps"], b["eps"])
check("MC tension cut-off: compression paths unchanged (sigma_t = 0)", same,
      "triaxial compression, extension and simple shear from 100 kPa")

# Out-of-range sigma_t is rejected
rejected = 0
for bad in (-1.0, 1.01 * ccot):
    try:
        models.mohr_coulomb(np.zeros(6), np.zeros(6), models.init_history(),
                            dict(mc, sigma_t=bad))
    except ValueError:
        rejected += 1
check("MC sigma_t outside [0, c cot(phi)] rejected", rejected == 2, "%d of 2" % rejected)


def tension_increment():
    mag = 10 ** rng.uniform(-6, -2)
    kind = rng.integers(4)
    if kind == 0:   # volumetric stretch: tension apex
        return mag * models.I_V * rng.uniform(0.5, 2)
    if kind == 1:   # two equal principal stretches: tension edge s1 = s2
        a, b = rng.uniform(0, 1), rng.uniform(-1, 1)
        return mag * np.array([a, a, b, 0, 0, 0])
    if kind == 2:   # uniaxial-like stretch: tension plane and mixed corners
        a, b = rng.uniform(0.5, 1), rng.uniform(-0.5, 0)
        return mag * np.array([a, b, b, 0, 0, 0])
    return mag * (rng.normal(size=6) + 0.5 * models.I_V)


for name, p in (("sigma_t = 0", dict(mc, sigma_t=0.0)),
                ("sigma_t = 0.5 c cot(phi)", dict(mc, sigma_t=0.5 * ccot)),
                ("sigma_t = 0.5 c cot(phi), associated", dict(mc, psi=phi, sigma_t=0.5 * ccot)),
                ("sigma_t = 0.2 c cot(phi), psi = 0", dict(mc, psi=0.0, sigma_t=0.2 * ccot))):
    worst, max_s1, n_cut, n_apex = -np.inf, -np.inf, 0, 0
    sig, hist = -pc * models.I_V, models.init_history()
    for i in range(n_per):
        if i % 50 == 0:  # restart from a random compressive state
            sig = -pc * rng.uniform(0.0, 1.0) * models.I_V
        deps = tension_increment() if rng.random() < 0.7 else random_increment()
        sig_new, hist_new = models.mohr_coulomb(deps, sig, hist, p)
        scale = max(np.max(np.abs(sig_new)), p["c"], 1.0)
        ft = models.rankine_yield(sig_new, p)
        worst = max(worst, models.mc_yield(sig_new, p) / scale, ft / scale)
        max_s1 = max(max_s1, ft + p["sigma_t"])
        if np.any(hist_new["eps_p"] != hist["eps_p"]) and abs(ft) < 1e-9 * scale:
            n_cut += 1
        n_apex += hist_new["apex"]
        sig, hist = sig_new, hist_new
    check("MC %s: 10000 random increments inside MC and Rankine surfaces" % name,
          worst < 1e-9 and n_cut > 1000 and n_apex > 100,
          "max f/scale = %.2e, max s1 = %.4g Pa, cut-off returns %d, apex %d"
          % (worst, max_s1, n_cut, n_apex))

# Associated flow: the return is the closest point in the energy norm (independent check:
# project onto every face of the full 9-plane surface in unsorted principal space and keep
# the nearest admissible point)
import itertools  # noqa: E402

p = dict(mc, psi=phi, sigma_t=0.4 * ccot)
K, G = models.bulk_shear(p)
Dp = np.full((3, 3), K - 2 * G / 3) + np.eye(3) * 2 * G
Cp = np.linalg.inv(Dp)
N, rr = [], []
for i, j in itertools.permutations(range(3), 2):
    n = np.zeros(3)
    n[i] += 1 + s
    n[j] += -1 + s
    N.append(n)
    rr.append(2 * c * co)
for i in range(3):
    N.append(np.eye(3)[i])
    rr.append(p["sigma_t"])
N, rr = np.array(N), np.array(rr)


def closest_point(w):
    if np.all(N @ w - rr <= 0):
        return w
    best = None
    for size in (1, 2, 3):
        for act in itertools.combinations(range(9), size):
            A, b = N[list(act)], rr[list(act)]
            if np.linalg.matrix_rank(A) < size:
                continue
            x = w - Dp @ A.T @ np.linalg.solve(A @ Dp @ A.T, A @ w - b)
            if np.all(N @ x - rr <= 1e-9 * c):
                d = (x - w) @ Cp @ (x - w)
                if best is None or d < best[0]:
                    best = (d, x)
    return best[1]


worst = 0.0
for _ in range(1000):
    w = np.sort(rng.normal(size=3) * 3e4 + rng.uniform(-3e4, 4e4))[::-1]
    if rng.random() < 0.2:
        w[1] = w[0]
    if rng.random() < 0.2:
        w[2] = w[1]
    sig, _ = models.mohr_coulomb(np.zeros(6), np.r_[w, 0, 0, 0], models.init_history(), p)
    got = np.sort(sig[:3])[::-1]
    want = np.sort(closest_point(w))[::-1]
    worst = max(worst, np.max(np.abs(got - want)) / max(np.max(np.abs(w)), c))
check("MC tension cut-off, associated: return = closest point (1000 trial states)",
      worst < 1e-12, "max relative difference %.2e" % worst)

print()
if failures:
    print("%d check(s) failed: %s" % (len(failures), ", ".join(failures)))
    sys.exit(1)
print("All checks passed.")
