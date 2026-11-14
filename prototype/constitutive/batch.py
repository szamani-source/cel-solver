"""Vectorised versions of the WP 1.4 models, for the 2D prototype (WP 1.3b).

The 2D prototype calls a model for every cell in every sweep (tens of
thousands of calls per step), which is too slow one point at a time. These
functions do the same return mappings as ``models.py`` on arrays of points:

    sig_new, eps_p_new, eps_p_eq_new, plastic = model(deps, sig, eps_p, eps_p_eq, params)

- deps, sig, eps_p: arrays (N, 6), Voigt order (xx, yy, zz, xy, yz, zx),
  engineering shear strains, tensor shear stresses (decision 0004).
- eps_p_eq: array (N,). plastic: boolean array (N,), True where the point
  yielded in this increment.

Same yield tolerance, same corrector formulas, same plastic strain recovery
(C^-1 : (trial - new)). ``check_batch()`` compares them with ``models.py`` on
random paths (run by tests/test_proto_2d_rotation.sh through lr2d_plastic.py).
Only the models the 2D prototype uses are here: linear elastic, von Mises,
Drucker-Prager, and (WP 3.1p) Mohr-Coulomb. Mohr-Coulomb does the main-plane,
edge and apex returns of ``models.mohr_coulomb`` on stacked eigen-decompositions;
the few points that need the tension cut-off return (WP 2.5c active-set search)
are handed to the point model one by one, so results match it in every case.
"""

import math

import numpy as np

import models

I_V = models.I_V


def _elastic_parts(params):
    K, G = models.bulk_shear(params)
    D = models.elastic_matrix(K, G)
    C = models.compliance_matrix(K, G)
    return K, G, D, C


def _mean(s):
    return (s[:, 0] + s[:, 1] + s[:, 2]) / 3.0


def _j2(s):
    d = s - _mean(s)[:, None] * I_V
    return 0.5 * (d[:, :3] ** 2).sum(axis=1) + (d[:, 3:] ** 2).sum(axis=1)


def _eq_increment(deps_p):
    """sqrt(2/3 e_dev : e_dev) of plastic strain increments (as models._finish)."""
    e = deps_p.copy()
    e[:, 3:] *= 0.5
    dev = e - ((e[:, 0] + e[:, 1] + e[:, 2]) / 3.0)[:, None] * I_V
    n2 = (dev[:, :3] ** 2).sum(axis=1) + 2.0 * (dev[:, 3:] ** 2).sum(axis=1)
    return np.sqrt(2.0 / 3.0 * n2)


def linear_elastic(deps, sig, eps_p, eps_p_eq, params):
    _, _, D, _ = _elastic_parts(params)
    return sig + deps @ D.T, eps_p, eps_p_eq, np.zeros(len(sig), dtype=bool)


def von_mises(deps, sig, eps_p, eps_p_eq, params):
    """Radial return with linear isotropic hardening H (models.von_mises)."""
    K, G, D, C = _elastic_parts(params)
    H = params.get("H", 0.0)
    trial = sig + deps @ D.T
    q = np.sqrt(3.0 * _j2(trial))
    sy = params["sigma_y"] + H * eps_p_eq
    f = q - sy
    plastic = f > models.YIELD_TOL * sy
    if not plastic.any():
        return trial, eps_p, eps_p_eq, plastic
    new = trial.copy()
    ep = eps_p.copy()
    epq = eps_p_eq.copy()
    t = trial[plastic]
    dgam = f[plastic] / (3.0 * G + H)
    p = _mean(t)
    dev = t - p[:, None] * I_V
    s = dev * (1.0 - 3.0 * G * dgam / q[plastic])[:, None]
    new[plastic] = s + p[:, None] * I_V
    ep[plastic] = eps_p[plastic] + (t - new[plastic]) @ C.T
    epq[plastic] = eps_p_eq[plastic] + dgam
    return new, ep, epq, plastic


def drucker_prager(deps, sig, eps_p, eps_p_eq, params):
    """Cone return, apex when the cone return overshoots (models.drucker_prager)."""
    K, G, D, C = _elastic_parts(params)
    eta, xi, eta_bar = models.dp_coefficients(params)
    c = params["c"]
    trial = sig + deps @ D.T
    p = _mean(trial)
    sqj = np.sqrt(_j2(trial))
    f = sqj + eta * p - xi * c
    scale = np.maximum(np.maximum(np.abs(p), sqj), max(c, 1.0))
    plastic = f > models.YIELD_TOL * scale
    if not plastic.any():
        return trial, eps_p, eps_p_eq, plastic
    new = trial.copy()
    t, pp, sq = trial[plastic], p[plastic], sqj[plastic]
    dgam = f[plastic] / (G + K * eta * eta_bar)
    cone = (sq - G * dgam >= 0.0) | (eta == 0.0)
    dev = t - pp[:, None] * I_V
    fac = np.where(sq > 0.0, 1.0 - G * dgam / np.where(sq > 0.0, sq, 1.0), 1.0)
    cone_new = dev * fac[:, None] + (pp - K * eta_bar * dgam)[:, None] * I_V
    apex_new = np.outer(np.full(len(t), xi * c / eta if eta != 0.0 else 0.0), I_V)
    new[plastic] = np.where(cone[:, None], cone_new, apex_new)
    deps_p = (t - new[plastic]) @ C.T
    ep = eps_p.copy()
    epq = eps_p_eq.copy()
    ep[plastic] = eps_p[plastic] + deps_p
    epq[plastic] = eps_p_eq[plastic] + _eq_increment(deps_p)
    return new, ep, epq, plastic


def _to_tensor(s):
    T = np.empty((len(s), 3, 3))
    T[:, 0, 0], T[:, 1, 1], T[:, 2, 2] = s[:, 0], s[:, 1], s[:, 2]
    T[:, 0, 1] = T[:, 1, 0] = s[:, 3]
    T[:, 1, 2] = T[:, 2, 1] = s[:, 4]
    T[:, 2, 0] = T[:, 0, 2] = s[:, 5]
    return T


def mohr_coulomb(deps, sig, eps_p, eps_p_eq, params):
    """Principal-space return of models.mohr_coulomb: main plane, the edge the
    main-plane return crosses, then the apex; same tolerances and branch order.
    Points whose Mohr-Coulomb return lies beyond the tension cut-off (only
    with params sigma_t) are returned by the point model."""
    K, G, D, C = _elastic_parts(params)
    c = params["c"]
    st = models.mc_tension_cutoff(params)
    sphi = math.sin(math.radians(params["phi"]))
    cphi = math.cos(math.radians(params["phi"]))
    spsi = math.sin(math.radians(params["psi"]))
    trial = sig + deps @ D.T
    n = len(trial)
    plastic = np.zeros(n, dtype=bool)
    if n == 0:
        return trial, eps_p, eps_p_eq, plastic
    w, V = np.linalg.eigh(_to_tensor(trial))
    w, V = w[:, ::-1], V[:, :, ::-1]                 # w0 >= w1 >= w2
    k = 2.0 * c * cphi
    scale = np.maximum(np.abs(w).max(axis=1), max(c, 1.0))
    f = (w[:, 0] - w[:, 2]) + (w[:, 0] + w[:, 2]) * sphi - k
    t = np.full(n, -np.inf) if st is None else w[:, 0] - st
    plastic = (f > models.YIELD_TOL * scale) | (t > models.YIELD_TOL * scale)
    if not plastic.any():
        return trial, eps_p, eps_p_eq, plastic

    Dp = np.full((3, 3), K - 2.0 * G / 3.0) + np.eye(3) * 2.0 * G
    main, edge_12, edge_23 = models._mc_principal_planes(sphi, spsi)
    idx = np.nonzero(plastic)[0]
    wp = w[idx]
    tol = 1e-12 * scale[idx]
    s_new = np.full_like(wp, np.nan)
    apex = np.zeros(len(idx), dtype=bool)
    mc = f[idx] > models.YIELD_TOL * scale[idx]

    def ret(rows, planes):
        A = np.array([[nf @ Dp @ ng for (_, ng) in planes] for (nf, _) in planes])
        b = np.stack([wp[rows] @ nf - k for (nf, _) in planes], axis=1)
        if len(planes) == 1:
            dg = b / A[0, 0]
        else:
            dg = np.linalg.solve(np.broadcast_to(A, (len(b),) + A.shape), b[..., None])[..., 0]
        corr = sum(dg[:, j:j + 1] * planes[j][1][None, :] for j in range(len(planes)))
        return wp[rows] - corr @ Dp.T, dg

    def ordered(s, tl):
        return (s[:, 0] >= s[:, 1] - tl) & (s[:, 1] >= s[:, 2] - tl)

    r = np.nonzero(mc)[0]
    if len(r):
        s1, _ = ret(r, [main])
        s_new[r] = s1
        bad = ~ordered(s1, tol[r])
        for use23 in (True, False):
            sel = bad & ((s1[:, 1] > s1[:, 0]) == use23)
            rr = r[sel]
            if not len(rr):
                continue
            s2, dg = ret(rr, [main, edge_23 if use23 else edge_12])
            ok = np.all(dg >= -tol[rr, None], axis=1) & ordered(s2, tol[rr])
            if not ok.all() and sphi == 0.0:
                raise RuntimeError("Tresca corner return failed")
            s_new[rr[ok]] = s2[ok]
            if sphi > 0.0:
                s_new[rr[~ok]] = c * cphi / sphi
                apex[rr[~ok]] = True
        if sphi > 0.0:
            beyond = (~apex[r]) & ((s_new[r, 0] + s_new[r, 2]) / 2.0 > c * cphi / sphi + tol[r])
            s_new[r[beyond]] = c * cphi / sphi
            apex[r[beyond]] = True

    new = trial.copy()
    Vp = V[idx]
    new[idx] = _from_principal(Vp, np.nan_to_num(s_new))
    ep = eps_p.copy()
    epq = eps_p_eq.copy()
    if st is not None:
        tens = (~mc) | (s_new[:, 0] > st + tol)
        for j in idx[tens]:
            h = {"eps_p": eps_p[j], "eps_p_eq": eps_p_eq[j], "apex": False}
            s, hn = models.mohr_coulomb(deps[j], sig[j], h, params)
            new[j] = s
    # Plastic strain from (trial - new) for every plastic point (models._finish)
    deps_p = (trial[idx] - new[idx]) @ C.T
    ep[idx] = eps_p[idx] + deps_p
    epq[idx] = eps_p_eq[idx] + _eq_increment(deps_p)
    return new, ep, epq, plastic


def _from_principal(V, s):
    """Voigt stresses of V diag(s) V^T for stacks V (N, 3, 3), s (N, 3)."""
    T = (V * s[:, None, :]) @ np.swapaxes(V, 1, 2)
    return np.stack([T[:, 0, 0], T[:, 1, 1], T[:, 2, 2], T[:, 0, 1], T[:, 1, 2],
                     T[:, 2, 0]], axis=1)


MODELS = {
    "linearElastic": linear_elastic,
    "vonMises": von_mises,
    "DruckerPrager": drucker_prager,
    "MohrCoulomb": mohr_coulomb,
}


def check_batch(npaths=200, nsteps=30, seed=1):
    """Run random strain paths through the batch and the point models.

    Returns the largest stress difference relative to the stress scale, per
    model. Paths are chosen to cross the yield surface (and the DP apex)."""
    rng = np.random.default_rng(seed)
    cases = {
        "linearElastic": dict(E=20.0e6, nu=0.3),
        "vonMises": dict(E=20.0e6, nu=0.3, sigma_y=10.0e3, H=2.0e5),
        "DruckerPrager": dict(E=20.0e6, nu=0.25, phi=30.0, c=5.0e3, psi=10.0,
                              fit="compression"),
    }
    out = {}
    for name, params in cases.items():
        point = models.MODELS[name]
        batch = MODELS[name]
        sig = np.zeros((npaths, 6))
        ep = np.zeros((npaths, 6))
        epq = np.zeros(npaths)
        ref = [(np.zeros(6), models.init_history()) for _ in range(npaths)]
        worst = 0.0
        for _ in range(nsteps):
            deps = rng.normal(0.0, 2.0e-4, (npaths, 6))
            deps[:, :3] += rng.normal(0.0, 1.0e-4, (npaths, 1))   # volumetric part
            sig, ep, epq, _ = batch(deps, sig, ep, epq, params)
            for k in range(npaths):
                s, h = point(deps[k], ref[k][0], ref[k][1], params)
                ref[k] = (s, h)
                scale = max(np.abs(s).max(), 1.0e3)
                d = max(np.abs(sig[k] - s).max() / scale,
                        np.abs(ep[k] - h["eps_p"]).max() / 1e-3,
                        abs(epq[k] - h["eps_p_eq"]) / 1e-3)
                worst = max(worst, d)
        out[name] = worst
    return out


MC_CASES = {
    # T0 bulk solid (docs/spec/target-problems.md): c = 0, psi = 0, cut-off 0
    "T0 c=0 phi=35 sigma_t=0": (dict(E=10.0e6, nu=0.3, phi=35.0, c=0.0, psi=0.0,
                                     sigma_t=0.0), -1.0e-4),
    "c=5kPa phi=30 psi=10": (dict(E=20.0e6, nu=0.25, phi=30.0, c=5.0e3, psi=10.0), 0.0),
    "c=5kPa phi=30 sigma_t=2kPa": (dict(E=20.0e6, nu=0.25, phi=30.0, c=5.0e3, psi=0.0,
                                        sigma_t=2.0e3), 0.0),
}


def check_batch_mc(npaths=200, nsteps=30, seed=2):
    """Batch Mohr-Coulomb against models.mohr_coulomb on random paths (WP 3.1p).

    Returns {case: (largest difference, counts)} with counts of elastic,
    main-plane, edge and apex returns of the point model (from its history
    and a re-check of the returned stress), so the test can require each branch
    to be exercised. The difference is measured as in check_batch."""
    rng = np.random.default_rng(seed)
    out = {}
    for name, (params, bias) in MC_CASES.items():
        sig = np.zeros((npaths, 6))
        ep = np.zeros((npaths, 6))
        epq = np.zeros(npaths)
        ref = [(np.zeros(6), models.init_history()) for _ in range(npaths)]
        worst = 0.0
        counts = {"elastic": 0, "plastic": 0, "apex": 0, "edge": 0}
        for _ in range(nsteps):
            deps = rng.normal(0.0, 2.0e-4, (npaths, 6))
            deps[:, :3] += rng.normal(bias, 1.0e-4, (npaths, 1))   # volumetric part
            sig, ep, epq, _ = mohr_coulomb(deps, sig, ep, epq, params)
            for k in range(npaths):
                s, h = models.mohr_coulomb(deps[k], ref[k][0], ref[k][1], params)
                yielded = h["eps_p_eq"] > ref[k][1]["eps_p_eq"]
                counts["plastic" if yielded else "elastic"] += 1
                if h["apex"]:
                    counts["apex"] += 1
                elif yielded:
                    p = np.sort(np.linalg.eigvalsh(models.stress_to_tensor(s)))
                    gap = min(p[2] - p[1], p[1] - p[0])
                    if gap <= 1e-9 * max(np.abs(p).max(), 1.0):
                        counts["edge"] += 1
                ref[k] = (s, h)
                scale = max(np.abs(s).max(), 1.0e3)
                d = max(np.abs(sig[k] - s).max() / scale,
                        np.abs(ep[k] - h["eps_p"]).max() / 1e-3,
                        abs(epq[k] - h["eps_p_eq"]) / 1e-3)
                worst = max(worst, d)
        out[name] = (worst, counts)
    return out


if __name__ == "__main__":
    for k, v in check_batch().items():
        print(f"{k}: largest difference {v:.2e}")
    for k, (v, n) in check_batch_mc().items():
        print(f"MohrCoulomb {k}: largest difference {v:.2e}, {n}")
