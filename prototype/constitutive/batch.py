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
Drucker-Prager. Mohr-Coulomb is not vectorised (eigen-decomposition per point).
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


MODELS = {
    "linearElastic": linear_elastic,
    "vonMises": von_mises,
    "DruckerPrager": drucker_prager,
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


if __name__ == "__main__":
    for k, v in check_batch().items():
        print(f"{k}: largest difference {v:.2e}")
