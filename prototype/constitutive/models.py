"""Single-point constitutive models for WP 1.4 (reference for the C++ port, WP 2.2/2.5).

Every model has the same signature

    sigma_new, history_new = model(deps, sigma, history, params)

- deps:    strain increment, Voigt order (xx, yy, zz, xy, yz, zx), engineering shear
           strains (gamma = 2 eps).
- sigma:   stress before the increment, same Voigt order (tensor shear components).
- history: dict from init_history(); never modified in place.
- params:  dict of material parameters (SI units, angles in degrees).

Stress is positive in tension. Small strain, perfect plasticity for Drucker-Prager and
Mohr-Coulomb (no softening or hardening), linear isotropic hardening for J2. The plastic
strain increment is recovered as C : (sigma_trial - sigma_new), which is exact for an
elastic predictor / plastic corrector scheme.
"""

import math

import numpy as np

I_V = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
_SHEAR = [(0, 1), (1, 2), (2, 0)]

# Relative tolerance on the trial yield function below which a step is taken as elastic.
YIELD_TOL = 1e-12


# ----------------------------------------------------------------------------- helpers

def init_history():
    return {"eps_p": np.zeros(6), "eps_p_eq": 0.0, "apex": False}


def bulk_shear(p):
    E, nu = p["E"], p["nu"]
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


def elastic_matrix(K, G):
    """6x6 stiffness mapping engineering strain to stress."""
    D = np.zeros((6, 6))
    D[:3, :3] = K - 2.0 * G / 3.0
    D[:3, :3] += np.eye(3) * 2.0 * G
    D[3:, 3:] = np.eye(3) * G
    return D


def compliance_matrix(K, G):
    return np.linalg.inv(elastic_matrix(K, G))


def stress_to_tensor(s):
    T = np.diag(s[:3]).astype(float)
    for k, (i, j) in enumerate(_SHEAR):
        T[i, j] = T[j, i] = s[3 + k]
    return T


def tensor_to_stress(T):
    return np.array([T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[1, 2], T[2, 0]])


def mean_stress(s):
    return (s[0] + s[1] + s[2]) / 3.0


def deviator(s):
    return s - mean_stress(s) * I_V


def j2_invariant(s):
    d = deviator(s)
    return 0.5 * (d[:3] @ d[:3]) + d[3:] @ d[3:]


def _finish(sig_trial, sig_new, history, K, G, apex=False):
    """Build the new history from the plastic corrector sig_trial - sig_new."""
    deps_p = compliance_matrix(K, G) @ (sig_trial - sig_new)
    e = deps_p.copy()
    e[3:] *= 0.5                     # tensor shear strains
    dev = e - (e[0] + e[1] + e[2]) / 3.0 * I_V
    dev_norm2 = dev[:3] @ dev[:3] + 2.0 * (dev[3:] @ dev[3:])
    h = dict(history)
    h["eps_p"] = history["eps_p"] + deps_p
    h["eps_p_eq"] = history["eps_p_eq"] + math.sqrt(2.0 / 3.0 * dev_norm2)
    h["apex"] = apex
    return sig_new, h


def _elastic(sig_trial, history):
    h = dict(history)
    h["apex"] = False
    return sig_trial, h


# ----------------------------------------------------------------------- linear elastic

def linear_elastic(deps, sigma, history, params):
    """params: E, nu."""
    K, G = bulk_shear(params)
    return _elastic(sigma + elastic_matrix(K, G) @ deps, history)


# ------------------------------------------------------------------------- von Mises J2

def von_mises(deps, sigma, history, params):
    """params: E, nu, sigma_y, H (linear isotropic hardening modulus, default 0).

    Radial return; exact for linear hardening.
    """
    K, G = bulk_shear(params)
    H = params.get("H", 0.0)
    trial = sigma + elastic_matrix(K, G) @ deps
    q = math.sqrt(3.0 * j2_invariant(trial))
    sy = params["sigma_y"] + H * history["eps_p_eq"]
    f = q - sy
    if f <= YIELD_TOL * sy:
        return _elastic(trial, history)
    dgam = f / (3.0 * G + H)
    p = mean_stress(trial)
    s = deviator(trial) * (1.0 - 3.0 * G * dgam / q)
    sig_new, h = _finish(trial, s + p * I_V, history, K, G)
    h["eps_p_eq"] = history["eps_p_eq"] + dgam   # exact value for J2
    return sig_new, h


def von_mises_yield(sigma, params, history=None):
    ep = 0.0 if history is None else history["eps_p_eq"]
    return math.sqrt(3.0 * j2_invariant(sigma)) - params["sigma_y"] - params.get("H", 0.0) * ep


# ----------------------------------------------------------------------- Drucker-Prager

def dp_coefficients(params):
    """Return (eta, xi, eta_bar) for f = sqrt(J2) + eta p - xi c, g = sqrt(J2) + eta_bar p.

    p is the mean stress (positive in tension). params["fit"] selects the match to
    Mohr-Coulomb: "compression" (outer cone, default), "extension" (inner cone) or
    "plane_strain" (de Souza Neto, Peric and Owen, section 6.4.3).
    """
    fit = params.get("fit", "compression")

    def coeffs(angle):
        a = math.radians(angle)
        sa, ca, ta = math.sin(a), math.cos(a), math.tan(a)
        if fit == "compression":
            d = math.sqrt(3.0) * (3.0 - sa)
            return 6.0 * sa / d, 6.0 * ca / d
        if fit == "extension":
            d = math.sqrt(3.0) * (3.0 + sa)
            return 6.0 * sa / d, 6.0 * ca / d
        if fit == "plane_strain":
            d = math.sqrt(9.0 + 12.0 * ta * ta)
            return 3.0 * ta / d, 3.0 / d
        raise ValueError("unknown Drucker-Prager fit '%s'" % fit)

    eta, xi = coeffs(params["phi"])
    eta_bar, _ = coeffs(params["psi"])
    return eta, xi, eta_bar


def drucker_prager(deps, sigma, history, params):
    """params: E, nu, phi, c, psi (degrees, Pa), fit (see dp_coefficients).

    Return mapping to the smooth cone, or to the apex when the cone return overshoots
    (de Souza Neto et al., Box 8.9, perfect plasticity).
    """
    K, G = bulk_shear(params)
    eta, xi, eta_bar = dp_coefficients(params)
    c = params["c"]
    trial = sigma + elastic_matrix(K, G) @ deps
    p, s = mean_stress(trial), deviator(trial)
    sqj = math.sqrt(j2_invariant(trial))
    f = sqj + eta * p - xi * c
    if f <= YIELD_TOL * max(abs(p), sqj, c, 1.0):
        return _elastic(trial, history)
    dgam = f / (G + K * eta * eta_bar)
    if sqj - G * dgam >= 0.0 or eta == 0.0:
        s_new = s * (1.0 - G * dgam / sqj) if sqj > 0.0 else s
        return _finish(trial, s_new + (p - K * eta_bar * dgam) * I_V, history, K, G)
    # Apex: zero deviator, mean stress at the cone tip (zero for cohesionless material)
    return _finish(trial, (xi * c / eta) * I_V, history, K, G, apex=True)


def dp_yield(sigma, params, history=None):
    eta, xi, _ = dp_coefficients(params)
    return math.sqrt(j2_invariant(sigma)) + eta * mean_stress(sigma) - xi * params["c"]


# ------------------------------------------------------------------------- Mohr-Coulomb

def _mc_principal_planes(sphi, spsi):
    """Normals (yield, flow) in sorted principal space s1 >= s2 >= s3 of the main plane
    and the two edge planes, for f = (smax - smin) + (smax + smin) sin(phi) - 2 c cos(phi).
    """
    main = (np.array([1 + sphi, 0.0, -1 + sphi]), np.array([1 + spsi, 0.0, -1 + spsi]))
    edge_12 = (np.array([1 + sphi, -1 + sphi, 0.0]), np.array([1 + spsi, -1 + spsi, 0.0]))
    edge_23 = (np.array([0.0, 1 + sphi, -1 + sphi]), np.array([0.0, 1 + spsi, -1 + spsi]))
    return main, edge_12, edge_23


def mohr_coulomb(deps, sigma, history, params):
    """params: E, nu, phi, c, psi (degrees, Pa).

    Return mapping in principal stress space (de Souza Neto et al., section 8.2.2):
    main plane, then the edge that the main-plane return crosses (triaxial compression
    or extension corner), then the apex. Because the surface is a set of planes and the
    material is perfectly plastic, each return is a small linear solve and is exact.
    """
    K, G = bulk_shear(params)
    c = params["c"]
    sphi = math.sin(math.radians(params["phi"]))
    cphi = math.cos(math.radians(params["phi"]))
    spsi = math.sin(math.radians(params["psi"]))
    trial = sigma + elastic_matrix(K, G) @ deps

    w, V = np.linalg.eigh(stress_to_tensor(trial))
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]          # w[0] >= w[1] >= w[2]
    k = 2.0 * c * cphi
    scale = max(np.max(np.abs(w)), c, 1.0)
    f_trial = (w[0] - w[2]) + (w[0] + w[2]) * sphi - k
    if f_trial <= YIELD_TOL * scale:
        return _elastic(trial, history)

    Dp = np.full((3, 3), K - 2.0 * G / 3.0) + np.eye(3) * 2.0 * G
    main, edge_12, edge_23 = _mc_principal_planes(sphi, spsi)

    def ret(planes):
        A = np.array([[nf @ Dp @ ng for (_, ng) in planes] for (nf, _) in planes])
        b = np.array([nf @ w - k for (nf, _) in planes])
        dg = np.linalg.solve(A, b)
        return w - Dp @ sum(d * ng for d, (_, ng) in zip(dg, planes)), dg

    tol = 1e-12 * scale
    apex = False
    s_new, dg = ret([main])
    if not (s_new[0] >= s_new[1] - tol and s_new[1] >= s_new[2] - tol):
        planes = [main, edge_23] if s_new[1] > s_new[0] else [main, edge_12]
        s_new, dg = ret(planes)
        ok = np.all(dg >= -tol) and s_new[0] >= s_new[1] - tol and s_new[1] >= s_new[2] - tol
        if not ok:
            if sphi == 0.0:
                raise RuntimeError("Tresca corner return failed")
            s_new = np.full(3, c * cphi / sphi)   # apex at c cot(phi)
            apex = True
    if not apex and sphi > 0.0 and (s_new[0] + s_new[2]) / 2.0 > c * cphi / sphi + tol:
        s_new = np.full(3, c * cphi / sphi)
        apex = True

    sig_new = tensor_to_stress(V @ np.diag(s_new) @ V.T)
    return _finish(trial, sig_new, history, K, G, apex=apex)


def mc_yield(sigma, params, history=None):
    w = np.linalg.eigvalsh(stress_to_tensor(sigma))
    phi = math.radians(params["phi"])
    return (w[-1] - w[0]) + (w[-1] + w[0]) * math.sin(phi) - 2.0 * params["c"] * math.cos(phi)


MODELS = {
    "linearElastic": linear_elastic,
    "vonMises": von_mises,
    "DruckerPrager": drucker_prager,
    "MohrCoulomb": mohr_coulomb,
}
