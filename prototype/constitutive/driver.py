"""Single-point driver for the WP 1.4 constitutive models.

Loads one material point along strain-driven or mixed strain/stress-controlled paths:
uniaxial strain, triaxial compression and extension, simple shear. Each step prescribes
some strain components and some stress components (Voigt order xx, yy, zz, xy, yz, zx);
the unknown strain components are found by Newton iteration with a finite-difference
Jacobian, so any model in models.MODELS can be driven without a tangent operator.

Usage (prints a short triaxial test):  python3 driver.py
"""

import numpy as np

import models


def mixed_step(model, params, sigma, history, deps_given, strain_ctrl, sig_target,
               tol=1e-9, max_iter=50):
    """One increment. strain_ctrl[i] True: deps[i] = deps_given[i]; False: the strain
    is free and sigma_new[i] must equal sig_target[i]."""
    deps = np.array(deps_given, dtype=float)
    free = ~np.asarray(strain_ctrl)
    if not free.any():
        return (*model(deps, sigma, history, params), deps)

    # Elastic predictor for the free components
    K, G = models.bulk_shear(params)
    D = models.elastic_matrix(K, G)
    ctrl = ~free
    rhs = sig_target[free] - sigma[free] - D[np.ix_(free, ctrl)] @ deps[ctrl]
    deps[free] = np.linalg.solve(D[np.ix_(free, free)], rhs)

    scale = max(np.max(np.abs(sig_target[free])), np.max(np.abs(sigma)), 1.0)
    for _ in range(max_iter):
        s, h = model(deps, sigma, history, params)
        r = s[free] - sig_target[free]
        if np.max(np.abs(r)) <= tol * scale:
            return s, h, deps
        J = np.zeros((free.sum(), free.sum()))
        step = 1e-9
        for j, idx in enumerate(np.flatnonzero(free)):
            d2 = deps.copy()
            d2[idx] += step
            J[:, j] = (model(d2, sigma, history, params)[0][free] - s[free]) / step
        deps[free] -= np.linalg.lstsq(J, r, rcond=None)[0]
    raise RuntimeError("mixed control did not converge, residual %.3g Pa" % np.max(np.abs(r)))


def run_path(model, params, sigma0, n_steps, deps_step, strain_ctrl, sig_target):
    """Apply n_steps equal increments. Returns dict with arrays eps (total strain) and
    sigma, shape (n_steps + 1, 6), and the list of history dicts."""
    sigma = np.array(sigma0, dtype=float)
    history = models.init_history()
    eps = np.zeros(6)
    out_eps, out_sig, out_h = [eps.copy()], [sigma.copy()], [history]
    for _ in range(n_steps):
        sigma, history, deps = mixed_step(model, params, sigma, history, deps_step,
                                          np.asarray(strain_ctrl), np.asarray(sig_target, float))
        eps = eps + deps
        out_eps.append(eps.copy())
        out_sig.append(sigma.copy())
        out_h.append(history)
    return {"eps": np.array(out_eps), "sigma": np.array(out_sig), "history": out_h}


def uniaxial_strain(model, params, eps_zz, n_steps, sigma0=None):
    """Oedometric loading: eps_zz prescribed, all other strains zero."""
    s0 = np.zeros(6) if sigma0 is None else sigma0
    d = np.zeros(6)
    d[2] = eps_zz / n_steps
    return run_path(model, params, s0, n_steps, d, [True] * 6, np.zeros(6))


def triaxial(model, params, p_conf, eps_axial, n_steps, extension=False):
    """Conventional triaxial test from an isotropic state -p_conf (p_conf > 0 is
    compression). Radial stresses stay at -p_conf; axial strain eps_axial (zz) is applied
    in n_steps. Compression: eps_axial < 0. Extension: eps_axial > 0."""
    if extension != (eps_axial > 0):
        raise ValueError("extension needs eps_axial > 0, compression eps_axial < 0")
    s0 = -p_conf * models.I_V
    d = np.zeros(6)
    d[2] = eps_axial / n_steps
    ctrl = [False, False, True, True, True, True]
    target = s0.copy()
    return run_path(model, params, s0, n_steps, d, ctrl, target)


def simple_shear(model, params, sigma_n, gamma, n_steps):
    """Simple shear in the x-y plane at constant normal stress sigma_yy = -sigma_n,
    from an isotropic state. gamma_xy is applied; eps_xx = eps_zz = 0 (plane strain);
    eps_yy (dilatancy) is free."""
    s0 = -sigma_n * models.I_V
    d = np.zeros(6)
    d[3] = gamma / n_steps
    ctrl = [True, False, True, True, True, True]
    return run_path(model, params, s0, n_steps, d, ctrl, s0.copy())


if __name__ == "__main__":
    mc = dict(E=20e6, nu=0.3, phi=30.0, c=10e3, psi=10.0)
    r = triaxial(models.mohr_coulomb, mc, 100e3, -0.05, 100)
    print("Mohr-Coulomb triaxial compression, p_conf = 100 kPa")
    print("  eps_zz        sigma_zz [Pa]   eps_v")
    for i in range(0, 101, 10):
        e = r["eps"][i]
        print("  %+.4f   %+.6e   %+.5f" % (e[2], r["sigma"][i, 2], e[0] + e[1] + e[2]))
