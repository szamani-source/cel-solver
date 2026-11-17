"""Rate-dependent friction mu(I) for Drucker-Prager and Mohr-Coulomb (WP 3.6a, decision 0037).

The friction law (Jop, Forterre and Pouliquen 2005, JFM 541, eq. 4.3):

    mu_J(I) = mu_s + (mu_2 - mu_s) / (I0 / I + 1),   I = gamma_dot d / sqrt(P / rho_s)

with the partial regularisation of Barker and Gray (2017, JFM 828, eqs. 5.2, 5.3) below
the first neutral point I_1^N of the well-posedness function (eq. 3.9)

    C(I) = 4 nu^2 - 4 nu + mu^2 (1 - nu / 2)^2,   nu = I mu'(I) / mu   (C > 0: ill posed)

    mu(I) = sqrt(alpha / ln(A / I))   for I <= I_1^N,   A = I_1^N exp(alpha / mu_J(I_1^N)^2)
    mu(I) = mu_J(I)                   for I >  I_1^N

alpha is not a material parameter: Barker and Gray require it close to and below 2 for
well-posedness and use 1.9 (their figure 5). The regularised mu goes to 0 as I -> 0, so
there is no yield stress: any shear stress creeps, at a rate that is negligible far below
mu_s (decision 0037). The high-I branch is Jop's (Barker and Gray's mu_inf term, eq. 6.3,
is left out).

How it enters the return mapping (decision 0037): mu sets a mobilised friction angle
phi_m = atan(mu). The Drucker-Prager or Mohr-Coulomb return is run with phi_m in place of
phi (same fit, same psi, c and sigma_t), so mu_s = tan(phi). The plastic shear strain
increment of the return, dgamma = sqrt(2 e:e) of the deviatoric plastic strain (tensor
components; equal to gamma_xy in simple shear), and the mean pressure after the return,
P = -tr(sigma)/3, give I = dgamma / dt d sqrt(rho_s / P). phi_m is the root of

    R(mu) = mu - mu_law(I(mu)) = 0.

R increases with mu (a larger cone gives less plastic flow, hence smaller I and smaller
mu_law), so the root is unique and bracketed; it is found by the Illinois (modified
regula falsi) method to 1e-14 mu_2. mu_2 = mu_s (or mu_2 absent) is the rate-independent
model, called unchanged (bit for bit).

Parameters (in addition to those of the base model, SI units):
    mu2    friction at large I [-]          (absent or equal to tan(phi): rate-independent)
    I0     inertial number scale [-]
    d      grain diameter [m]
    rho_s  grain (solid) density [kg/m^3]
    muI_alpha  Barker-Gray alpha [-], in (0, 2); absent or 0: no regularisation (plain
               Jop law: elastic below mu_s, as the present model). Barker and Gray use
               1.9. The regularisation exists only if the Jop law has a well-posed
               range of I; for mu_s = tan 35 deg it has none for any mu_2 - mu_s up to
               0.8 (decision 0037), so T0 runs the plain Jop law.
    dt     time step of the increment [s] (set by the caller for each increment)
"""

import math

import numpy as np

import models

# Barker and Gray (2017) figure 5
ALPHA_BARKER_GRAY = 1.9


class FrictionLaw:
    """mu(I): Jop law, optionally with the Barker-Gray low-I branch (alpha > 0)."""

    def __init__(self, mu_s, mu_2, I0, alpha=0.0):
        if not (mu_2 > mu_s > 0.0 and I0 > 0.0):
            raise ValueError("mu(I) needs mu_2 > mu_s > 0 and I0 > 0")
        if alpha < 0.0 or alpha >= 2.0:
            raise ValueError("Barker-Gray alpha must be in (0, 2) (0: no regularisation)")
        self.mu_s, self.mu_2, self.I0, self.alpha = mu_s, mu_2, I0, alpha
        self.I1N = None
        self.A = None
        if alpha > 0.0:
            self.I1N = self.first_neutral_point()
            self.A = self.I1N * math.exp(alpha / self.mu_jop(self.I1N) ** 2)

    def mu_jop(self, I):
        if math.isinf(I):
            return self.mu_2
        return self.mu_s + (self.mu_2 - self.mu_s) * I / (self.I0 + I)

    def well_posedness(self, I):
        """C(I) of Barker and Gray (2017) eq. 3.9 for the Jop law (C > 0: ill posed)."""
        mu = self.mu_jop(I)
        dmu = (self.mu_2 - self.mu_s) * self.I0 / (self.I0 + I) ** 2
        nu = I * dmu / mu
        return 4.0 * nu * nu - 4.0 * nu + mu * mu * (1.0 - 0.5 * nu) ** 2

    def first_neutral_point(self):
        """Smallest I > 0 with C(I) = 0: scan log10 I in [-12, 3] by 0.01, then bisect."""
        n = 1500
        prev = 10.0 ** -12
        if self.well_posedness(prev) <= 0.0:
            raise ValueError("mu(I) already well posed at I = 1e-12: no regularisation needed")
        for k in range(1, n + 1):
            I = 10.0 ** (-12 + 15.0 * k / n)
            if self.well_posedness(I) <= 0.0:
                lo, hi = math.log(prev), math.log(I)
                for _ in range(200):
                    mid = 0.5 * (lo + hi)
                    if self.well_posedness(math.exp(mid)) > 0.0:
                        lo = mid
                    else:
                        hi = mid
                    if hi - lo <= 1e-15:
                        break
                return math.exp(hi)
            prev = I
        raise ValueError("mu(I) parameters give no well-posed range of I (Barker-Gray "
                         "regularisation impossible); use muI_alpha 0")

    def mu(self, I):
        if I <= 0.0:
            return 0.0 if self.alpha > 0.0 else self.mu_s
        if self.alpha > 0.0 and I <= self.I1N:
            return math.sqrt(self.alpha / math.log(self.A / I))
        return self.mu_jop(I)


def plastic_shear_increment(sig_trial, sig_new, K, G):
    """sqrt(2 e:e) of the deviatoric plastic strain C (trial - new), tensor components."""
    deps_p = models.compliance_matrix(K, G) @ (sig_trial - sig_new)
    e = deps_p.copy()
    e[3:] *= 0.5
    dev = e - (e[0] + e[1] + e[2]) / 3.0 * models.I_V
    return math.sqrt(2.0 * (dev[:3] @ dev[:3] + 2.0 * (dev[3:] @ dev[3:])))


def friction_law(params):
    """FrictionLaw from params, or None if rate-independent (mu2 absent or = tan phi)."""
    if "mu2" not in params:
        return None
    mu_s = math.tan(math.radians(params["phi"]))
    if params["mu2"] == mu_s:
        return None
    if params["mu2"] < mu_s:
        raise ValueError("mu2 = %g must not be below mu_s = tan(phi) = %g"
                         % (params["mu2"], mu_s))
    return FrictionLaw(mu_s, params["mu2"], params["I0"],
                       params.get("muI_alpha", 0.0))


def inertial_number(dgam, sig_new, params):
    if dgam == 0.0:
        return 0.0
    P = -models.mean_stress(sig_new)
    if not P > 1e-300 * params["rho_s"]:   # as celMuI.C (overflow guard, WP 3.6s)
        return math.inf
    return dgam / params["dt"] * params["d"] * math.sqrt(params["rho_s"] / P)


def rate_dependent(base):
    """Wrap a base model (drucker_prager or mohr_coulomb) with mu(I) friction."""

    def model(deps, sigma, history, params):
        law = friction_law(params)
        if law is None:
            return base(deps, sigma, history, params)
        K, G = models.bulk_shear(params)
        trial = sigma + models.elastic_matrix(K, G) @ deps

        def at(mu):
            """Base return at friction mu: (R, sigma, history, I)."""
            if mu == law.mu_s:
                p = params                      # the rate-independent model, unchanged
            else:
                p = dict(params, phi=math.degrees(math.atan(mu)))
                if p.get("sigma_t") is not None:
                    # A larger phi_m lowers the apex c cot(phi_m): keep the cut-off
                    # at or below it (decision 0037)
                    ph = math.radians(p["phi"])
                    p["sigma_t"] = min(p["sigma_t"], p["c"] * math.cos(ph) / math.sin(ph))
            s, h = base(deps, sigma, history, p)
            I = inertial_number(plastic_shear_increment(trial, s, K, G), s, params)
            return mu - law.mu(I), s, h, I

        ra, sa, ha, Ia = at(law.mu_s)
        if law.alpha == 0.0 and Ia == 0.0:
            return sa, ha                       # elastic below mu_s (Jop law)
        if ra == 0.0:
            return sa, ha
        if ra < 0.0:
            a, b = law.mu_s, law.mu_2
            rb, sb, hb, _ = at(b)
        else:
            a, b = 1e-3 * law.mu_s, law.mu_s
            rb, sb, hb = ra, sa, ha
            ra, sa, ha, _ = at(a)
            if ra >= 0.0:
                return sa, ha                   # creep below double precision
        if rb <= 0.0:
            return sb, hb                       # root at mu_2 (I infinite, e.g. P = 0)
        return _illinois(at, a, b, ra, rb, sa, ha, sb, hb, 1e-14 * law.mu_2)

    model.__name__ = "mu_i_" + base.__name__
    return model


def _illinois(at, a, b, fa, fb, sa, ha, sb, hb, tol, max_iter=200):
    """Illinois root of R on [a, b] with fa < 0 < fb; returns (sigma, history)."""
    side = 0
    for _ in range(max_iter):
        c = (a * fb - b * fa) / (fb - fa)
        if not (a < c < b):
            c = 0.5 * (a + b)
        fc, sc, hc, _ = at(c)
        if fc == 0.0:
            return sc, hc
        if fc < 0.0:
            a, fa, sa, ha = c, fc, sc, hc
            if side == -1:
                fb *= 0.5
            side = -1
        else:
            b, fb, sb, hb = c, fc, sc, hc
            if side == 1:
                fa *= 0.5
            side = 1
        if b - a <= tol:
            break
    else:
        raise RuntimeError("mu(I) return did not converge")
    # Take the end with the smaller residual
    return (sa, ha) if abs(fa) < abs(fb) else (sb, hb)


mu_i_drucker_prager = rate_dependent(models.drucker_prager)
mu_i_mohr_coulomb = rate_dependent(models.mohr_coulomb)


def sand_placeholder():
    """Placeholder T0 sand (decision 0037): mu_s = tan 35 deg (T0), mu_2 - mu_s and I0
    from Jop, Forterre and Pouliquen (2005) glass beads (tan 32.76 - tan 20.9, 0.279),
    d = 1 mm and rho_s = 2650 kg/m^3 (quartz). Not measured on Sia's material."""
    mu_s = math.tan(math.radians(35.0))
    dmu = math.tan(math.radians(32.76)) - math.tan(math.radians(20.9))
    return dict(phi=35.0, mu2=mu_s + dmu, I0=0.279, d=1e-3, rho_s=2650.0)


if __name__ == "__main__":
    p = sand_placeholder()
    law = friction_law(dict(p, muI_alpha=0.0))
    print("placeholder sand (plain Jop law): mu_s %.6f mu_2 %.6f I0 %.3f"
          % (law.mu_s, law.mu_2, law.I0))
    for I in (1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1.0):
        print("  I %.1e  mu %.6f  mu_J %.6f  C %+.3e" % (I, law.mu(I), law.mu_jop(I),
                                                       law.well_posedness(I)))
