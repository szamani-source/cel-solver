"""Nonlocal granular fluidity (NGF), WP N.1, decision 0111.

Kamrin and Koval (2012), Henann and Kamrin (2013, PNAS 110, 6730), in the
elasto-viscoplastic form of Henann and Kamrin (2014, IJP 60, 145):

    mu = |tau| / P,  gamma_dot_p = g mu                       (flow rule)
    g_loc(mu, P) = sqrt(P / (rho_s d^2)) (mu - mu_s) / (b mu)  for mu > mu_s, else 0
    g - xi^2 lap(g) = g_loc,  xi = A d / sqrt(|mu - mu_s|)     (nonlocal equation)

The steady homogeneous limit is the linear local law mu = mu_s + b I. Here
b = (mu_2 - mu_s) / I0 from the solver's mu(I) entries (the small-I slope of
the Jop law), A from the entry nonlocalA.

The nonlocal equation is written without dividing by |mu - mu_s|:

    k g - lap(g) = k g_loc,  k = |mu - mu_s| / (A d)^2

so mu = mu_s needs no cap (k = 0 there: pure diffusion), with a floor k_min
so a layer at mu = mu_s everywhere stays solvable.

Point update with g fixed over the step (elastic predictor, then the
deviatoric stress relaxed along its own direction at fixed P):

    tau = tau_trial - G dt gamma_dot_p = tau_trial - G dt g tau / P
    =>  tau = tau_trial / (1 + G dt g / P)

(The C++ kernel, celNonlocal.C, takes the local part of g at the new
stress, implicit, and lags only g - g_loc(old): the explicit form here is
unstable once G dt mu g'(mu)/P > 2. Both have the same steady state.)

tau here is the equivalent shear stress sqrt(s:s / 2) and gamma_dot_p the
equivalent plastic shear rate sqrt(2 e_p:e_p) (both equal to the xy
components in simple shear).
"""

import math

import numpy as np

# Henann and Kamrin (2013) glass beads (to check against the paper, N.1)
GLASS = dict(muS=0.3819, mu2=0.6435, I0=0.278, A=0.48, d=1.0e-3, rhoS=2450.0)


class NGF:
    def __init__(self, muS, mu2, I0, A, d, rhoS, kmin_rel=1e-8):
        self.muS, self.A, self.d, self.rhoS = muS, A, d, rhoS
        self.b = (mu2 - muS) / I0
        self.kmin = kmin_rel / (A * d) ** 2

    def g_loc(self, mu, P):
        mu, P = np.asarray(mu, float), np.asarray(P, float)
        out = np.zeros(np.broadcast(mu, P).shape)
        m = (mu > self.muS) & (P > 0)
        mm, pp = np.broadcast_to(mu, out.shape)[m], np.broadcast_to(P, out.shape)[m]
        out[m] = np.sqrt(pp / (self.rhoS * self.d ** 2)) * (mm - self.muS) / (self.b * mm)
        return out

    def k(self, mu):
        return np.maximum(np.abs(np.asarray(mu, float) - self.muS) / (self.A * self.d) ** 2,
                          self.kmin)

    def local_mu(self, I):
        return self.muS + self.b * I

    def point_update(self, tau_trial, P, G, dt, g):
        return tau_trial / (1.0 + G * dt * g / P)


def solve_layer(ngf, y, mu, P, bc="neumann"):
    """Steady g across a 1D layer (cell centres y, uniform spacing) for the
    given mu(y), P(y): k g - g'' = k g_loc, zero gradient (or g = 0) at both
    ends. Second-order finite volumes, dense tridiagonal solve."""
    n = len(y)
    h = y[1] - y[0]
    k = ngf.k(mu)
    rhs = k * ngf.g_loc(mu, P)
    M = np.zeros((n, n))
    for i in range(n):
        M[i, i] = k[i]
        for j in (i - 1, i + 1):
            if 0 <= j < n:
                M[i, i] += 1 / h ** 2
                M[i, j] -= 1 / h ** 2
            elif bc == "dirichlet":
                M[i, i] += 2 / h ** 2
    return np.linalg.solve(M, rhs)


def homogeneous_check(ngf, P=1000.0, G=1.0e6, rate=1.0, dt=1e-6, steps=200000):
    """Simple shear at fixed P and shear rate, g = g_loc (uniform, so the
    nonlocal term vanishes): returns mu at the end and mu_s + b I."""
    tau = 0.0
    for _ in range(steps):
        mu = tau / P
        g = float(ngf.g_loc(mu, P))
        tau = ngf.point_update(tau + G * rate * dt, P, G, dt, g)
    I = rate * ngf.d / math.sqrt(P / ngf.rhoS)
    return tau / P, ngf.local_mu(I)


if __name__ == "__main__":
    n = NGF(**GLASS)
    for rate in (1.0, 10.0, 50.0):
        print(rate, homogeneous_check(n, rate=rate))
