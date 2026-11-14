#!/usr/bin/env python3
"""WP 1.2: 1D elastic-plastic Lagrange-plus-remap prototype.

Reuses the WP 1.1 bar (``lr1d.Bar1D``: cell-centred MUSCL-Hancock Godunov
Lagrange step, conservative MUSCL remap, MC limiter) and adds:

- Uniaxial strain. The bar is laterally confined, so the strain increment of a
  cell is (d eps_xx, 0, 0, 0, 0, 0). Each cell keeps the full 6-component
  stress (Voigt order xx, yy, zz, xy, yz, zx, decision 0004) and its history
  (``eps_p`` 6 components, ``eps_p_eq``). A WP 1.4 model
  (``prototype/constitutive/models.py``) maps the increment to the new stress.
  Only sigma_xx enters the momentum equation and the Riemann solver.
- Wave speeds for the Riemann solver, the characteristic slopes and the time
  step are the elastic ones, c = sqrt((K + 4G/3)/rho). A plastic wave is
  slower, so the elastic impedance is an upper bound and the face solve stays
  dissipative (decision 0008).
- Remap of all 6 stress components and all history variables, mass-weighted:
  each is reconstructed in mass coordinates and carried by the mass flux, the
  same way WP 1.1 carries velocity. So sum(m * q) is conserved exactly and q
  stays inside the range of its neighbours. Stress uses the MC limiter like
  the other fields. History uses the more compressive superbee limiter: with
  MC a plastic-strain edge kept widening (5 cells 10-90 % after 60 cells of
  travel, 6 after 120); with superbee it stays 2 cells (decision 0008).

SI units, stress positive in tension. Periodic bar, area 1 m^2.

Run ``python3 lr1d_plastic.py --verify`` for the WP 1.2 acceptance checks
(used by tests/test_proto_1d_plastic.sh), ``--dp`` for the Drucker-Prager run.
"""

import argparse
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "constitutive"))

from lr1d import Bar1D, front_position, left, limiter, right  # noqa: E402
import models  # noqa: E402


def superbee(dl, dr):
    """Superbee limiter (Roe 1985): max(min(2|dl|, |dr|), min(|dl|, 2|dr|))
    with the common sign, zero at extrema. TVD and the most compressive of the
    standard limiters, so a discontinuity keeps a fixed width of about two
    cells instead of slowly widening. It squares off smooth peaks, so it is
    used only for history variables (decision 0008)."""
    a, b = np.abs(dl), np.abs(dr)
    mag = np.maximum(np.minimum(2.0 * a, b), np.minimum(a, 2.0 * b))
    return np.where(dl * dr > 0.0, np.sign(dl) * mag, 0.0)


def mass_weighted_flux(q, m, mass_flux, lim=limiter):
    """Flux of m*q through each face, for a remap that moves mass_flux
    (positive left to right, face j between cell j-1 and cell j).

    q is reconstructed linearly in mass coordinates with a limited slope and
    averaged over the swept mass, as lr1d.Bar1D.remap does for velocity.
    q may be (n,) or (n, k); the columns are limited independently.
    lim is the slope limiter (MC by default).
    """
    if q.ndim == 2:
        return np.stack([mass_weighted_flux(q[:, k], m, mass_flux, lim)
                         for k in range(q.shape[1])], axis=1)
    slope = lim((q - left(q)) / (0.5 * (m + left(m))),
                    (right(q) - q) / (0.5 * (m + right(m))))
    am = np.abs(mass_flux)
    q_from_left = left(q) + left(slope) * 0.5 * (left(m) - am)
    q_from_right = q - slope * 0.5 * (m - am)
    return mass_flux * np.where(mass_flux > 0.0, q_from_left, q_from_right)


class PlasticBar1D(Bar1D):
    """Periodic bar in uniaxial strain with a WP 1.4 constitutive model."""

    def __init__(self, n, length, rho, model, params):
        K, G = models.bulk_shear(params)
        super().__init__(n, length, rho, K + 4.0 * G / 3.0)   # elastic 1D modulus
        self.K, self.G = K, G
        self.model = models.MODELS[model] if isinstance(model, str) else model
        self.params = params
        self.sig6 = np.zeros((n, 6))
        self.eps_p = np.zeros((n, 6))
        self.eps_p_eq = np.zeros(n)
        # Mass-weighted eps_p_eq made by the constitutive model (not by the
        # remap), summed over the run. The streaming check needs it to be zero.
        self.plastic_made = 0.0
        self._lag = None

    def set_stress(self, sig6):
        self.sig6 = np.array(sig6, dtype=float)
        self.sig = self.sig6[:, 0].copy()

    # --- Lagrange step: stress from the constitutive model --------------------

    def update_stress(self, h, h_new):
        deps = np.zeros(6)
        de = (h_new - h) / (0.5 * (h + h_new))
        sig6 = np.empty_like(self.sig6)
        eps_p = np.empty_like(self.eps_p)
        eps_p_eq = np.empty_like(self.eps_p_eq)
        for j in range(self.n):
            deps[0] = de[j]
            hist = {"eps_p": self.eps_p[j], "eps_p_eq": self.eps_p_eq[j], "apex": False}
            sig6[j], hn = self.model(deps, self.sig6[j], hist, self.params)
            eps_p[j] = hn["eps_p"]
            eps_p_eq[j] = hn["eps_p_eq"]
        self.plastic_made += (self.rho * h * (eps_p_eq - self.eps_p_eq)).sum()
        self._lag = (sig6, eps_p, eps_p_eq)
        return sig6[:, 0].copy()

    # --- remap ------------------------------------------------------------------

    def remap(self, d, h, rho, sig, v):
        sig6, eps_p, eps_p_eq = self._lag
        m = rho * h
        mass_flux = d * self.swept_average(rho, h, d)
        m_new = m + mass_flux - right(mass_flux)

        def carry(q, lim=limiter):
            f = mass_weighted_flux(q, m, mass_flux, lim)
            mq = (m * q.T).T if q.ndim == 2 else m * q
            return ((mq + f - np.roll(f, -1, axis=0)).T / m_new).T

        self.v = carry(v)
        self.sig6 = carry(sig6)
        self.eps_p = carry(eps_p, superbee)
        self.eps_p_eq = carry(eps_p_eq, superbee)
        self.sig = self.sig6[:, 0].copy()
        self.rho = m_new / self.dx


# ----------------------------------------------------------------------------
# Materials and the impact problem
# ----------------------------------------------------------------------------

RHO = 2000.0                    # kg/m^3
VM = dict(E=20.0e6, nu=0.3, sigma_y=10.0e3, H=0.0)       # von Mises, no hardening
# Drucker-Prager for the --dp run. With nu 0.3 and phi 30 the elastic
# uniaxial-strain path (K0 = nu/(1-nu) = 0.43) stays inside the cone and never
# yields in compression, so the run uses nu 0.2 (K0 = 0.25) and phi 20.
DP = dict(E=20.0e6, nu=0.2, phi=20.0, c=5.0e3, psi=0.0, fit="compression")
LENGTH = 1.0


def impact_bar(n, params, model, V, v0):
    """Periodic bar: left half moves at v0 + V, right half at v0 - V.

    At x = L/2 the halves collide: a compression structure runs each way.
    At x = 0 (= L) they separate: a tension structure runs each way. Both
    stay apart for t < L/(4 c_e).
    """
    bar = PlasticBar1D(n, LENGTH, RHO, model, params)
    bar.v = np.where(bar.xc < 0.5 * LENGTH, v0 + V, v0 - V)
    return bar


def largest_reversal(q, increasing):
    """Largest amount by which q goes back against its expected direction
    (0 for a monotone sequence)."""
    q = np.asarray(q) if increasing else -np.asarray(q)
    return float(np.max(np.maximum.accumulate(q) - q))


# ----------------------------------------------------------------------------
# Acceptance checks for WP 1.2
# ----------------------------------------------------------------------------

def uniaxial_response(model, params, rho):
    """Probe a model in uniaxial-strain compression from zero stress.

    Returns (hel, c_p): the axial stress magnitude at first yield (found by
    bisection on a single elastic step) and the plastic wave speed from the
    tangent d sigma_xx / d eps_xx over the next 2 yield strains. For von Mises
    these equal Y (K + 4G/3)/(2G) and sqrt(K/rho); for Drucker-Prager they are
    the reference values of the --dp run.
    """
    K, G = models.bulk_shear(params)
    M = K + 4.0 * G / 3.0
    h0 = models.init_history()
    z6 = np.zeros(6)

    def plastic(e):
        _, h = model(np.array([-e, 0, 0, 0, 0, 0.0]), z6, h0, params)
        return h["eps_p_eq"] > 0.0

    lo, hi = 0.0, 1e-6
    while not plastic(hi):
        lo, hi = hi, 2.0 * hi
        if hi > 1.0:
            raise ValueError("no yield in uniaxial-strain compression")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if plastic(mid) else (mid, hi)
    e_y = hi
    sig, h = model(np.array([-e_y, 0, 0, 0, 0, 0.0]), z6, h0, params)
    s0 = sig[0]
    for _ in range(200):
        sig, h = model(np.array([-e_y / 100.0, 0, 0, 0, 0, 0.0]), sig, h, params)
    tangent = (sig[0] - s0) / (-2.0 * e_y)
    return M * e_y, math.sqrt(tangent / rho)


def step_loading(params, model, n, v0_frac, hel, c_p, startup_cells=10.0):
    """Step loading above the HEL by a symmetric impact (impact_bar).

    hel and c_p are the reference precursor amplitude and plastic wave speed.
    The impact velocity is 3 times the precursor particle-velocity jump, so
    the final stress is about 2.6 HEL (von Mises). Extrema are tracked over
    the whole run but reported in two windows: the start-up (until the
    elastic wave has run startup_cells cells) and the rest.
    Returns a dict of measured values.
    """
    K, G = models.bulk_shear(params)
    c_e = math.sqrt((K + 4.0 * G / 3.0) / RHO)
    v0 = v0_frac * c_e
    v1 = hel / (RHO * c_e)
    V = 3.0 * v1
    sig_final = -(hel + RHO * c_p * (V - v1))

    bar = impact_bar(n, params, model, V, v0)
    m0, p0 = bar.total_mass(), bar.total_momentum()
    pscale = (bar.cell_mass() * np.abs(bar.v)).sum()
    keys = ("sig", "v", "bound_sig", "bound_v")
    worst = {"start": dict.fromkeys(keys, 0.0), "main": dict.fromkeys(keys, 0.0)}

    def region(br):
        """Cells from the impact plane to well ahead of the right-going
        precursor, ordered left to right (the impact plane moves with v0)."""
        x_imp = (0.5 * LENGTH + v0 * br.t) % LENGTH
        reach = min(1.02 * c_e * br.t + 10 * br.dx, 0.24 * LENGTH)
        idx = np.arange(br.n)
        rel = (br.xc - x_imp) % LENGTH
        sel = idx[rel <= reach]
        return sel[np.argsort(rel[sel])]

    def track(br, k):
        w = worst["start" if c_e * br.t < startup_cells * br.dx else "main"]
        sel = region(br)
        w["sig"] = max(w["sig"], largest_reversal(br.sig[sel], True))
        w["v"] = max(w["v"], largest_reversal(br.v[sel], False))
        w["bound_sig"] = max(w["bound_sig"], np.abs(br.sig).max() - abs(sig_final))
        w["bound_v"] = max(w["bound_v"], np.abs(br.v - v0).max() - V)

    def fronts(t):
        x_imp = 0.5 * LENGTH + v0 * t
        xe = front_position(bar.xc, -bar.sig, 0.5 * hel, x_imp + (c_e - V) * t, LENGTH)
        xp = front_position(bar.xc, -bar.sig, hel + 0.5 * (-sig_final - hel),
                            x_imp + c_p * t, LENGTH)
        return xe, xp

    t1 = 0.08 * LENGTH / c_e
    t2 = 0.22 * LENGTH / c_e        # before the compression and tension waves meet
    bar.run(t1, callback=track)
    xe1, xp1 = fronts(t1)
    bar.run(t2 - t1, callback=track)
    xe2, xp2 = fronts(t2)

    # Lab speeds minus the velocity of the material just ahead of each front.
    u_e = (xe2 - xe1) / (t2 - t1) - (v0 - V)
    u_p = (xp2 - xp1) / (t2 - t1) - (v0 - V + v1)

    # Precursor plateau: middle third between the two fronts.
    rel = (bar.xc - xp2) % LENGTH
    gap = (xe2 - xp2) % LENGTH
    plateau = bar.sig[(rel > gap / 3.0) & (rel < 2.0 * gap / 3.0)]

    return dict(c_e=c_e, c_p=c_p, hel=hel, sig_final=sig_final, V=V, v0=v0,
                u_e=u_e, u_p=u_p, amp=-np.median(plateau),
                plateau_cells=plateau.size, worst=worst,
                dm=abs(bar.total_mass() - m0) / m0,
                dp=abs(bar.total_momentum() - p0) / pscale,
                finite=bool(np.all(np.isfinite(bar.sig6))),
                plastic_cells=int(np.sum(bar.eps_p_eq > 0.0)),
                apex_tension=float(bar.sig.max()), bar=bar)


def streaming(params, model, n, v0_frac, cells_moved):
    """A pre-strained patch (eps_p_eq 0.02, deviatoric eps_p) at zero stress,
    carried through the mesh at v0 with no loading. Returns a dict."""
    K, G = models.bulk_shear(params)
    c_e = math.sqrt((K + 4.0 * G / 3.0) / RHO)
    v0 = v0_frac * c_e
    bar = PlasticBar1D(n, LENGTH, RHO, model, params)
    bar.v[:] = v0
    a, b = 0.2 * LENGTH, 0.4 * LENGTH
    inside = (bar.xc >= a) & (bar.xc < b)
    epq = 0.02
    bar.eps_p_eq = np.where(inside, epq, 0.0)
    bar.eps_p = np.outer(bar.eps_p_eq, [-1.0, 0.5, 0.5, 0.0, 0.0, 0.0])
    m0 = bar.total_mass()
    tot0 = (bar.cell_mass() * bar.eps_p_eq).sum()
    tot0_p = (bar.cell_mass()[:, None] * bar.eps_p).sum(axis=0)
    bar.run(cells_moved * bar.dx / v0)

    tot = (bar.cell_mass() * bar.eps_p_eq).sum()
    tot_p = (bar.cell_mass()[:, None] * bar.eps_p).sum(axis=0)
    q = bar.eps_p_eq / epq
    # Signed distance from the expected patch centre (periodic)
    rel = bar.xc - (0.5 * (a + b) + v0 * bar.t)
    rel = (rel + 0.5 * LENGTH) % LENGTH - 0.5 * LENGTH

    def width(lo, hi):
        mid = (lo < q) & (q < hi)
        return int(np.sum(mid & (rel < 0.0))), int(np.sum(mid & (rel >= 0.0)))

    far = np.abs(rel) > 0.5 * (b - a) + 4.0 * bar.dx
    return dict(v0=v0, cells_moved=v0 * bar.t / bar.dx,
                drift=abs(tot - tot0) / tot0,
                drift_p=np.abs(tot_p - tot0_p).max() / tot0,
                made=abs(bar.plastic_made) / tot0,
                far_max=float(np.abs(q[far]).max()),
                neg=float(max(0.0, -bar.eps_p_eq.min())),
                w10=width(0.1, 0.9), w01=width(0.01, 0.99),
                dm=abs(bar.total_mass() - m0) / m0,
                sig_max=float(np.abs(bar.sig6).max()))


def verify(v0_frac=0.25, n=1000):
    ok = True

    def check(name, passed, detail):
        nonlocal ok
        print(f"  {name}: {detail}  {'PASS' if passed else 'FAIL'}")
        ok &= passed

    K, G = models.bulk_shear(VM)
    hel = VM["sigma_y"] * (K + 4.0 * G / 3.0) / (2.0 * G)
    c_p = math.sqrt(K / RHO)
    print(f"Material: rho {RHO} kg/m^3, E {VM['E']:.3e} Pa, nu {VM['nu']}, "
          f"von Mises Y {VM['sigma_y']:.3e} Pa, no hardening")
    print(f"Bar {LENGTH} m, {n} cells, material streaming at {v0_frac} c_e.")

    print("\n[1] Step loading above the HEL (von Mises): two waves")
    r = step_loading(VM, "vonMises", n, v0_frac, hel, c_p)
    print(f"  impact velocity +-{r['V']:.4f} m/s, exact final stress "
          f"{r['sig_final']:.4e} Pa, HEL {r['hel']:.4e} Pa")
    e_e = abs(r["u_e"] - r["c_e"]) / r["c_e"]
    e_p = abs(r["u_p"] - r["c_p"]) / r["c_p"]
    e_a = abs(r["amp"] - r["hel"]) / r["hel"]
    check("precursor speed", e_e <= 0.01,
          f"{r['u_e']:.3f} m/s vs sqrt((K+4G/3)/rho) {r['c_e']:.3f}, "
          f"error {100 * e_e:.3f} % (tol 1 %)")
    check("plastic wave speed", e_p <= 0.02,
          f"{r['u_p']:.3f} m/s vs sqrt(K/rho) {r['c_p']:.3f}, "
          f"error {100 * e_p:.3f} % (tol 2 %)")
    check("precursor amplitude", e_a <= 0.02,
          f"{r['amp']:.4e} Pa vs Y(K+4G/3)/(2G) {r['hel']:.4e} Pa, "
          f"error {100 * e_a:.3f} % ({r['plateau_cells']} plateau cells, tol 2 %)")

    # Extrema. "Reversal": the profile from the impact plane to ahead of the
    # precursor should rise monotonically in stress and fall in velocity;
    # this is the largest step back. "Beyond range": |stress| above the exact
    # final stress, |v - v0| above the impact velocity, anywhere in the bar.
    # Tolerance 1 % of the jump (final stress, impact velocity).
    js, jv = abs(r["sig_final"]), r["V"]
    w, s0 = r["worst"]["main"], r["worst"]["start"]
    print("\n[2] No new extrema in the step loading run, after the start-up")
    print("  (start-up = until the elastic wave has run 10 cells; decision 0008)")
    check("stress", max(w["sig"], w["bound_sig"]) <= 0.01 * js,
          f"reversal {100 * w['sig'] / js:.3f} %, beyond range "
          f"{100 * w['bound_sig'] / js:.3f} % of the jump (tol 1 %)")
    check("velocity", max(w["v"], w["bound_v"]) <= 0.01 * jv,
          f"reversal {100 * w['v'] / jv:.3f} %, beyond range "
          f"{100 * w['bound_v'] / jv:.3f} % of the jump (tol 1 %)")
    print(f"  start-up, for the record (known issue, not checked): stress reversal "
          f"{100 * s0['sig'] / js:.3f} %, beyond range {100 * s0['bound_sig'] / js:.3f} %; "
          f"velocity {100 * s0['v'] / jv:.3f} %, {100 * s0['bound_v'] / jv:.3f} %")

    print("\n[3] Conservation in the step loading run")
    check("mass", r["dm"] <= 1e-12, f"relative change {r['dm']:.3e}")
    check("momentum", r["dp"] <= 1e-12, f"relative change {r['dp']:.3e}")

    print("\n[4] Streaming: pre-strained patch (eps_p_eq 0.02, zero stress), no loading")
    widths = []
    for cells in (60, 120):
        s = streaming(VM, "vonMises", 400, v0_frac, cells)
        widths.append(max(s["w10"]))
        print(f"  patch carried {s['cells_moved']:.1f} cells at {s['v0']:.3f} m/s:")
        check("total eps_p_eq kept", s["drift"] <= 1e-12 and s["drift_p"] <= 1e-12,
              f"relative change {s['drift']:.3e}; eps_p components {s['drift_p']:.3e}")
        # Stress stays at round-off: the strain increments are rounding noise
        # of the cell lengths (about 1e-16), times E.
        check("no plastic strain made", s["made"] == 0.0 and s["sig_max"] <= 1e-8 * VM["sigma_y"],
              f"model-made eps_p_eq {s['made']:.1e}, largest stress {s['sig_max']:.1e} Pa")
        check("none far from the patch", s["far_max"] <= 0.01 and s["neg"] == 0.0,
              f"largest eps_p_eq more than 4 cells outside the edges "
              f"{100 * s['far_max']:.3f} % of the patch value, most negative {s['neg']:.1e}")
        check("edges within 4 cells", max(s["w10"]) <= 4,
              f"cells between 10 % and 90 %: {s['w10'][0]} and {s['w10'][1]} "
              f"(between 1 % and 99 %: {s['w01'][0]} and {s['w01'][1]})")
        check("mass", s["dm"] <= 1e-12, f"relative change {s['dm']:.3e}")
    check("edges stop widening", widths[1] <= widths[0],
          f"10-90 % width {widths[0]} cells after 60, {widths[1]} after 120")

    print("\nWP 1.2 checks:", "ALL PASS" if ok else "FAILED")
    return ok


def dp_run(v0_frac=0.25, n=1000):
    """Drucker-Prager step loading, recorded in the WP 1.2 log (not a pass
    criterion). The tension structure at x = 0 reaches the cone apex."""
    params = DP
    hel, c_p = uniaxial_response(models.drucker_prager, params, RHO)
    K, G = models.bulk_shear(params)
    print(f"Drucker-Prager {params}")
    print(f"Probe: HEL {hel:.4e} Pa, plastic wave speed {c_p:.3f} m/s "
          f"(sqrt(K/rho) {math.sqrt(K / RHO):.3f}, c_e {math.sqrt((K + 4 * G / 3) / RHO):.3f})")
    r = step_loading(params, "DruckerPrager", n, v0_frac, hel, c_p)
    js, jv = abs(r["sig_final"]), r["V"]
    w, s0 = r["worst"]["main"], r["worst"]["start"]
    print(f"impact +-{r['V']:.4f} m/s, reference final stress {r['sig_final']:.4e} Pa")
    print(f"precursor speed {r['u_e']:.3f} (c_e {r['c_e']:.3f}), plastic wave "
          f"{r['u_p']:.3f} (probe {r['c_p']:.3f}), precursor {r['amp']:.4e} (HEL {hel:.4e})")
    print(f"after start-up: stress reversal {100 * w['sig'] / js:.3f} %, beyond "
          f"{100 * w['bound_sig'] / js:.3f} %; velocity {100 * w['v'] / jv:.3f} %, "
          f"{100 * w['bound_v'] / jv:.3f} %")
    print(f"start-up: stress {100 * max(s0['sig'], s0['bound_sig']) / js:.3f} %, "
          f"velocity {100 * max(s0['v'], s0['bound_v']) / jv:.3f} %")
    print(f"mass {r['dm']:.1e}, momentum {r['dp']:.1e}, finite {r['finite']}, "
          f"plastic cells {r['plastic_cells']}, largest tension {r['apex_tension']:.4e} Pa")
    return 0 if r["finite"] else 1


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--verify", action="store_true", help="run the WP 1.2 checks")
    p.add_argument("--dp", action="store_true", help="Drucker-Prager step loading run")
    p.add_argument("--v0", type=float, default=0.25,
                   help="background velocity as a fraction of c_e (default 0.25)")
    p.add_argument("--cells", type=int, default=1000)
    args = p.parse_args(argv)
    if args.dp:
        return dp_run(args.v0, args.cells)
    return 0 if verify(args.v0, args.cells) else 1


if __name__ == "__main__":
    sys.exit(main())
