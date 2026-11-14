#!/usr/bin/env python3
"""WP 1.1: 1D elastic Lagrange-plus-remap prototype.

A periodic 1D bar on a fixed uniform (Eulerian) mesh. Each cycle is

1. Lagrange step: the mesh moves with the material. Cell-centred
   second-order Godunov step (MUSCL-Hancock): limited linear reconstruction
   in characteristic variables, a half-step predictor inside each cell, then
   an acoustic Riemann solver at each face gives the face velocity v* and
   face stress s*. Faces move with v*; cell momentum changes by the s* jump.
2. Remap: the faces go back to their fixed positions. Mass, stress (volume
   weighted) and momentum are moved between cells with second-order monotone
   MUSCL fluxes (MC limiter), integrated exactly over the swept region.

All fields (density, velocity, stress) live at cell centres, as in the
OpenFOAM solver this prototype prepares for. Decision 0004 records why the
first, staggered-grid version was replaced. Cell j spans face j (left) to
face j+1 (right). Area is 1 m^2.

SI units. Stress positive in tension. M is the 1D (constrained) modulus.

Run ``python3 lr1d.py --help`` for the command line, or ``--verify`` for the
WP 1.1 acceptance checks (used by tests/test_proto_1d_elastic.sh).
"""

import argparse
import sys

import numpy as np


# ----------------------------------------------------------------------------
# Small helpers (all arrays periodic)
# ----------------------------------------------------------------------------

def right(a):
    """Value of the right neighbour: right(a)[j] = a[j+1]."""
    return np.roll(a, -1)


def left(a):
    """Value of the left neighbour: left(a)[j] = a[j-1]."""
    return np.roll(a, 1)


def limiter(dl, dr):
    """Monotonized central (MC) limiter of two one-sided differences
    (van Leer 1977): min(2|dl|, 2|dr|, |dl + dr|/2) with the common sign,
    zero at extrema. A linear reconstruction with this slope stays inside the
    range of the neighbouring cell values, so the schemes below are TVD.
    MC clips smooth peaks less than the harmonic van Leer limiter, which kept
    the smooth-wave order under 1.8 (decision 0004).
    """
    mag = np.minimum(np.minimum(2.0 * np.abs(dl), 2.0 * np.abs(dr)),
                     0.5 * np.abs(dl + dr))
    return np.where(dl * dr > 0.0, np.sign(dl) * mag, 0.0)


def limited_slope(q, h):
    """Limited slope (per metre) of cell values q on cells of length h."""
    dl = (q - left(q)) / (0.5 * (h + left(h)))
    dr = (right(q) - q) / (0.5 * (h + right(h)))
    return limiter(dl, dr)


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------

class Bar1D:
    """Periodic elastic bar, cell-centred Lagrange step plus MUSCL remap."""

    def __init__(self, n, length, rho, M):
        self.n = n
        self.length = length
        self.dx = length / n
        self.M = M
        self.rho0 = rho
        self.c = np.sqrt(M / rho)

        # Fixed Eulerian mesh
        self.xf = self.dx * np.arange(n)          # face positions (left face of cell j)
        self.xc = self.xf + 0.5 * self.dx         # cell centres

        # State on the Eulerian mesh, all at cell centres
        self.rho = np.full(n, rho)                # density
        self.sig = np.zeros(n)                    # stress
        self.v = np.zeros(n)                      # velocity
        self.t = 0.0

    # --- integral quantities ------------------------------------------------

    def cell_mass(self):
        return self.rho * self.dx

    def total_mass(self):
        return self.cell_mass().sum()

    def total_momentum(self):
        return (self.cell_mass() * self.v).sum()

    # --- Lagrange step --------------------------------------------------------

    def lagrange_step(self, dt):
        """Move the mesh with the material for one step.

        Returns face displacements d, Lagrangian cell lengths h and the
        Lagrangian cell density, stress and velocity at the new time.
        """
        h = np.full(self.n, self.dx)
        rho, sig, v = self.rho, self.sig, self.v
        m = rho * h
        z = np.sqrt(self.M * rho)                 # acoustic impedance rho*c

        # Limited slopes in characteristic variables. For rho v_t = s_x,
        # s_t = M v_x, dwR = ds - z dv is carried by right-going waves and
        # dwL = ds + z dv by left-going ones. Differences use the cell's own
        # z, so a uniform background velocity does not leak in (Galilean
        # invariance).
        gs_l = (sig - left(sig)) / (0.5 * (h + left(h)))
        gs_r = (right(sig) - sig) / (0.5 * (h + right(h)))
        gv_l = (v - left(v)) / (0.5 * (h + left(h)))
        gv_r = (right(v) - v) / (0.5 * (h + right(h)))
        wr = limiter(gs_l - z * gv_l, gs_r - z * gv_r)
        wl = limiter(gs_l + z * gv_l, gs_r + z * gv_r)
        ds = 0.5 * (wr + wl) * h                  # stress change across the cell
        dv = 0.5 * (wl - wr) / z * h              # velocity change across the cell

        # Hancock predictor: advance the cell-centre values by dt/2.
        sig_p = sig + 0.5 * dt * self.M * dv / h
        v_p = v + 0.5 * dt * ds / m

        # Face states: face j has cell j-1 on its left and cell j on its right.
        sL = left(sig_p + 0.5 * ds)
        vL = left(v_p + 0.5 * dv)
        zL = left(z)
        sR = sig_p - 0.5 * ds
        vR = v_p - 0.5 * dv
        zR = z

        # Acoustic Riemann solver (stress positive in tension).
        zs = zL + zR
        v_star = (zL * vL + zR * vR + (sR - sL)) / zs
        s_star = (zR * sL + zL * sR + zL * zR * (vR - vL)) / zs

        # Update. Momentum change telescopes, so it is conserved exactly.
        v_new = v + dt * (right(s_star) - s_star) / m
        d = dt * v_star                           # face displacement
        h_new = h + right(d) - d
        h_mid = 0.5 * (h + h_new)
        sig_new = sig + self.M * (h_new - h) / h_mid   # hypoelastic
        rho_new = m / h_new
        return d, h_new, rho_new, sig_new, v_new

    # --- remap ---------------------------------------------------------------

    @staticmethod
    def swept_average(q, h, d):
        """Average of the limited linear reconstruction of q over the region
        each face sweeps when it moves back from X_j + d_j to X_j.

        Face j sits between cell j-1 (left) and cell j (right). If d_j > 0 the
        region [X_j, X_j + d_j] is the right end of Lagrangian cell j-1 and
        goes to Eulerian cell j; if d_j < 0 it is the left end of Lagrangian
        cell j and goes to Eulerian cell j-1.
        """
        s = limited_slope(q, h)
        a = np.abs(d)
        from_left = left(q) + left(s) * 0.5 * (left(h) - a)
        from_right = q - s * 0.5 * (h - a)
        return np.where(d > 0.0, from_left, from_right)

    def remap(self, d, h, rho, sig, v):
        """Return Lagrangian fields to the fixed mesh (conservative)."""
        m = rho * h

        # Fluxes through face j, positive left-to-right (cell j-1 -> cell j).
        mass_flux = d * self.swept_average(rho, h, d)
        sigvol_flux = d * self.swept_average(sig, h, d)

        # Momentum: mass flux times the limited velocity of the swept mass,
        # reconstructed in mass coordinates so v stays within its neighbours.
        vslope_m = limiter((v - left(v)) / (0.5 * (m + left(m))),
                            (right(v) - v) / (0.5 * (m + right(m))))
        am = np.abs(mass_flux)
        v_from_left = left(v) + left(vslope_m) * 0.5 * (left(m) - am)
        v_from_right = v - vslope_m * 0.5 * (m - am)
        mom_flux = mass_flux * np.where(d > 0.0, v_from_left, v_from_right)

        m_new = m + mass_flux - right(mass_flux)
        sigvol_new = sig * h + sigvol_flux - right(sigvol_flux)
        p_new = m * v + mom_flux - right(mom_flux)

        self.rho = m_new / self.dx
        self.sig = sigvol_new / self.dx
        self.v = p_new / m_new

    # --- driver ----------------------------------------------------------------

    def stable_dt(self, cfl):
        return cfl * self.dx / (self.c + np.abs(self.v).max())

    def step(self, dt):
        self.remap(*self.lagrange_step(dt))
        self.t += dt

    def run(self, t_end, cfl=0.5, callback=None):
        """Advance by t_end with a constant step chosen to land on it."""
        nsteps = int(np.ceil(t_end / self.stable_dt(cfl)))
        dt = t_end / nsteps
        for k in range(nsteps):
            self.step(dt)
            if callback is not None:
                callback(self, k + 1)
        return nsteps


# ----------------------------------------------------------------------------
# Test problems
# ----------------------------------------------------------------------------

# Default material: a dense sand-like solid (WP 1.1 only needs elasticity).
RHO = 2000.0        # kg/m^3
MODULUS = 50.0e6    # Pa, so c = 158.1 m/s
LENGTH = 1.0        # m, periodic


def gaussian_cell_average(xf, dx, x0, w, period):
    """Exact cell averages of the periodic Gaussian on cells [xf, xf+dx]."""
    from math import erf, sqrt, pi
    erfv = np.vectorize(erf)
    a = (xf - x0 + 0.5 * period) % period - 0.5 * period
    b = a + dx
    return 0.5 * sqrt(pi) * w * (erfv(b / w) - erfv(a / w)) / dx


def init_right_going(bar, sig_cells, amp, v0):
    """Right-going wave: sig = amp*shape, v = v0 - sig/(rho c).

    For rho v_t = s_x, s_t = M v_x a right-going wave f(x - ct) has
    v = -f/(rho c). v0 is a uniform background velocity, which makes the
    material stream through the fixed mesh so the remap does real work.
    """
    bar.sig = amp * sig_cells
    bar.v = v0 - bar.sig / (bar.rho0 * bar.c)


def smooth_case(n, v0, t_end, amp=1.0e2, w=0.05, x0=0.3, cfl=0.5):
    """Gaussian pulse. Returns the bar and the L2 stress error at t_end
    (cell averages against exact cell averages)."""
    bar = Bar1D(n, LENGTH, RHO, MODULUS)
    init_right_going(bar, gaussian_cell_average(bar.xf, bar.dx, x0, w, LENGTH), amp, v0)
    bar.run(t_end, cfl)
    x_end = x0 + (v0 + bar.c) * t_end
    exact = amp * gaussian_cell_average(bar.xf, bar.dx, x_end, w, LENGTH)
    err = np.sqrt(np.mean((bar.sig - exact) ** 2)) / amp
    return bar, err


def square(x, a, b):
    return np.where((x >= a) & (x < b), 1.0, 0.0)


def front_position(xc, sig, level, guess, period):
    """Position nearest to guess where sig falls through level going right
    (linear interpolation between cell centres, periodic)."""
    n = len(xc)
    dx = xc[1] - xc[0]
    best = None
    for a in range(n):
        b = (a + 1) % n
        if sig[a] >= level > sig[b]:
            x = xc[a] + dx * (sig[a] - level) / (sig[a] - sig[b])
            dist = abs((x - guess + 0.5 * period) % period - 0.5 * period)
            if best is None or dist < best[0]:
                best = (dist, x)
    if best is None:
        raise RuntimeError("front not found")
    return best[1]


# ----------------------------------------------------------------------------
# Acceptance checks for WP 1.1
# ----------------------------------------------------------------------------

def verify(v0_frac=0.25):
    """Run the WP 1.1 checks. Returns True if all pass."""
    ok = True
    c = np.sqrt(MODULUS / RHO)
    v0 = v0_frac * c
    print(f"Material: rho = {RHO} kg/m^3, M = {MODULUS:.3e} Pa, c = {c:.4f} m/s")
    print(f"Background velocity v0 = {v0:.4f} m/s ({v0_frac} c): the material")
    print("streams through the fixed mesh, so the remap moves real amounts.")

    # 1. Smooth Gaussian pulse: convergence order of the L2 stress error.
    print("\n[1] Smooth Gaussian pulse (width 0.05 m), L2 stress error")
    t_end = 0.4 * LENGTH / (c + v0)        # pulse travels 0.4 of the bar
    levels = [200, 400, 800]
    errs = []
    for n in levels:
        _, err = smooth_case(n, v0, t_end)
        errs.append(err)
        print(f"  cells {n:5d}  L2 error {err:.4e}")
    orders = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    print("  observed orders:", ", ".join(f"{o:.3f}" for o in orders))
    pass1 = min(orders) >= 1.8
    print(f"  order >= 1.8: {'PASS' if pass1 else 'FAIL'}")
    ok &= pass1

    # 2. Square (step) pulse: monotone and correct speed.
    # The linear scheme is TVD, so the only "new extrema" allowed are the
    # ones the physics and arithmetic make: (a) the hypoelastic bar is
    # slightly nonlinear, so the exact solution leaves the linear bounds by
    # O(strain) relative (about 5e-7 of the jump at 100 Pa); (b) round-off in
    # the cell lengths gives stress noise of order M * 1e-16 per step (about
    # 1e-6 Pa over this run). A 10 Pa pulse (strain 2e-7) keeps both near
    # 1e-7 of the jump; the tolerance is 1e-6 of the jump. For scale, the
    # first, staggered version overshot by 5e-2 (decision 0004).
    print("\n[2] Square stress pulse (10 Pa), monotonicity and wave speed")
    n = 800
    amp = 10.0
    a, b = 0.2, 0.4
    bar = Bar1D(n, LENGTH, RHO, MODULUS)
    init_right_going(bar, square(bar.xc, a, b), amp, v0)
    lo_s, hi_s = bar.sig.min(), bar.sig.max()
    lo_v, hi_v = bar.v.min(), bar.v.max()
    m0, p0 = bar.total_mass(), bar.total_momentum()
    pscale = (bar.cell_mass() * np.abs(bar.v)).sum()

    over = {"sig": 0.0, "v": 0.0}

    def track(br, k):
        over["sig"] = max(over["sig"], br.sig.max() - hi_s, lo_s - br.sig.min())
        over["v"] = max(over["v"], br.v.max() - hi_v, lo_v - br.v.min())

    t1 = 0.1 * LENGTH / (c + v0)
    bar.run(t1, callback=track)
    x1 = front_position(bar.xc, bar.sig, 0.5 * amp, b + (c + v0) * t1, LENGTH)
    t2 = 0.4 * LENGTH / (c + v0)
    bar.run(t2 - t1, callback=track)
    x2 = front_position(bar.xc, bar.sig, 0.5 * amp, b + (c + v0) * t2, LENGTH)
    speed = (x2 - x1) / (t2 - t1) - v0
    speed_err = abs(speed - c) / c
    tol_s = 1e-6 * (hi_s - lo_s)
    tol_v = 1e-6 * (hi_v - lo_v)
    print(f"  largest new extremum over the run: stress {over['sig']:.3e} Pa "
          f"(tol {tol_s:.1e}), velocity {over['v']:.3e} m/s (tol {tol_v:.1e})")
    pass2a = over["sig"] <= tol_s and over["v"] <= tol_v
    print(f"  monotone: {'PASS' if pass2a else 'FAIL'}")
    print(f"  front speed relative to the material {speed:.4f} m/s, "
          f"sqrt(M/rho) = {c:.4f} m/s, error {100 * speed_err:.3f} %")
    pass2b = speed_err <= 0.01
    print(f"  speed within 1 %: {'PASS' if pass2b else 'FAIL'}")
    ok &= pass2a and pass2b

    # 3. Conservation over the square pulse run.
    print("\n[3] Conservation over the square pulse run")
    dm = abs(bar.total_mass() - m0) / m0
    dp = abs(bar.total_momentum() - p0) / pscale
    print(f"  relative mass change     {dm:.3e}")
    print(f"  relative momentum change {dp:.3e}")
    pass3 = dm <= 1e-12 and dp <= 1e-12
    print(f"  round-off (<= 1e-12): {'PASS' if pass3 else 'FAIL'}")
    ok &= pass3

    print("\nWP 1.1 checks:", "ALL PASS" if ok else "FAILED")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--verify", action="store_true", help="run the WP 1.1 checks")
    p.add_argument("--v0", type=float, default=0.25,
                   help="background velocity as a fraction of c (default 0.25)")
    p.add_argument("--cells", type=int, default=400, help="cells for a single smooth run")
    args = p.parse_args(argv)
    if args.verify:
        return 0 if verify(args.v0) else 1
    c = np.sqrt(MODULUS / RHO)
    v0 = args.v0 * c
    bar, err = smooth_case(args.cells, v0, 0.4 * LENGTH / (c + v0))
    print(f"cells {args.cells}, t = {bar.t:.4e} s, L2 stress error {err:.4e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
