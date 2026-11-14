#!/usr/bin/env python3
"""WP 1.1: 1D elastic Lagrange-plus-remap prototype.

A periodic 1D bar on a fixed uniform (Eulerian) mesh. Each cycle is

1. Lagrange step: the mesh moves with the material. Central-difference time
   integration in velocity-Verlet form, so velocity, stress and density are all
   known at whole time levels when the remap starts.
2. Remap: the nodes go back to their fixed positions. Mass, stress (volume
   weighted) and nodal momentum are moved between cells with second-order
   monotone MUSCL fluxes (van Leer limiter). Momentum is remapped on the dual
   mesh with mass fluxes from the half-index shift (Benson 1992, section 4).

Staggered grid: velocity at nodes, density and stress at cell centres.
Cell j spans node j (left) to node j+1 (right). Area is 1 m^2.

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


def van_leer_slope(q, h):
    """Limited slope (per metre) of cell values q on cells of length h.

    Harmonic (van Leer) limiter: zero at extrema, so the reconstruction
    never leaves the range of the neighbouring cell values.
    """
    dl = (q - left(q)) / (0.5 * (h + left(h)))
    dr = (right(q) - q) / (0.5 * (h + right(h)))
    prod = dl * dr
    s = np.zeros_like(q)
    mask = prod > 0.0
    s[mask] = 2.0 * prod[mask] / (dl[mask] + dr[mask])
    return s


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------

class Bar1D:
    """Periodic elastic bar, Lagrange step plus MUSCL remap."""

    def __init__(self, n, length, rho, M, alpha_q=1.0):
        self.n = n
        self.length = length
        self.dx = length / n
        self.M = M
        self.rho0 = rho
        self.c = np.sqrt(M / rho)
        self.alpha_q = alpha_q

        # Fixed Eulerian mesh
        self.xn = self.dx * np.arange(n)          # node positions
        self.xc = self.xn + 0.5 * self.dx         # cell centres

        # State on the Eulerian mesh
        self.rho = np.full(n, rho)                # cell density
        self.sig = np.zeros(n)                    # cell stress
        self.v = np.zeros(n)                      # node velocity
        self.t = 0.0

    # --- integral quantities ------------------------------------------------

    def cell_mass(self):
        return self.rho * self.dx

    def node_mass(self):
        m = self.cell_mass()
        return 0.5 * (left(m) + m)

    def total_mass(self):
        return self.cell_mass().sum()

    def total_momentum(self):
        return (self.node_mass() * self.v).sum()

    # --- Lagrange step pieces -------------------------------------------------

    def viscosity(self, v, h, rho):
        """Limited linear artificial viscosity (Christensen-type).

        q = alpha * rho * c * h * strain_rate * (1 - psi). psi is a limiter on the
        ratio of neighbouring strain rates: psi -> 1 where the velocity field is
        smooth (q vanishes to second order), psi = 0 at a jump (full damping).
        """
        if self.alpha_q == 0.0:
            return np.zeros_like(h)
        erate = (right(v) - v) / h
        tiny = 1e-300
        safe = np.where(np.abs(erate) > tiny, erate, tiny)
        rl = left(erate) / safe
        rr = right(erate) / safe
        psi = np.clip(np.minimum.reduce([0.5 * (rl + rr), 2.0 * rl, 2.0 * rr]), 0.0, 1.0)
        return self.alpha_q * rho * self.c * h * erate * (1.0 - psi)

    @staticmethod
    def acceleration(s, m_node):
        """Node acceleration from cell stress s (stress + viscosity)."""
        return (s - left(s)) / m_node

    def lagrange_step(self, dt):
        """Move the mesh with the material for one step.

        Returns node displacements d, Lagrangian cell lengths h, and the
        Lagrangian cell density, stress and node velocity at the new time.
        """
        mc = self.cell_mass()
        mn = self.node_mass()
        h0 = np.full(self.n, self.dx)

        s0 = self.sig + self.viscosity(self.v, h0, self.rho)
        vh = self.v + 0.5 * dt * self.acceleration(s0, mn)

        d = dt * vh                                  # node displacement
        h = self.dx + right(d) - d                   # new cell length
        hmid = self.dx + 0.5 * (right(d) - d)        # mid-step cell length
        sig = self.sig + self.M * (right(d) - d) / hmid   # hypoelastic update
        rho = mc / h

        s1 = sig + self.viscosity(vh, hmid, rho)
        v = vh + 0.5 * dt * self.acceleration(s1, mn)
        return d, h, rho, sig, v

    # --- remap ---------------------------------------------------------------

    @staticmethod
    def swept_average(q, h, d):
        """Average of the MUSCL reconstruction of q over the region swept by
        each node when it moves back from x_i = X_i + d_i to X_i.

        Node i sits between cell i-1 (left) and cell i (right). If d_i > 0 the
        region [X_i, X_i + d_i] is the right end of Lagrangian cell i-1, which
        belongs to Eulerian cell i after the remap; if d_i < 0 it is the left
        end of Lagrangian cell i, which belongs to Eulerian cell i-1.
        """
        s = van_leer_slope(q, h)
        a = np.abs(d)
        from_left = left(q) + left(s) * 0.5 * (left(h) - a)
        from_right = q - s * 0.5 * (h - a)
        return np.where(d > 0.0, from_left, from_right)

    def remap(self, d, h, rho, sig, v):
        """Return Lagrangian fields to the fixed mesh (conservative)."""
        mc = rho * h
        mn = 0.5 * (left(mc) + mc)

        # Fluxes through node i, positive left-to-right (cell i-1 -> cell i).
        mass_flux = d * self.swept_average(rho, h, d)
        sigvol_flux = d * self.swept_average(sig, h, d)

        mc_new = mc + mass_flux - right(mass_flux)
        sigvol_new = sig * h + sigvol_flux - right(sigvol_flux)

        # Momentum on the dual mesh. Dual face j sits at the centre of cell j,
        # between node j and node j+1. Half-index shift of the mass flux.
        dual_flux = 0.5 * (mass_flux + right(mass_flux))
        # MUSCL value of v at dual face j from the donor node, mass coordinates.
        dmn_l = 0.5 * (left(mn) + mn)
        dmn_r = 0.5 * (mn + right(mn))
        dvl = (v - left(v)) / dmn_l
        dvr = (right(v) - v) / dmn_r
        prod = dvl * dvr
        slope = np.zeros_like(v)
        mask = prod > 0.0
        slope[mask] = 2.0 * prod[mask] / (dvl[mask] + dvr[mask])
        a = np.abs(dual_flux)
        v_from_left = v + slope * 0.5 * (mn - a)                       # donor node j
        v_from_right = right(v) - right(slope) * 0.5 * (right(mn) - a)  # donor node j+1
        v_face = np.where(dual_flux > 0.0, v_from_left, v_from_right)
        mom_flux = dual_flux * v_face

        p_new = mn * v + left(mom_flux) - mom_flux
        mn_new = 0.5 * (left(mc_new) + mc_new)

        self.rho = mc_new / self.dx
        self.sig = sigvol_new / self.dx
        self.v = p_new / mn_new

    # --- driver ----------------------------------------------------------------

    def stable_dt(self, cfl):
        return cfl * self.dx / (self.c + np.abs(self.v).max())

    def step(self, dt):
        self.remap(*self.lagrange_step(dt))
        self.t += dt

    def run(self, t_end, cfl=0.5, callback=None):
        """Advance to t_end with a constant step chosen to land on t_end."""
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


def gaussian(x, x0, w, period):
    """Periodic Gaussian (sum of the nearest images)."""
    xi = (x - x0 + 0.5 * period) % period - 0.5 * period
    return np.exp(-(xi / w) ** 2)


def init_right_going(bar, shape, amp, v0):
    """Right-going wave: sig = amp*shape, v = v0 - sig/(rho c).

    For rho v_t = sig_x, sig_t = M v_x a right-going wave f(x - ct) has
    v = -f/(rho c). v0 is a uniform background velocity, which makes the
    material stream through the fixed mesh so the remap does real work.
    The node velocity is the average of the two adjacent cell values.
    """
    bar.sig = amp * shape(bar.xc)
    sig_nodes = amp * shape(bar.xn)
    bar.v = v0 - sig_nodes / (bar.rho0 * bar.c)


def smooth_case(n, v0, t_end, amp=1.0e2, w=0.05, x0=0.3, cfl=0.5):
    """Gaussian pulse. Returns the bar and the L2 stress error at t_end."""
    bar = Bar1D(n, LENGTH, RHO, MODULUS)
    shape = lambda x: gaussian(x, x0, w, LENGTH)
    init_right_going(bar, shape, amp, v0)
    m0, p0 = bar.total_mass(), bar.total_momentum()
    bar.run(t_end, cfl)
    exact = amp * gaussian(bar.xc, x0 + (v0 + bar.c) * t_end, w, LENGTH)
    err = np.sqrt(np.mean((bar.sig - exact) ** 2)) / amp
    return bar, err, m0, p0


def square(x, a, b):
    return np.where((x >= a) & (x < b), 1.0, 0.0)


def front_position(xc, sig, amp, guess, period):
    """Position where sig first drops below amp/2 going right from guess
    (linear interpolation between cell centres)."""
    n = len(xc)
    dx = xc[1] - xc[0]
    j = int(np.floor(((guess - xc[0]) % period) / dx))
    # search a window around the guess for the falling crossing
    for k in range(-n // 4, n // 4):
        a = (j + k) % n
        b = (a + 1) % n
        if sig[a] >= 0.5 * amp > sig[b]:
            frac = (sig[a] - 0.5 * amp) / (sig[a] - sig[b])
            return xc[a] + frac * dx + (period if b < a and False else 0.0)
    raise RuntimeError("front not found")


# ----------------------------------------------------------------------------
# Acceptance checks for WP 1.1
# ----------------------------------------------------------------------------

def verify(v0_frac=0.25):
    """Run the three WP 1.1 checks. Returns True if all pass."""
    ok = True
    c = np.sqrt(MODULUS / RHO)
    v0 = v0_frac * c
    print(f"Material: rho = {RHO} kg/m^3, M = {MODULUS:.3e} Pa, c = {c:.4f} m/s")
    print(f"Background velocity v0 = {v0:.4f} m/s ({v0_frac} c)")

    # 1. Smooth Gaussian pulse: convergence order of the L2 stress error.
    print("\n[1] Smooth Gaussian pulse, L2 stress error")
    t_end = 0.4 * LENGTH / (c + v0)        # pulse travels 0.4 of the bar
    levels = [200, 400, 800]
    errs = []
    for n in levels:
        bar, err, m0, p0 = smooth_case(n, v0, t_end)
        errs.append(err)
        print(f"  cells {n:5d}  L2 error {err:.4e}")
    orders = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    print("  observed orders:", ", ".join(f"{o:.3f}" for o in orders))
    pass1 = min(orders) >= 1.8
    print(f"  order >= 1.8: {'PASS' if pass1 else 'FAIL'}")
    ok &= pass1

    # 2. Step (square) wave: monotone and correct speed.
    print("\n[2] Square stress pulse, monotonicity and wave speed")
    n = 800
    amp = 1.0e2
    a, b = 0.2, 0.4
    bar = Bar1D(n, LENGTH, RHO, MODULUS)
    init_right_going(bar, lambda x: square(x, a, b), amp, v0)
    lo_s, hi_s = bar.sig.min(), bar.sig.max()
    lo_v, hi_v = bar.v.min(), bar.v.max()
    m0, p0 = bar.total_mass(), bar.total_momentum()
    pscale = (bar.node_mass() * np.abs(bar.v)).sum()

    over = {"sig": 0.0, "v": 0.0}

    def track(br, k):
        over["sig"] = max(over["sig"], br.sig.max() - hi_s, lo_s - br.sig.min())
        over["v"] = max(over["v"], br.v.max() - hi_v, lo_v - br.v.min())

    t1 = 0.1 * LENGTH / (c + v0)
    bar.run(t1, callback=track)
    x1 = front_position(bar.xc, bar.sig, amp, b + (c + v0) * t1, LENGTH)
    t2 = 0.4 * LENGTH / (c + v0)
    bar.run(t2 - t1, callback=track)
    x2 = front_position(bar.xc, bar.sig, amp, b + (c + v0) * t2, LENGTH)
    speed = (x2 - x1) / (t2 - t1) - v0
    speed_err = abs(speed - c) / c
    tol_s = 1e-10 * amp
    tol_v = 1e-10 * (hi_v - lo_v)
    print(f"  max new extremum: stress {over['sig']:.3e} Pa (tol {tol_s:.1e}), "
          f"velocity {over['v']:.3e} m/s (tol {tol_v:.1e})")
    pass2a = over["sig"] <= tol_s and over["v"] <= tol_v
    print(f"  monotone: {'PASS' if pass2a else 'FAIL'}")
    print(f"  front speed relative to material {speed:.4f} m/s, "
          f"sqrt(M/rho) = {c:.4f}, error {100 * speed_err:.3f} %")
    pass2b = speed_err <= 0.01
    print(f"  speed within 1 %: {'PASS' if pass2b else 'FAIL'}")
    ok &= pass2a and pass2b

    # 3. Conservation (on the square pulse run, the harder case).
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
    bar, err, m0, p0 = smooth_case(args.cells, v0, 0.4 * LENGTH / (c + v0))
    print(f"cells {args.cells}, t = {bar.t:.4e} s, L2 stress error {err:.4e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
