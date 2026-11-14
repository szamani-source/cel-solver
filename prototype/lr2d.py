#!/usr/bin/env python3
"""WP 1.3a: 2D plane-strain elastic Lagrange-plus-remap prototype.

The WP 1.1 scheme (``lr1d.Bar1D``) extended to 2D on a fixed uniform
Cartesian mesh, by dimensional splitting. Each sweep is the 1D cycle along
one axis, applied to every row (x sweep) or column (y sweep):

1. Lagrange step: faces normal to the sweep direction move with the material.
   Cell-centred MUSCL-Hancock Godunov step. In the sweep direction n (with
   tangential direction t) the plane-strain equations split into
   - a P system (v_n, s_nn) with impedance rho c_p, c_p = sqrt(M/rho),
     M = K + 4G/3,
   - an S system (v_t, s_nt) with impedance rho c_s, c_s = sqrt(G/rho),
   - s_tt and s_zz, which only change with the normal strain (stationary in
     the material frame).
   Slopes are limited in characteristic variables (MC limiter). An acoustic
   Riemann solver per system gives the face velocities v_n*, v_t* and the
   face stresses s_nn*, s_nt*. The normal strain increment is
   (h_new - h)/h_mid, as in 1D; the shear strain increment is
   dt (v_t*_right - v_t*_left)/h_mid. Stress update with the elastic matrix
   (no rotation correction yet: WP 1.3b).
2. Remap: the faces go back. Mass, both momentum components (mass-weighted,
   in mass coordinates) and the four stresses (volume-weighted) are moved
   with MUSCL fluxes over the swept region, exactly as in 1D.

Each sweep is conservative, so the 2D step is too. The two sweep orders
(x then y, y then x) are run from the same state and their conserved
quantities averaged ("symmetric" splitting, the default), so a problem
mirrored about the diagonal gives mirrored results. ``split="alternate"``
instead alternates the order each step (half the cost, not mirror-exact).
Decision 0012 records the choices.

Fields are arrays of shape (ny, nx): index [j, i] is the cell at row j
(y) and column i (x). SI units, stress positive in tension, thickness 1 m.
Boundaries per axis: "periodic" or "wall" (fixed, no-slip: both velocity
components are zero on the wall).

Run ``python3 lr2d.py --verify`` for the WP 1.3a acceptance checks
(used by tests/test_proto_2d_elastic.sh).
"""

import argparse
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from lr1d import Bar1D, limiter  # noqa: E402

NG = 2  # ghost cells on each side of a sweep line


# ----------------------------------------------------------------------------
# One sweep along the last axis of (rows, n) arrays
# ----------------------------------------------------------------------------

def pad(q, bc, odd=False):
    """Add NG ghost cells at both ends of the last axis.

    periodic: wrap around. wall: mirror image; odd=True flips the sign
    (velocities), so a Riemann solve at the wall gives zero velocity."""
    if bc == "periodic":
        return np.pad(q, ((0, 0), (NG, NG)), mode="wrap")
    p = np.pad(q, ((0, 0), (NG, NG)), mode="symmetric")
    if odd:
        p[:, :NG] *= -1.0
        p[:, -NG:] *= -1.0
    return p


def gradients(q, h):
    """One-sided differences per metre for padded cells 1 .. n+2
    (as lr1d: (q - left q) / mean width)."""
    gl = (q[:, 1:-1] - q[:, :-2]) / (0.5 * (h[:, 1:-1] + h[:, :-2]))
    gr = (q[:, 2:] - q[:, 1:-1]) / (0.5 * (h[:, 2:] + h[:, 1:-1]))
    return gl, gr


def swept_average(q, h, d):
    """lr1d.Bar1D.swept_average on padded arrays. q, h are padded (n+4
    columns), d the face displacements of the n+1 faces of the real cells.
    Returns the average of the limited reconstruction over each swept
    region (n+1 values)."""
    gl, gr = gradients(q, h)
    s = limiter(gl, gr)                      # padded cells 1 .. n+2
    qc, hc = q[:, 1:-1], h[:, 1:-1]
    a = np.abs(d)
    from_left = qc[:, :-1] + s[:, :-1] * 0.5 * (hc[:, :-1] - a)
    from_right = qc[:, 1:] - s[:, 1:] * 0.5 * (hc[:, 1:] - a)
    return np.where(d > 0.0, from_left, from_right)


def mass_flux_of(q, m, mass_flux):
    """Flux of m*q for a mass flux, q reconstructed in mass coordinates
    (lr1d.Bar1D.remap's velocity flux). q, m padded."""
    gl, gr = gradients(q, m)
    s = limiter(gl, gr)
    qc, mc = q[:, 1:-1], m[:, 1:-1]
    am = np.abs(mass_flux)
    from_left = qc[:, :-1] + s[:, :-1] * 0.5 * (mc[:, :-1] - am)
    from_right = qc[:, 1:] - s[:, 1:] * 0.5 * (mc[:, 1:] - am)
    return mass_flux * np.where(mass_flux > 0.0, from_left, from_right)


def sweep_1d(f, h0, dt, mat, bc):
    """One Lagrange step plus remap along the last axis.

    f: dict with rho, vn, vt, snn, stt, szz, snt (arrays (rows, n)); returns
    a new dict. h0: cell width in the sweep direction. mat: (M, lam, G).
    """
    M, lam, G = mat
    rows, n = f["rho"].shape
    h = np.full((rows, n + 2 * NG), h0)

    # --- Lagrange step (padded arrays; cells 1 .. n+2 get slopes) ---------
    rho = pad(f["rho"], bc)
    vn = pad(f["vn"], bc, odd=True)
    vt = pad(f["vt"], bc, odd=True)
    snn = pad(f["snn"], bc)
    snt = pad(f["snt"], bc)
    m = rho * h
    zp = np.sqrt(M * rho)
    zs = np.sqrt(G * rho)

    def characteristic(s, v, z):
        """Limited changes of s and v across each cell (cells 1 .. n+2)."""
        zc, hc = z[:, 1:-1], h[:, 1:-1]
        gs_l, gs_r = gradients(s, h)
        gv_l, gv_r = gradients(v, h)
        wr = limiter(gs_l - zc * gv_l, gs_r - zc * gv_r)
        wl = limiter(gs_l + zc * gv_l, gs_r + zc * gv_r)
        return 0.5 * (wr + wl) * hc, 0.5 * (wl - wr) / zc * hc

    def riemann(s, v, z, ds, dv, mod):
        """Hancock predictor, then the acoustic Riemann solver at the n+1
        faces of the real cells. Returns v*, s*."""
        mc, hc, zc = m[:, 1:-1], h[:, 1:-1], z[:, 1:-1]
        s_p = s[:, 1:-1] + 0.5 * dt * mod * dv / hc
        v_p = v[:, 1:-1] + 0.5 * dt * ds / mc
        sL = (s_p + 0.5 * ds)[:, :-1]
        vL = (v_p + 0.5 * dv)[:, :-1]
        zL = zc[:, :-1]
        sR = (s_p - 0.5 * ds)[:, 1:]
        vR = (v_p - 0.5 * dv)[:, 1:]
        zR = zc[:, 1:]
        zsum = zL + zR
        v_star = (zL * vL + zR * vR + (sR - sL)) / zsum
        s_star = (zR * sL + zL * sR + zL * zR * (vR - vL)) / zsum
        return v_star, s_star

    ds_n, dv_n = characteristic(snn, vn, zp)
    vn_star, snn_star = riemann(snn, vn, zp, ds_n, dv_n, M)
    ds_t, dv_t = characteristic(snt, vt, zs)
    vt_star, snt_star = riemann(snt, vt, zs, ds_t, dv_t, G)

    # Update the real cells. Face j is the left face of cell j.
    hr = h[:, NG:-NG]
    mr = m[:, NG:-NG]
    vn_new = f["vn"] + dt * (snn_star[:, 1:] - snn_star[:, :-1]) / mr
    vt_new = f["vt"] + dt * (snt_star[:, 1:] - snt_star[:, :-1]) / mr
    d = dt * vn_star
    h_new = hr + d[:, 1:] - d[:, :-1]
    # Same operation order as lr1d.Bar1D.update_stress: h_new - hr cancels,
    # so a different rounding order would show up at 1e-10 relative.
    hmid = 0.5 * (hr + h_new)
    snn_new = f["snn"] + M * (h_new - hr) / hmid
    stt_new = f["stt"] + lam * (h_new - hr) / hmid
    szz_new = f["szz"] + lam * (h_new - hr) / hmid
    snt_new = f["snt"] + G * (dt * (vt_star[:, 1:] - vt_star[:, :-1])) / hmid
    rho_new = rho[:, NG:-NG] - rho[:, NG:-NG] * (h_new - hr) / h_new

    # --- Remap back to the fixed faces ---------------------------------------
    hp = pad(h_new, bc)
    rp = pad(rho_new, bc)
    mp = rp * hp
    mass_flux = d * swept_average(rp, hp, d)
    # Increments of the Lagrangian values, in the same operation order as
    # lr1d.Bar1D.remap (an exact identity when nothing moves).
    out = {}
    dh = h_new - h0
    dmf = mass_flux[:, :-1] - mass_flux[:, 1:]
    out["rho"] = rho_new + (rho_new * dh + dmf) / h0
    m_new = out["rho"] * h0
    for key, val, odd in (("vn", vn_new, True), ("vt", vt_new, True)):
        flux = mass_flux_of(pad(val, bc, odd), mp, mass_flux)
        out[key] = val + (flux[:, :-1] - flux[:, 1:] - val * dmf) / m_new
    for key, val in (("snn", snn_new), ("stt", stt_new), ("szz", szz_new),
                     ("snt", snt_new)):
        flux = d * swept_average(pad(val, bc), hp, d)
        out[key] = val + (val * dh + flux[:, :-1] - flux[:, 1:]) / h0
    return out


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------

class Box2D:
    """Plane-strain elastic box, dimensionally split Lagrange plus remap."""

    KEYS = ("rho", "vx", "vy", "sxx", "syy", "szz", "sxy")

    def __init__(self, nx, ny, lx, ly, rho, K, G, bc=("periodic", "periodic"),
                 split="symmetric"):
        self.nx, self.ny = nx, ny
        self.lx, self.ly = lx, ly
        self.dx, self.dy = lx / nx, ly / ny
        self.K, self.G = K, G
        self.M = K + 4.0 * G / 3.0
        self.lam = K - 2.0 * G / 3.0
        self.rho0 = rho
        self.cp = math.sqrt(self.M / rho)
        self.cs = math.sqrt(G / rho)
        self.bc = bc
        self.split = split
        self.xc = (np.arange(nx) + 0.5) * self.dx
        self.yc = (np.arange(ny) + 0.5) * self.dy
        shape = (ny, nx)
        self.rho = np.full(shape, rho)
        for k in self.KEYS[1:]:
            setattr(self, k, np.zeros(shape))
        self.t = 0.0
        self.nstep = 0

    # --- state handling ----------------------------------------------------

    def state(self):
        return {k: getattr(self, k).copy() for k in self.KEYS}

    def set_state(self, s):
        for k in self.KEYS:
            setattr(self, k, s[k])

    def sweep(self, s, dt, axis):
        """One sweep on state dict s along axis 0 (x) or 1 (y)."""
        mat = (self.M, self.lam, self.G)
        if axis == 0:
            f = {"rho": s["rho"], "vn": s["vx"], "vt": s["vy"], "snn": s["sxx"],
                 "stt": s["syy"], "szz": s["szz"], "snt": s["sxy"]}
            g = sweep_1d(f, self.dx, dt, mat, self.bc[0])
            return {"rho": g["rho"], "vx": g["vn"], "vy": g["vt"], "sxx": g["snn"],
                    "syy": g["stt"], "szz": g["szz"], "sxy": g["snt"]}
        f = {"rho": s["rho"].T, "vn": s["vy"].T, "vt": s["vx"].T, "snn": s["syy"].T,
             "stt": s["sxx"].T, "szz": s["szz"].T, "snt": s["sxy"].T}
        g = sweep_1d(f, self.dy, dt, mat, self.bc[1])
        return {"rho": g["rho"].T, "vx": g["vt"].T, "vy": g["vn"].T, "sxx": g["stt"].T,
                "syy": g["snn"].T, "szz": g["szz"].T, "sxy": g["snt"].T}

    # --- integral quantities ----------------------------------------------

    def cell_mass(self):
        return self.rho * self.dx * self.dy

    def total_mass(self):
        return self.cell_mass().sum()

    def total_momentum(self):
        m = self.cell_mass()
        return (m * self.vx).sum(), (m * self.vy).sum()

    # --- driver --------------------------------------------------------------

    def stable_dt(self, cfl):
        vmax = max(np.abs(self.vx).max(), np.abs(self.vy).max())
        return cfl * min(self.dx, self.dy) / (self.cp + vmax)

    def step(self, dt):
        s = self.state()
        if self.split == "symmetric":
            a = self.sweep(self.sweep(s, dt, 0), dt, 1)
            b = self.sweep(self.sweep(s, dt, 1), dt, 0)
            rho = 0.5 * (a["rho"] + b["rho"])
            new = {"rho": rho}
            # Mass-weighted mean (ra va + rb vb)/(ra + rb), written so that it
            # is symmetric in a and b bit-for-bit and exact when va == vb.
            dr = a["rho"] - b["rho"]
            for k in ("vx", "vy"):
                new[k] = 0.5 * (a[k] + b[k]) + 0.25 * dr * (a[k] - b[k]) / rho
            for k in ("sxx", "syy", "szz", "sxy"):
                new[k] = 0.5 * (a[k] + b[k])
        else:
            order = (0, 1) if self.nstep % 2 == 0 else (1, 0)
            new = self.sweep(self.sweep(s, dt, order[0]), dt, order[1])
        self.set_state(new)
        self.t += dt
        self.nstep += 1

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
# Plane waves
# ----------------------------------------------------------------------------

# Material: the WP 1.1 sand-like solid, M = 50 MPa, with Poisson ratio 0.3.
RHO = 2000.0
MODULUS = 50.0e6                       # M = K + 4G/3, as lr1d.MODULUS
NU = 0.3
G_MOD = MODULUS * (1.0 - 2.0 * NU) / (2.0 * (1.0 - NU))
K_MOD = MODULUS - 4.0 * G_MOD / 3.0
LENGTH = 1.0


def gauss_second_antiderivative(u, w):
    """F with F'' = exp(-(u/w)^2)."""
    erfv = np.vectorize(math.erf)
    return 0.5 * math.sqrt(math.pi) * w * (u * erfv(u / w)
                                           + w / math.sqrt(math.pi) * np.exp(-(u / w) ** 2))


def gauss_first_antiderivative(u, w):
    erfv = np.vectorize(math.erf)
    return 0.5 * math.sqrt(math.pi) * w * erfv(u / w)


def wave_profile(box, kx, ky, s0, w):
    """Exact cell averages of the periodic profile g(s) = sum_k
    exp(-((s - s0 - k P)/w)^2), s = kx x + ky y (kx, ky in {0, 1}),
    period P = LENGTH. Pulse width along the wave normal is w/|k|."""
    P = LENGTH
    X, Y = np.meshgrid(np.arange(box.nx) * box.dx, np.arange(box.ny) * box.dy)
    s = kx * X + ky * Y
    u = (s - s0 + 0.5 * P) % P - 0.5 * P      # nearest image near 0
    out = np.zeros_like(u)
    for k in range(-2, 3):
        uk = u - k * P
        if kx and ky:
            h = box.dx                         # square cells, dx = dy
            F = gauss_second_antiderivative
            out += (F(uk + 2 * h, w) - 2 * F(uk + h, w) + F(uk, w)) / h ** 2
        else:
            h = box.dx if kx else box.dy
            F = gauss_first_antiderivative
            out += (F(uk + h, w) - F(uk, w)) / h
    return out


def set_plane_wave(box, kind, angle_deg, eps, v0=(0.0, 0.0)):
    """Right-going plane wave along n = (cos a, sin a), a in {0, 45, 90}.

    kind "P": displacement along n, strain eps n n. kind "S": displacement
    along t = (-sin a, cos a), engineering shear strain eps. eps is an
    (ny, nx) array (the profile times the amplitude). v0 is a uniform
    streaming velocity."""
    a = math.radians(angle_deg)
    n = np.round(np.array([math.cos(a), math.sin(a)]), 15)
    t = np.array([-n[1], n[0]])
    if kind == "P":
        c = box.cp
        sig = box.lam * np.eye(2) + 2.0 * box.G * np.outer(n, n)
        szz = box.lam
        pol = n
    else:
        c = box.cs
        sig = box.G * (np.outer(n, t) + np.outer(t, n))
        szz = 0.0
        pol = t
    box.sxx = sig[0, 0] * eps
    box.syy = sig[1, 1] * eps
    box.sxy = sig[0, 1] * eps
    box.szz = szz * eps
    box.vx = v0[0] - c * eps * pol[0]
    box.vy = v0[1] - c * eps * pol[1]
    return n, t, c


def measured(box, kind, n, t):
    """The stress component a plane wave carries: s_nn (P) or s_nt (S)."""
    S = [[box.sxx, box.sxy], [box.sxy, box.syy]]
    b = n if kind == "P" else t
    return sum(n[i] * S[i][j] * b[j] for i in range(2) for j in range(2))


def phase_shift(prof, period, nmode=1):
    """Position of the profile from the phase of its first Fourier mode."""
    k = np.arange(len(prof))
    c = (prof * np.exp(-2j * np.pi * nmode * k / len(prof))).sum()
    return -np.angle(c) / (2 * np.pi * nmode) * period


def smooth_case(kind, angle, n, v0_frac=0.25, amp=1.0e-6, w=0.1, travel=0.4,
                split="symmetric"):
    """Smooth Gaussian plane wave. Returns (L2 error, wave speed relative to
    the material, exact speed, run time per cell per step, box).

    Axis-aligned waves use an n x 4 strip; the 45 degree wave an n x n box.
    The material streams at v0_frac * c_p along the diagonal."""
    if angle == 45:
        nx, ny, kx, ky = n, n, 1, 1
    elif angle == 0:
        nx, ny, kx, ky = n, 4, 1, 0
    else:
        nx, ny, kx, ky = 4, n, 0, 1
    box = Box2D(nx, ny, nx * LENGTH / n, ny * LENGTH / n, RHO, K_MOD, G_MOD, split=split)
    knorm = math.hypot(kx, ky)
    wv = w * knorm                              # width in s = kx x + ky y
    vd = v0_frac * box.cp / math.sqrt(2.0)
    v0 = (vd, vd)
    prof0 = wave_profile(box, kx, ky, 0.3, wv)
    nvec, tvec, c = set_plane_wave(box, kind, angle, amp * prof0, v0)
    vn0 = v0[0] * nvec[0] + v0[1] * nvec[1]
    period_n = LENGTH / knorm                   # period along the normal
    t_end = travel * period_n / (c + vn0)
    t0 = time.perf_counter()
    steps = box.run(t_end)
    cost = (time.perf_counter() - t0) / (steps * nx * ny)
    shift_s = (c + vn0) * t_end * knorm         # shift in s
    exact_sigma = measured_exact(box, kind, nvec, tvec, amp)
    exact = exact_sigma * amp * wave_profile(box, kx, ky, 0.3 + shift_s, wv)
    num = measured(box, kind, nvec, tvec)
    err = np.sqrt(np.mean((num - exact) ** 2)) / (abs(exact_sigma) * amp)
    # Speed from the phase of the first Fourier mode of the 1D profile,
    # profile taken along s (diagonals i + j for 45 degrees).
    prof = line_profile(num, kx, ky)
    ref = line_profile(exact_sigma * prof0 * amp, kx, ky)
    ds = phase_shift(prof, LENGTH) - phase_shift(ref, LENGTH)
    ds = (ds + 0.5 * LENGTH) % LENGTH - 0.5 * LENGTH
    # shift_s may be more than one period: add the whole periods back
    whole = shift_s - ((shift_s + 0.5 * LENGTH) % LENGTH - 0.5 * LENGTH)
    speed = (whole + ds) / knorm / t_end - vn0
    return err, speed, c, cost, box


def measured_exact(box, kind, n, t, amp):
    """s_nn per unit strain (P: M) or s_nt per unit shear strain (S: G)."""
    return box.M if kind == "P" else box.G


def line_profile(q, kx, ky):
    """Average of q along lines of constant s = kx x + ky y."""
    if kx and ky:
        ny, nx = q.shape
        J, I = np.indices(q.shape)
        k = (I + J) % nx
        return np.bincount(k.ravel(), q.ravel(), nx) / np.bincount(k.ravel(), None, nx)
    return q.mean(axis=0) if kx else q.mean(axis=1)


# ----------------------------------------------------------------------------
# Acceptance checks for WP 1.3a
# ----------------------------------------------------------------------------

def check_convergence(kind, angle, levels, report):
    errs, speeds, costs = [], [], []
    for n in levels:
        err, speed, c, cost, _ = smooth_case(kind, angle, n)
        errs.append(err)
        speeds.append(speed)
        costs.append(cost)
        print(f"  cells {n:4d}  L2 error {err:.4e}  speed {speed:.4f} m/s")
    orders = [math.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    serr = abs(speeds[-1] - c) / c
    ok_o = min(orders) >= 1.7
    ok_s = serr <= 0.01
    print(f"  orders {', '.join(f'{o:.3f}' for o in orders)} (>= 1.7: "
          f"{'PASS' if ok_o else 'FAIL'}); speed error {100 * serr:.3f} % of "
          f"{c:.4f} m/s (<= 1 %: {'PASS' if ok_s else 'FAIL'})")
    report.append((kind, angle, errs, orders, costs[-1]))
    return ok_o and ok_s


def check_reduces_to_1d():
    """A P-wave along x on an n x 2 mesh against lr1d.Bar1D, same dt."""
    from lr1d import gaussian_cell_average, init_right_going
    n = 200
    c = math.sqrt(MODULUS / RHO)
    v0 = 0.25 * c
    bar = Bar1D(n, LENGTH, RHO, MODULUS)
    init_right_going(bar, gaussian_cell_average(bar.xf, bar.dx, 0.3, 0.05, LENGTH),
                     100.0, v0)
    box = Box2D(n, 2, LENGTH, 2 * LENGTH / n, RHO, K_MOD, G_MOD)
    box.sxx = np.tile(bar.sig, (2, 1))
    box.syy = box.lam / box.M * box.sxx
    box.szz = box.syy.copy()
    box.vx = np.tile(bar.v, (2, 1))
    t_end = 0.4 * LENGTH / (c + v0)
    nsteps = int(np.ceil(t_end / bar.stable_dt(0.5)))
    dt = t_end / nsteps
    for _ in range(nsteps):
        bar.step(dt)
        box.step(dt)
    ds = np.abs(box.sxx - bar.sig).max() / np.abs(bar.sig).max()
    dv = np.abs(box.vx - bar.v).max() / np.abs(bar.v).max()
    dr = np.abs(box.rho - bar.rho).max() / RHO
    side = max(np.abs(box.vy).max(), np.abs(box.sxy).max())
    worst = max(ds, dv, dr)
    ok = worst <= 1e-10 and side == 0.0
    print(f"  {nsteps} steps; largest relative difference: stress {ds:.2e}, "
          f"velocity {dv:.2e}, density {dr:.2e}; vy, sxy stay {side:.1e}")
    print(f"  matches lr1d to 1e-10: {'PASS' if ok else 'FAIL'}")
    return ok


def check_mirror():
    """A P pulse off the diagonal, run and mirrored about the diagonal."""
    n = 64
    c = math.sqrt(MODULUS / RHO)

    def make(mirror):
        box = Box2D(n, n, LENGTH, LENGTH, RHO, K_MOD, G_MOD)
        X, Y = np.meshgrid(box.xc, box.yc)
        if mirror:
            X, Y = Y, X
        r2 = (X - 0.3) ** 2 + (Y - 0.6) ** 2
        p = np.exp(-r2 / 0.08 ** 2)
        box.sxx = 100.0 * p
        box.syy = 60.0 * p
        box.szz = 40.0 * p
        box.sxy = 30.0 * np.exp(-((X - 0.5) ** 2 + (Y - 0.4) ** 2) / 0.1 ** 2)
        box.vx = 0.3 * c - 0.02 * p
        box.vy = 0.1 * c + 0.01 * np.exp(-((X - 0.7) ** 2 + Y ** 2) / 0.1 ** 2)
        if mirror:
            box.sxx, box.syy = box.syy, box.sxx
            box.vx, box.vy = box.vy, box.vx
        return box

    a, b = make(False), make(True)
    t_end = 0.3 * LENGTH / c
    a.run(t_end)
    b.run(t_end)
    pairs = (("rho", "rho"), ("vx", "vy"), ("vy", "vx"), ("sxx", "syy"),
             ("syy", "sxx"), ("szz", "szz"), ("sxy", "sxy"))
    worst = 0.0
    for ka, kb in pairs:
        qa, qb = getattr(a, ka), getattr(b, kb).T
        scale = max(np.abs(qa).max(), 1e-300)
        worst = max(worst, np.abs(qa - qb).max() / scale)
    ok = worst <= 1e-12
    print(f"  largest relative difference from the mirror image {worst:.2e}")
    print(f"  mirror-symmetric to 1e-12: {'PASS' if ok else 'FAIL'}")
    return ok


def check_step(angle, n, startup_cells=0.0):
    """Square P pulse along angle, material streaming at 0.25 c_p along the
    diagonal. No new extrema above 1 % of each field's jump; conservation.

    Extrema in the first startup_cells cells of wave travel (along the
    normal) are reported, not checked. Along x there is no start-up
    overshoot (startup_cells 0). At 45 degrees a sharp step overshoots by
    about 2.7 % in the first cell of travel and below 0.3 % after three
    (decision 0013; the 13 % of decision 0012 came from a jagged initial
    step, see below)."""
    box = Box2D(n, n, LENGTH, LENGTH, RHO, K_MOD, G_MOD)
    kx, ky = (1, 1) if angle == 45 else (1, 0)
    if angle == 45:
        # s = x + y from integer indices. Computed from the cell centres in
        # floating point, cells of one diagonal fall on both sides of the
        # threshold, the front is jagged at the cell scale and radiates
        # genuine 2D waves (most of the 13 % in decision 0012).
        J, I = np.indices((n, n))
        s = ((I + J + 1) % n) * box.dx
    else:
        X, Y = np.meshgrid(box.xc, box.yc)
        s = (kx * X) % LENGTH
    eps = 1.0e-7 * np.where((s >= 0.2) & (s < 0.5), 1.0, 0.0)
    vd = 0.25 * box.cp / math.sqrt(2.0)
    set_plane_wave(box, "P", angle, eps, (vd, vd))
    keys = ("vx", "vy", "sxx", "syy", "szz", "sxy")
    lo = {k: getattr(box, k).min() for k in keys}
    hi = {k: getattr(box, k).max() for k in keys}
    over = {k: 0.0 for k in keys}
    over_start = {k: 0.0 for k in keys}
    m0 = box.total_mass()
    px0, py0 = box.total_momentum()
    pscale = (box.cell_mass() * np.hypot(box.vx, box.vy)).sum()
    hn = box.dx / math.hypot(kx, ky)          # cell size along the normal

    def track(b, k):
        w = over_start if box.cp * b.t < startup_cells * hn else over
        for q in keys:
            a = getattr(b, q)
            w[q] = max(w[q], a.max() - hi[q], lo[q] - a.min())

    t_end = 0.4 * LENGTH / math.hypot(kx, ky) / (box.cp + 0.25 * box.cp)
    box.run(t_end, callback=track)
    worst = 0.0
    worst_start = 0.0
    final = 0.0
    for q in keys:
        jump = hi[q] - lo[q]
        if jump > 0.0:
            worst = max(worst, over[q] / jump)
            worst_start = max(worst_start, over_start[q] / jump)
            a = getattr(box, q)
            final = max(final, (a.max() - hi[q]) / jump, (lo[q] - a.min()) / jump)
    ok_e = worst <= 0.01
    dm = abs(box.total_mass() - m0) / m0
    px, py = box.total_momentum()
    dpx, dpy = abs(px - px0) / pscale, abs(py - py0) / pscale
    ok_c = dm <= 1e-12 and dpx <= 1e-12 and dpy <= 1e-12
    verdict = "PASS" if ok_e else "FAIL"
    if startup_cells > 0.0:
        print(f"  {angle} deg, {n} x {n}: start-up (first {startup_cells:g} cells of "
              f"travel) {100 * worst_start:.3f} % of the jump (reported only)")
    after = " after the start-up" if startup_cells > 0.0 else " over the run"
    print(f"  {angle} deg, {n} x {n}: largest new extremum{after} "
          f"{100 * worst:.3f} % of the jump, at the end {100 * final:.3f} % "
          f"(<= 1 %: {verdict})")
    print(f"  mass change {dm:.1e}, momentum change x {dpx:.1e}, y {dpy:.1e} "
          f"(<= 1e-12: {'PASS' if ok_c else 'FAIL'})")
    return ok_e and ok_c


def check_wall():
    """Not a WP criterion: a wall box keeps its mass and stops at the wall."""
    n = 40
    box = Box2D(n, n, LENGTH, LENGTH, RHO, K_MOD, G_MOD, bc=("wall", "wall"))
    X, Y = np.meshgrid(box.xc, box.yc)
    box.vx = 1e-3 * box.cp * np.exp(-((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 0.1 ** 2)
    m0 = box.total_mass()
    box.run(1.2 * LENGTH / box.cp)
    dm = abs(box.total_mass() - m0) / m0
    ok = dm <= 1e-12 and np.isfinite(box.sxx).all()
    print(f"  wall box: mass change {dm:.1e} ({'PASS' if ok else 'FAIL'})")
    return ok


def verify(quick=False):
    ok = True
    print(f"Material: rho {RHO} kg/m^3, K {K_MOD:.4e} Pa, G {G_MOD:.4e} Pa "
          f"(nu {NU}), c_p {math.sqrt(MODULUS / RHO):.4f} m/s, "
          f"c_s {math.sqrt(G_MOD / RHO):.4f} m/s")
    print("Smooth waves: Gaussian width 0.1 m along the normal, material streaming")
    print("at 0.25 c_p along the diagonal; meshes 100, 200, 400 cells per metre.")
    levels = [50, 100, 200] if quick else [100, 200, 400]
    report = []
    for kind, angle in (("P", 0), ("P", 90), ("P", 45), ("S", 0), ("S", 45)):
        print(f"\n[{kind}-wave, {angle} deg]")
        ok &= check_convergence(kind, angle, levels, report)
    print("\n[Reduces to 1D: P-wave along x vs lr1d.Bar1D, 200 cells, 0.25 c streaming]")
    ok &= check_reduces_to_1d()
    print("\n[Mirror symmetry about the diagonal, 64 x 64]")
    ok &= check_mirror()
    print("\n[Step P-wave along x, material streaming 0.25 c_p diagonally]")
    ok &= check_step(0, 100 if quick else 200)
    print("\n[Step P-wave at 45 degrees, same streaming; start-up of 10 cells "
          "reported (decision 0013)]")
    ok &= check_step(45, 100 if quick else 200, startup_cells=10.0)
    print("\n[Fixed walls]")
    ok &= check_wall()
    print("\nFor the log: 45 degree vs axis-aligned L2 errors on the finest mesh")
    for kind, angle, errs, orders, cost in report:
        print(f"  {kind} {angle:2d} deg: error {errs[-1]:.3e}, "
              f"{1e9 * cost:.0f} ns per cell per step")
    print("\nWP 1.3a checks:", "ALL PASS" if ok else "FAILED")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--verify", action="store_true", help="run the WP 1.3a checks")
    p.add_argument("--quick", action="store_true", help="coarser meshes (not the test)")
    args = p.parse_args(argv)
    if args.verify:
        return 0 if verify(args.quick) else 1
    for kind, angle in (("P", 0), ("P", 45), ("S", 45)):
        err, speed, c, cost, _ = smooth_case(kind, angle, 100)
        print(f"{kind} {angle} deg, 100 cells: L2 error {err:.3e}, speed {speed:.3f} "
              f"(exact {c:.3f}), {1e9 * cost:.0f} ns per cell per step")
    return 0


if __name__ == "__main__":
    sys.exit(main())
