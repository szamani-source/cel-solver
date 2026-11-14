#!/usr/bin/env python3
"""WP 3.1p: 2D tank under gravity with smooth or frictional walls (settling).

Stage 1 of target T0 (docs/spec/target-problems.md): a closed rectangular tank,
completely full of cohesionless Mohr-Coulomb material, is loaded by gravity and
settles to rest. The scheme is the WP 1.3b one (``lr2d_plastic.py``: symmetric
dimensional splitting, cell-centred MUSCL-Hancock Godunov Lagrange step with
acoustic Riemann solvers, Jaumann/Hughes-Winget stress update, mass-weighted
MUSCL remap). This file adds, per sweep:

1. **Gravity, well-balanced.** Gravity acts in -y. In the y sweep each column
   gets the discrete weight profile G (G at a face = g times the mass per unit
   area below it, built from the cell densities at the start of the sweep;
   G at a centre = face below + half the cell's own weight). The Lagrange
   step reconstructs, predicts and solves the Riemann problems for
   ``s_hat = s_yy - G`` instead of s_yy, with no separate source term:

       dv/dt = (1/m) d(s_yy)/dy - g = (1/m) d(s_yy - G)/dy,

   which is exact because G jumps by rho g h across a cell. A column in
   equilibrium has s_hat constant, so the slopes, the Hancock predictor and the
   Riemann solutions are all exactly uniform and nothing moves (round-off
   only). Walls reflect s_hat evenly; the free top is the zero-traction
   condition s_yy = 0, that is s_hat = -G(top). The Riemann solutions give the
   face stress s_hat*, and the actual face stress is s_hat* + G(face).
2. **Boundaries.** Each end of a sweep is one of
   - ``"smooth"``: rigid wall, zero normal velocity, zero shear traction;
   - ``"friction"``: rigid wall, zero normal velocity, Coulomb friction: the
     sticking solution (zero tangential velocity) is kept if its shear
     traction is at most tan(delta_w) times the wall pressure, otherwise the
     shear is set to that limit with the sign of the sticking shear and the
     tangential velocity follows from the interior characteristic, so the
     friction opposes the sliding;
   - ``"free"``: zero traction (the top of a completely full tank).
   At a wall whose sticking normal stress would be tensile the material
   separates: zero normal traction (and zero shear) instead. Ghost cells are
   mirrors with the matching parity (only the slopes of the cells next to the
   boundary use them); the boundary face's Riemann solution is then set from
   the interior state alone with the one-sided relation
   s* = s_I + z_I (v* - v_I) (sign by side). The remap has no mass flux
   through any boundary face, so mass is exact. At the free top this means the
   top cells lose a little density as the column settles (about 1 mm of
   settlement over 1.2 m); WP 3.2p replaces this with a real free surface.
3. **Mohr-Coulomb** from ``constitutive/batch.py`` (vectorised in this
   package, checked against ``constitutive/models.py``).
4. **Settling controls.** Gravity is ramped from 0 to g over ``t_ramp`` with
   a smooth (1 - cos) profile, and velocity-proportional damping
   ``v <- v / (1 + alpha dt)`` after each step (``damping = alpha``; set 0
   for discharge). Settled: kinetic energy at most ``ke_tol`` times the
   elastic strain energy, after the ramp.

``python3 silo2d.py --verify`` runs the WP 3.1p checks
(tests/test_proto_silo_settle.sh). SI units, stress positive in tension,
thickness 1 m, arrays (ny, nx) as in lr2d.py.
"""

import argparse
import math
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "constitutive"))

import lr2d_plastic as lp  # noqa: E402
from lr1d import limiter  # noqa: E402
from lr1d_plastic import superbee  # noqa: E402
from lr2d_plastic import NG, PERM_Y, gradients, swept_average, mass_flux_of  # noqa: E402
from lr2d_plastic import rotate_inplane  # noqa: E402
import batch  # noqa: E402
import models  # noqa: E402

KINDS = ("smooth", "friction", "free")


# ----------------------------------------------------------------------------
# Boundary helpers
# ----------------------------------------------------------------------------

def parity(kind, system):
    """Ghost parities (stress, velocity), +1 even, -1 odd, for the normal
    ("n") or tangential ("t") wave system at a boundary of this kind. The
    frictional wall uses the sticking-wall ghosts; the friction law is applied
    to the boundary face afterwards."""
    if kind == "free":
        return -1, 1
    if system == "n" or kind == "friction":
        return 1, -1
    return -1, 1                      # smooth wall, tangential system


def pad_ends(q, par_lo, par_hi, c_lo=0.0, c_hi=0.0):
    """Mirror ghosts along axis 1; odd parity reflects about c (per row)."""
    p = mirror(q)
    if par_lo < 0:
        p[:, :NG] = 2.0 * np.reshape(c_lo, (-1, 1)) - p[:, :NG]
    if par_hi < 0:
        p[:, -NG:] = 2.0 * np.reshape(c_hi, (-1, 1)) - p[:, -NG:]
    return p


def mirror(q):
    """np.pad(q, NG along axis 1, mode="symmetric"), without np.pad's overhead."""
    return np.concatenate([q[:, NG - 1::-1], q, q[:, :-NG - 1:-1]], axis=1)


def weight_profile(rho, h, g):
    """Discrete weight profile along axis 1 of (rows, n) densities: G at the
    n + 1 faces (0 at the first) and at the n centres (Pa)."""
    w = g * rho * h
    Gf = np.zeros((rho.shape[0], rho.shape[1] + 1))
    Gf[:, 1:] = np.cumsum(w, axis=1)
    Gc = Gf[:, :-1] + 0.5 * w
    return Gf, Gc


def boundary_face(kind, sgn, sI, vI, zI, Gface, normal=None, mu=0.0):
    """Riemann solution (s*, v*) at a boundary face from the interior state.

    sgn = +1 for the face after the last cell (interior on its left), -1 for
    the face before the first. One-sided relation s* = sI + sgn zI (v* - vI).
    normal=None: the normal system (s is s_hat, Gface the weight at the face).
    Otherwise the tangential system with the actual normal face stress
    ``normal`` for the friction limit."""
    def v_of(s):
        return vI + sgn * (s - sI) / zI
    if normal is None:
        if kind == "free":
            s = -Gface + 0.0 * sI
            return s, v_of(s)
        s = sI - sgn * zI * vI                       # v* = 0
        sep = s + Gface > 0.0                        # would pull on the wall
        s = np.where(sep, -Gface, s)
        return s, np.where(sep, v_of(s), 0.0)
    if kind in ("free", "smooth"):
        s = np.zeros_like(sI)
        return s, v_of(s)
    stick = sI - sgn * zI * vI
    tmax = mu * np.maximum(0.0, -normal)
    slip = np.abs(stick) > tmax
    s = np.where(slip, np.sign(stick) * tmax, stick)
    return s, np.where(slip, v_of(s), 0.0)


# ----------------------------------------------------------------------------
# One sweep
# ----------------------------------------------------------------------------

def silo_sweep(f, h0, dt, mat, ends, model, params, grav=0.0, mu=0.0,
               rotation="jaumann"):
    """Lagrange step plus remap along axis 1, with boundaries and gravity.

    As lr2d_plastic.sweep_1d (f: rho, vn, vt (rows, n); s, ep (rows, n, 6,
    local Voigt order); epq) plus: ends = (kind at index 0, kind at index n),
    grav = gravity acting towards index 0 (y sweep) or 0 (x sweep), mu =
    tan(delta_w). Returns (new state, model-made eps_p_eq mass sum,
    boundary record {"lo"/"hi": (s_nn, s_nt, v_n, v_t) at the face})."""
    M, G = mat
    rows, n = f["rho"].shape
    lo, hi = ends
    h = np.full((rows, n + 2 * NG), h0)
    rho = mirror(f["rho"])
    m = rho * h
    hr = h[:, NG:-NG]
    mr = m[:, NG:-NG]
    Gf, Gc = weight_profile(f["rho"], h0, grav)

    psn_lo, pvn_lo = parity(lo, "n")
    psn_hi, pvn_hi = parity(hi, "n")
    pst_lo, pvt_lo = parity(lo, "t")
    pst_hi, pvt_hi = parity(hi, "t")
    snn = pad_ends(f["s"][..., 0] - Gc, psn_lo, psn_hi, -Gf[:, 0], -Gf[:, -1])
    vn = pad_ends(f["vn"], pvn_lo, pvn_hi)
    snt = pad_ends(f["s"][..., 3], pst_lo, pst_hi)
    vt = pad_ends(f["vt"], pvt_lo, pvt_hi)
    zp = np.sqrt(M * rho)
    zs = np.sqrt(G * rho)

    def characteristic(s, v, z):
        zc, hc = z[:, 1:-1], h[:, 1:-1]
        gs_l, gs_r = gradients(s, h)
        gv_l, gv_r = gradients(v, h)
        wr = limiter(gs_l - zc * gv_l, gs_r - zc * gv_r)
        wl = limiter(gs_l + zc * gv_l, gs_r + zc * gv_r)
        return 0.5 * (wr + wl) * hc, 0.5 * (wl - wr) / zc * hc

    def riemann(s, v, z, ds, dv, mod):
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
        return v_star, s_star, (sL, vL, zL, sR, vR, zR)

    ds_n, dv_n = characteristic(snn, vn, zp)
    vn_star, snn_star, fn = riemann(snn, vn, zp, ds_n, dv_n, M)
    ds_t, dv_t = characteristic(snt, vt, zs)
    vt_star, snt_star, ft = riemann(snt, vt, zs, ds_t, dv_t, G)

    rec = {}
    for side, kind, j, sgn, k in (("lo", lo, 0, -1, 0), ("hi", hi, -1, 1, -1)):
        # interior state at the face: right state (index 0) or left state
        o = 3 if sgn < 0 else 0
        s_b, v_b = boundary_face(kind, sgn, fn[o][:, k], fn[o + 1][:, k], fn[o + 2][:, k],
                                 Gf[:, j])
        snn_star[:, j], vn_star[:, j] = s_b, v_b
        normal = s_b + Gf[:, j]
        t_b, w_b = boundary_face(kind, sgn, ft[o][:, k], ft[o + 1][:, k], ft[o + 2][:, k],
                                 None, normal=normal, mu=mu)
        snt_star[:, j], vt_star[:, j] = t_b, w_b
        rec[side] = (normal, t_b.copy(), v_b.copy(), w_b.copy())

    vn_new = f["vn"] + dt * (snn_star[:, 1:] - snn_star[:, :-1]) / mr
    vt_new = f["vt"] + dt * (snt_star[:, 1:] - snt_star[:, :-1]) / mr

    d = dt * vn_star
    h_new = hr + d[:, 1:] - d[:, :-1]
    hmid = 0.5 * (hr + h_new)
    gam = (dt * (vt_star[:, 1:] - vt_star[:, :-1])) / hmid
    deps = np.zeros((rows, n, 6))
    deps[..., 0] = (h_new - hr) / hmid
    deps[..., 3] = gam
    s_old, ep_old = f["s"], f["ep"]
    if rotation == "jaumann":
        a = 0.5 * gam
        s_old = rotate_inplane(s_old, a)
        ep_old = rotate_inplane(ep_old, a, 0.5)
    act = (deps[..., 0] != 0.0) | (gam != 0.0)
    s_new, ep_new, epq_new = s_old.copy(), ep_old.copy(), f["epq"].copy()
    if act.any():
        s_new[act], ep_new[act], epq_new[act], _ = model(
            deps[act], s_old[act], ep_old[act], f["epq"][act], params)
    made = (mr * (epq_new - f["epq"])).sum()
    rho_new = rho[:, NG:-NG] - rho[:, NG:-NG] * (h_new - hr) / h_new

    # --- remap, no mass flux through the boundary faces ---------------------
    hp = mirror(h_new)
    rp = mirror(rho_new)
    mp = rp * hp
    mass_flux = d * swept_average(rp, hp, d)
    mass_flux[:, 0] = 0.0
    mass_flux[:, -1] = 0.0
    out = {}
    dh = h_new - h0
    dmf = mass_flux[:, :-1] - mass_flux[:, 1:]
    out["rho"] = rho_new + (rho_new * dh + dmf) / h0
    m_new = out["rho"] * h0

    def carry(val, par=(1, 1), lim=limiter):
        if val.ndim == 3:
            pv = mirror(val)
        else:
            pv = pad_ends(val, *par)
        flux = mass_flux_of(pv, mp, mass_flux, lim)
        if val.ndim == 3:
            return val + (flux[:, :-1] - flux[:, 1:] - val * dmf[..., None]) / m_new[..., None]
        return val + (flux[:, :-1] - flux[:, 1:] - val * dmf) / m_new

    out["vn"] = carry(vn_new, (pvn_lo, pvn_hi))
    out["vt"] = carry(vt_new, (pvt_lo, pvt_hi))
    out["s"] = carry(s_new)
    out["ep"] = carry(ep_new, lim=superbee)
    out["epq"] = carry(epq_new, lim=superbee)
    return out, made, rec


# ----------------------------------------------------------------------------
# Tank
# ----------------------------------------------------------------------------

SIDES = ("left", "right", "bottom", "top")


class SiloBox(lp.PlasticBox2D):
    """Closed tank [0, lx] x [0, ly], gravity -y, boundaries per side.

    walls: dict side -> "smooth" | "friction" | "free" (default smooth walls,
    free top). delta_w: wall friction angle (degrees, all frictional sides).
    t_ramp: gravity ramp time (0: full gravity at once). damping: alpha (1/s)."""

    def __init__(self, nx, ny, lx, ly, rho, model, params, walls=None, delta_w=0.0,
                 g=9.81, t_ramp=0.0, damping=0.0):
        super().__init__(nx, ny, lx, ly, rho, model, params, bc=("wall", "wall"))
        self.walls = {"left": "smooth", "right": "smooth", "bottom": "smooth", "top": "free"}
        self.walls.update(walls or {})
        for k, v in self.walls.items():
            if k not in SIDES or v not in KINDS:
                raise ValueError(f"bad boundary {k}: {v}")
        self.delta_w = delta_w
        self.mu = math.tan(math.radians(delta_w))
        self.g = g
        self.t_ramp = t_ramp
        self.damping = damping
        self.C = models.compliance_matrix(self.K, self.G)
        self.face = {}                 # side -> (s_nn, s_nt, v_n, v_t) at the faces
        self.slip = {"left": np.zeros(ny), "right": np.zeros(ny),
                     "bottom": np.zeros(nx), "top": np.zeros(nx)}
        self._g_now = g

    def g_at(self, t):
        if self.t_ramp <= 0.0 or t >= self.t_ramp:
            return self.g
        return self.g * 0.5 * (1.0 - math.cos(math.pi * t / self.t_ramp))

    def sweep(self, st, dt, axis):
        mat = (self.M, self.G)
        if axis == 0:
            fl = {"rho": st["rho"], "vn": st["vx"], "vt": st["vy"], "s": st["s"],
                  "ep": st["ep"], "epq": st["epq"]}
            g, made, rec = silo_sweep(fl, self.dx, dt, mat,
                                      (self.walls["left"], self.walls["right"]),
                                      self.model, self.params, 0.0, self.mu, self.rotation)
            self._rec.append(("left", rec["lo"]))
            self._rec.append(("right", rec["hi"]))
            return {"rho": g["rho"], "vx": g["vn"], "vy": g["vt"], "s": g["s"],
                    "ep": g["ep"], "epq": g["epq"]}, made * self.dy
        tr = lambda a: np.swapaxes(a, 0, 1)                      # noqa: E731
        fl = {"rho": tr(st["rho"]), "vn": tr(st["vy"]), "vt": tr(st["vx"]),
              "s": tr(st["s"][..., PERM_Y]), "ep": tr(st["ep"][..., PERM_Y]),
              "epq": tr(st["epq"])}
        g, made, rec = silo_sweep(fl, self.dy, dt, mat,
                                  (self.walls["bottom"], self.walls["top"]),
                                  self.model, self.params, self._g_now, self.mu,
                                  self.rotation)
        self._rec.append(("bottom", rec["lo"]))
        self._rec.append(("top", rec["hi"]))
        return {"rho": tr(g["rho"]), "vx": tr(g["vt"]), "vy": tr(g["vn"]),
                "s": tr(g["s"])[..., PERM_Y], "ep": tr(g["ep"])[..., PERM_Y],
                "epq": tr(g["epq"])}, made * self.dx

    def step(self, dt):
        self._g_now = self.g_at(self.t + 0.5 * dt)
        self._rec = []
        super().step(dt)
        if self.damping > 0.0:
            fac = 1.0 / (1.0 + self.damping * dt)
            self.vx = self.vx * fac
            self.vy = self.vy * fac
        acc = {}
        for side, r in self._rec:
            acc.setdefault(side, []).append(r)
        for side, rs in acc.items():
            self.face[side] = tuple(np.mean([r[i] for r in rs], axis=0) for i in range(4))
            self.slip[side] += np.abs(self.face[side][3]) * dt

    # --- set-up and diagnostics ---------------------------------------------

    def set_k0_state(self, K0):
        """Exact discrete equilibrium under full gravity: s_yy from the weight
        profile the y sweep builds (s_yy = 0 at the top face), s_xx = s_zz =
        K0 s_yy, velocities zero."""
        Gf, Gc = weight_profile(self.rho.T, self.dy, self.g)
        syy = (Gc - Gf[:, -1:]).T
        self.s = np.zeros((self.ny, self.nx, 6))
        self.s[..., 1] = syy
        self.s[..., 0] = K0 * syy
        self.s[..., 2] = K0 * syy
        self.vx = np.zeros((self.ny, self.nx))
        self.vy = np.zeros((self.ny, self.nx))

    def kinetic_energy(self):
        return 0.5 * (self.cell_mass() * (self.vx ** 2 + self.vy ** 2)).sum()

    def strain_energy(self):
        s = self.s
        return 0.5 * ((s @ self.C.T) * s).sum() * self.dx * self.dy

    def weight(self):
        return self.total_mass() * self._g_now

    def forces(self):
        """Upward forces on the material per metre depth (N/m): bottom
        (-s_yy at the bottom faces), each side wall (shear traction)."""
        fb = -(self.face["bottom"][0]).sum() * self.dx
        fl = -(self.face["left"][1]).sum() * self.dy
        fr = (self.face["right"][1]).sum() * self.dy
        return fb, fl, fr


# ----------------------------------------------------------------------------
# Material, geometry, theory
# ----------------------------------------------------------------------------

RHO = 1800.0
T0_MC = dict(E=10.0e6, nu=0.3, phi=35.0, c=0.0, psi=0.0, sigma_t=0.0)
G_ACC = 9.81
W_TANK = 0.40
DELTA_W = 20.0


def janssen(z, W, K, mu, rho=RHO, g=G_ACC):
    """Mean vertical stress magnitude (Pa) at depth z in a plane bin of width
    W (two walls): W dsv/dz = rho g W - 2 mu K sv, sv(0) = 0."""
    a = 2.0 * mu * K / W
    return rho * g / a * (1.0 - np.exp(-a * z))


def k0_elastic(nu):
    return nu / (1.0 - nu)


def settle(box, cfl=0.5, ke_tol=1e-6, t_max=1.0, every=20, log=None):
    """Run until the ramp is over and KE <= ke_tol * strain energy (checked
    every `every` steps). Returns (settled, wall seconds, history list)."""
    dt = box.stable_dt(cfl)
    hist = []
    t0 = time.perf_counter()
    while box.t < t_max:
        box.step(dt)
        if box.nstep % every == 0:
            ke, u = box.kinetic_energy(), box.strain_energy()
            hist.append((box.t, ke, u))
            if log:
                log(box, ke, u)
            if box.t >= box.t_ramp and u > 0.0 and ke <= ke_tol * u:
                return True, time.perf_counter() - t0, hist
    return False, time.perf_counter() - t0, hist


def wall_extrapolate(q, side):
    """Linear extrapolation of a cell field (ny, nx) to the left or right wall."""
    if side == "left":
        return 1.5 * q[:, 0] - 0.5 * q[:, 1]
    return 1.5 * q[:, -1] - 0.5 * q[:, -2]


def continue_undamped(box, frac=0.1, cfl=0.5):
    """Run a further `frac` of the steps taken so far with damping off.
    Returns (largest stress change / rho g H, KE / strain energy at the end,
    largest speed during the continuation)."""
    s0 = box.s.copy()
    alpha, box.damping = box.damping, 0.0
    dt = box.stable_dt(cfl)
    n = max(1, int(math.ceil(frac * box.nstep)))
    vmax = 0.0
    for _ in range(n):
        box.step(dt)
        vmax = max(vmax, np.abs(box.vx).max(), np.abs(box.vy).max())
    box.damping = alpha
    scale = RHO * G_ACC * box.ly
    return (np.abs(box.s - s0).max() / scale,
            box.kinetic_energy() / box.strain_energy(), vmax, n)


# ----------------------------------------------------------------------------
# Cases
# ----------------------------------------------------------------------------

def make_tank(nx, ny, W, H, friction, t_ramp, damping):
    walls = {"left": "smooth", "right": "smooth", "bottom": "smooth", "top": "free"}
    if friction:
        walls.update(left="friction", right="friction", bottom="friction")
    return SiloBox(nx, ny, W, H, RHO, "MohrCoulomb", T0_MC, walls=walls,
                   delta_w=DELTA_W if friction else 0.0, t_ramp=t_ramp, damping=damping)


def lake_at_rest(nx=10, ny=120, W=0.1, H=1.2, nsteps=1000):
    box = make_tank(nx, ny, W, H, False, 0.0, 0.0)
    box.set_k0_state(k0_elastic(T0_MC["nu"]))
    s0 = box.s.copy()
    dt = box.stable_dt(0.5)
    vmax = 0.0
    for _ in range(nsteps):
        box.step(dt)
        vmax = max(vmax, np.abs(box.vx).max(), np.abs(box.vy).max())
    return vmax, np.abs(box.s - s0).max() / (RHO * G_ACC * H), dt


def settled_case(friction, nx, ny, W, H, t_ramp, damping, ke_tol=1e-6, quiet=True):
    box = make_tank(nx, ny, W, H, friction, t_ramp, damping)
    m0 = box.total_mass()
    ok, sec, hist = settle(box, ke_tol=ke_tol)
    res = dict(box=box, settled=ok, seconds=sec, steps=box.nstep, t=box.t,
               dt=box.stable_dt(0.5), dm=abs(box.total_mass() - m0) / m0)
    fb, fl, fr = box.forces()
    res["forces"] = (fb, fl, fr, box.weight())
    res["snapshot"] = box.state()
    res["face"] = {k: tuple(a.copy() for a in v) for k, v in box.face.items()}
    t1 = time.perf_counter()
    res["cont"] = continue_undamped(box)
    res["seconds_cont"] = time.perf_counter() - t1
    return res


def janssen_fit(box, face):
    """Mean vertical stress against depth, the K measured at the walls and
    the K that fits best (least largest relative error), below z = W."""
    W, H, mu = box.lx, box.ly, box.mu
    s = box.s
    z = H - box.yc
    sv = -s[..., 1].mean(axis=1)
    sh = 0.5 * (-wall_extrapolate(s[..., 0], "left") - wall_extrapolate(s[..., 0], "right"))
    deep = z >= W
    K_wall = sh[deep].sum() / sv[deep].sum()
    err = lambda K: np.abs(janssen(z, W, K, mu)[deep] / sv[deep] - 1.0)   # noqa: E731
    Ks = np.linspace(0.2, 0.8, 601)
    K_fit = Ks[np.argmin([err(K).max() for K in Ks])]
    return dict(z=z, sv=sv, deep=deep, K_wall=K_wall, e_wall=err(K_wall).max(),
                K_fit=K_fit, e_fit=err(K_fit).max(),
                hydro=sv / (RHO * G_ACC * z))


# ----------------------------------------------------------------------------
# Acceptance checks for WP 3.1p
# ----------------------------------------------------------------------------

# Settling controls used by the checks (decision 0020): gravity ramp 0.1 s
# (about two periods of the column's first mode, pi c / (2 H) = 113 rad/s for
# H = 1.2 m) and damping alpha = 250 1/s (about critical for that mode), so the
# column loads monotonically and the walls stay at the friction limit.
T_RAMP = 0.1
ALPHA = 250.0
KE_TOL = 1e-6
CELL = 0.01                     # 10 mm cells (the spec's 5 mm is the production mesh)


def report_settle(r, label):
    print(f"  {label}: {r['steps']} steps of dt = {1e6 * r['dt']:.1f} us to t = "
          f"{r['t']:.3f} s, {r['seconds']:.0f} s wall clock; settled (KE <= "
          f"{KE_TOL:g} strain energy): {r['settled']}")
    ds, ke, vmax, n = r["cont"]
    good = r["settled"] and ds <= 0.02 and ke <= 10.0 * KE_TOL
    print(f"  {n} more steps without damping: largest stress change {100 * ds:.2f} % of "
          f"rho g H (<= 2 %), KE / strain energy {ke:.1e} (<= {10 * KE_TOL:g}), largest "
          f"speed {vmax:.1e} m/s; mass change {r['dm']:.0e}: {'PASS' if good else 'FAIL'}")
    return good and r["dm"] <= 1e-12


def verify():
    ok = True
    t_all = time.perf_counter()
    print("[1] Vectorised Mohr-Coulomb (constitutive/batch.py) against models.py,")
    print("    200 random paths of 30 increments")
    for name, (d, n) in batch.check_batch_mc().items():
        good = d <= 1e-12 and n["plastic"] > 0 and n["edge"] > 0 and n["apex"] > 0
        ok &= good
        print(f"  {name:28s}: largest difference {d:.1e} (<= 1e-12); returns: "
              f"{n['plastic']} plastic, {n['edge']} edge, {n['apex']} apex: "
              f"{'PASS' if good else 'FAIL'}")

    print("\n[2] Lake at rest: K0 column (K0 = nu / (1 - nu)) between smooth walls,")
    print("    0.1 m x 0.6 m, 10 mm cells, full gravity, no damping, 1000 steps")
    vmax, ds, dt = lake_at_rest(nx=10, ny=60, W=0.1, H=0.6)
    good = vmax <= 1e-10
    ok &= good
    print(f"  largest speed {vmax:.1e} m/s (<= 1e-10), largest stress change "
          f"{ds:.1e} of rho g H: {'PASS' if good else 'FAIL'}")

    H = 3.0 * W_TANK
    nyc = int(round(H / CELL))
    print(f"\n[3] Smooth walls, settled from zero stress: column 0.1 m x {H:.1f} m, "
          f"{int(round(0.1 / CELL))} x {nyc} cells")
    r = settled_case(False, int(round(0.1 / CELL)), nyc, 0.1, H, T_RAMP, ALPHA, KE_TOL)
    ok &= report_settle(r, "settling")
    box = r["box"]
    box.set_state(r["snapshot"])
    weight_per_width = box.total_mass() * G_ACC / box.lx           # rho g H, exact mass
    sb = -r["face"]["bottom"][0].mean()
    e_b = abs(sb / weight_per_width - 1.0)
    K0 = k0_elastic(T0_MC["nu"])
    z = H - box.yc
    below = z >= 0.1 * H
    kx = box.s[below, :, 0] / box.s[below, :, 1]
    kz = box.s[below, :, 2] / box.s[below, :, 1]
    e_k = np.abs(kx / K0 - 1.0).max()
    good = e_b <= 0.02 and e_k <= 0.05
    ok &= good
    print(f"  bottom s_yy = {sb:.1f} Pa against rho g H = {weight_per_width:.1f} Pa: "
          f"{100 * e_b:.3f} % (<= 2 %)")
    print(f"  expected s_xx / s_yy = K0 = nu / (1 - nu) = {K0:.4f} (uniaxial strain, "
          f"elastic: inside the Mohr-Coulomb cone, active limit "
          f"{(1 - math.sin(math.radians(35))) / (1 + math.sin(math.radians(35))):.3f})")
    print(f"  s_xx / s_yy from {kx.min():.4f} to {kx.max():.4f} below the top 10 %: "
          f"{100 * e_k:.2f} % from K0 (<= 5 %); s_zz / s_yy {kz.min():.4f} to "
          f"{kz.max():.4f}; plastic strain max {box.epq.max():.1e}: "
          f"{'PASS' if good else 'FAIL'}")

    nx = int(round(W_TANK / CELL))
    print(f"\n[4] Frictional walls and bottom (delta_w {DELTA_W:g} deg), W = {W_TANK} m, "
          f"H = 3 W, {nx} x {nyc} cells")
    r = settled_case(True, nx, nyc, W_TANK, H, T_RAMP, ALPHA, KE_TOL)
    ok &= report_settle(r, "settling")
    box = r["box"]
    box.set_state(r["snapshot"])
    j = janssen_fit(box, r["face"])
    # WP 3.0 (decision 0019): the fitted K must lie between the wall-yield and
    # the elastic value, tested as 0.27 to 0.45
    good = j["e_wall"] <= 0.10 and 0.27 <= j["K_wall"] <= 0.45 and 0.27 <= j["K_fit"] <= 0.45
    ok &= good
    print(f"  K measured at the walls (wall s_xx / mean s_yy, depth >= W): "
          f"{j['K_wall']:.4f}; best-fitting K {j['K_fit']:.4f} (error {100 * j['e_fit']:.2f} %); both within 0.27 to 0.45 (decision 0019)")
    print(f"  mean vertical stress against Janssen with the measured K, depth >= W: "
          f"largest error {100 * j['e_wall']:.2f} % (<= 10 %): {'PASS' if good else 'FAIL'}")
    sv, zz = j["sv"], j["z"]
    for frac in (0.25, 0.5, 0.75, 1.0):
        i = int(np.argmin(np.abs(zz - frac * H)))
        print(f"    depth {zz[i]:.3f} m: s_v {sv[i]:8.1f} Pa, Janssen "
              f"{janssen(zz[i], W_TANK, j['K_wall'], box.mu):8.1f} Pa, hydrostatic "
              f"{RHO * G_ACC * zz[i]:8.1f} Pa")
    mu = box.mu
    worst = 0.0
    for side in ("left", "right"):
        ratio = np.abs(wall_extrapolate(box.s[..., 3], side)
                       / wall_extrapolate(box.s[..., 0], side)) / mu
        sl = box.slip[side]
        zone = (sl >= 0.1 * sl.max()) & j["deep"]
        e = np.abs(ratio[zone] - 1.0).max()
        worst = max(worst, e)
        print(f"  {side} wall: slips (slip >= 10 % of the largest) on "
              f"{int((sl >= 0.1 * sl.max()).sum())} of {box.ny} faces; shear / normal "
              f"stress (cells extrapolated to the wall) / tan(delta_w) from "
              f"{ratio[zone].min():.4f} to {ratio[zone].max():.4f} where it slips, "
              f"depth >= W")
    good = worst <= 0.05
    ok &= good
    print(f"  wall ratio within {100 * worst:.2f} % of tan(delta_w) (<= 5 %): "
          f"{'PASS' if good else 'FAIL'}")
    fb, fl, fr, wt = r["forces"]
    e_f = abs(fb + fl + fr - wt) / wt
    good = e_f <= 0.01
    ok &= good
    print(f"  weight {wt:.1f} N/m = bottom {fb:.1f} + walls {fl:.1f} + {fr:.1f} N/m, "
          f"imbalance {100 * e_f:.3f} % (<= 1 %): {'PASS' if good else 'FAIL'}")

    print(f"\nWP 3.1p checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--verify", action="store_true", help="run the WP 3.1p checks")
    args = p.parse_args(argv)
    if args.verify:
        return 0 if verify() else 1
    p.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
