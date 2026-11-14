#!/usr/bin/env python3
"""WP 1.3b: 2D plane-strain Lagrange-plus-remap with rotation and plasticity.

Extends the WP 1.3a scheme (``lr2d.py``: dimensionally split, cell-centred
MUSCL-Hancock Godunov Lagrange step plus MUSCL remap, symmetric splitting)
with what a solid needs beyond small elastic waves:

1. **Objective stress update (Jaumann rate, Hughes-Winget increment).** In a
   sweep along n (tangent t) the sweep knows dv_n/dn and dv_t/dn. The strain
   increment is (d eps_nn, d gamma_nt) = ((h_new - h)/h_mid,
   dt (v_t*_R - v_t*_L)/h_mid) as in WP 1.3a, and the spin increment is
   dW_tn = d gamma_nt / 2 (the n-derivative half of W = skew(grad v)). The
   stress is first rotated with Q = (I - dW/2)^-1 (I + dW/2), an exact
   rotation (Hughes and Winget 1980), then the model adds the response to the
   strain increment. The x and y sweeps add up to the full spin. Plastic
   strain (a tensor history) is rotated the same way. ``rotation="none"``
   switches the rotation off (for comparison only).
2. **Constitutive models.** Each cell keeps the 6-component stress and the
   WP 1.4 history (eps_p, eps_p_eq). A WP 1.4 model, vectorised in
   ``constitutive/batch.py`` and checked against ``constitutive/models.py``,
   maps the local strain increment (d eps_nn, 0, 0, d gamma_nt, 0, 0) to the
   new stress. The models are isotropic, so they work in the sweep frame.
   The Riemann solvers, slopes and time step use the elastic impedances
   (decision 0008).
3. **Remap as WP 1.2.** Velocity, all stress components and the history are
   remapped mass-weighted (in mass coordinates, carried by the mass flux), so
   sum(m q) is exact. Stress uses the MC limiter, history superbee
   (decision 0008). The symmetric split averages them mass-weighted, like
   velocity.
4. **Prescribed velocity (kinematic) mode.** ``Box.prescribed = f(x, y)``
   returning (vx, vy) replaces the Riemann solvers by the given face
   velocities and skips the momentum update; the cell velocities are set from
   f after each step. Used for the rigid-rotation and simple-shear checks.

Fields are arrays (ny, nx) (stress and plastic strain (ny, nx, 6)), as in
lr2d.py. SI units, stress positive in tension, thickness 1 m.

``python3 lr2d_plastic.py --verify`` runs the WP 1.3b checks
(tests/test_proto_2d_rotation.sh). ``--study`` runs the formulation-decision
measurements recorded in docs/decisions/0013-eulerian-formulation.md.
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

from lr1d import limiter  # noqa: E402
from lr1d_plastic import superbee  # noqa: E402
import batch  # noqa: E402
import models  # noqa: E402

NG = 2
# Local Voigt order in a sweep along n: (nn, tt, zz, nt, tz, zn).
# x sweep: n = x, t = y, the global order. y sweep: n = y, t = x.
PERM_Y = [1, 0, 2, 3, 5, 4]


# ----------------------------------------------------------------------------
# Helpers on (rows, n[, k]) arrays, sweep along axis 1
# ----------------------------------------------------------------------------

def pad(q, bc, odd=False):
    width = [(0, 0), (NG, NG)] + [(0, 0)] * (q.ndim - 2)
    if bc == "periodic":
        return np.pad(q, width, mode="wrap")
    p = np.pad(q, width, mode="symmetric")
    if odd:
        p[:, :NG] *= -1.0
        p[:, -NG:] *= -1.0
    return p


def _expand(w, q):
    return w if w.ndim == q.ndim else w[..., None]


def gradients(q, h):
    h = _expand(h, q)
    gl = (q[:, 1:-1] - q[:, :-2]) / (0.5 * (h[:, 1:-1] + h[:, :-2]))
    gr = (q[:, 2:] - q[:, 1:-1]) / (0.5 * (h[:, 2:] + h[:, 1:-1]))
    return gl, gr


def swept_average(q, h, d):
    """Volume-weighted swept-region average (lr2d.swept_average)."""
    gl, gr = gradients(q, h)
    s = limiter(gl, gr)
    qc, hc = q[:, 1:-1], h[:, 1:-1]
    a = np.abs(d)
    from_left = qc[:, :-1] + s[:, :-1] * 0.5 * (hc[:, :-1] - a)
    from_right = qc[:, 1:] - s[:, 1:] * 0.5 * (hc[:, 1:] - a)
    return np.where(d > 0.0, from_left, from_right)


def mass_flux_of(q, m, mass_flux, lim=limiter):
    """Flux of m q through the faces, q reconstructed in mass coordinates
    (lr2d.mass_flux_of with a choice of limiter; q may have a last axis)."""
    gl, gr = gradients(q, m)
    s = lim(gl, gr)
    qc = q[:, 1:-1]
    mc = _expand(m[:, 1:-1], q)
    mf = _expand(mass_flux, q)
    am = np.abs(mf)
    from_left = qc[:, :-1] + s[:, :-1] * 0.5 * (mc[:, :-1] - am)
    from_right = qc[:, 1:] - s[:, 1:] * 0.5 * (mc[:, 1:] - am)
    return mf * np.where(mf > 0.0, from_left, from_right)


def rotate_inplane(t, a, shear_factor=1.0):
    """Rotate the in-plane part of Voigt arrays t (rows, n, 6) by the
    Hughes-Winget rotation for the spin increment dW = [[0, -a], [a, 0]]
    (in the local (n, t) frame). shear_factor 0.5 for engineering shear
    strains. a == 0 leaves t unchanged bit-for-bit."""
    den = 1.0 + 0.25 * a * a
    c = (1.0 - 0.25 * a * a) / den
    s = a / den
    nn, tt = t[..., 0], t[..., 1]
    nt = t[..., 3] * shear_factor
    out = t.copy()
    out[..., 0] = c * c * nn - 2.0 * c * s * nt + s * s * tt
    out[..., 1] = s * s * nn + 2.0 * c * s * nt + c * c * tt
    out[..., 3] = (c * s * (nn - tt) + (c * c - s * s) * nt) / shear_factor
    # Out-of-plane shears (tz, zn) rotate as a vector; zero in plane strain.
    tz, zn = t[..., 4], t[..., 5]
    out[..., 5] = c * zn - s * tz
    out[..., 4] = s * zn + c * tz
    return out


# ----------------------------------------------------------------------------
# One sweep
# ----------------------------------------------------------------------------

def sweep_1d(f, h0, dt, mat, bc, model, params, rotation="jaumann", kin=None):
    """Lagrange step plus remap along axis 1 of (rows, n) arrays.

    f: dict rho, vn, vt (rows, n), s, ep (rows, n, 6, local Voigt order),
    epq (rows, n). mat: (M, G) elastic. kin: None, or (vn_star, vt_star) face
    velocities (rows, n+1) for the kinematic mode. Returns (new dict,
    model-made eps_p_eq mass sum)."""
    M, G = mat[:2]
    Gz = mat[2] if len(mat) > 2 else G     # modulus for the S impedance (study)
    rows, n = f["rho"].shape
    h = np.full((rows, n + 2 * NG), h0)
    rho = pad(f["rho"], bc)
    m = rho * h
    hr = h[:, NG:-NG]
    mr = m[:, NG:-NG]

    if kin is None:
        vn = pad(f["vn"], bc, odd=True)
        vt = pad(f["vt"], bc, odd=True)
        snn = pad(f["s"][..., 0], bc)
        snt = pad(f["s"][..., 3], bc)
        zp = np.sqrt(M * rho)
        zs = np.sqrt(Gz * rho)

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
            return v_star, s_star

        ds_n, dv_n = characteristic(snn, vn, zp)
        vn_star, snn_star = riemann(snn, vn, zp, ds_n, dv_n, M)
        ds_t, dv_t = characteristic(snt, vt, zs)
        vt_star, snt_star = riemann(snt, vt, zs, ds_t, dv_t, G)
        vn_new = f["vn"] + dt * (snn_star[:, 1:] - snn_star[:, :-1]) / mr
        vt_new = f["vt"] + dt * (snt_star[:, 1:] - snt_star[:, :-1]) / mr
    else:
        vn_star, vt_star = kin
        vn_new, vt_new = f["vn"], f["vt"]

    d = dt * vn_star
    h_new = hr + d[:, 1:] - d[:, :-1]
    hmid = 0.5 * (hr + h_new)
    gam = (dt * (vt_star[:, 1:] - vt_star[:, :-1])) / hmid
    # Strain increment in the local frame (same arithmetic as lr1d_plastic)
    deps = np.zeros((rows, n, 6))
    deps[..., 0] = (h_new - hr) / hmid
    deps[..., 3] = gam
    s_old, ep_old = f["s"], f["ep"]
    if rotation == "jaumann":
        a = 0.5 * gam                         # spin increment dW_tn
        s_old = rotate_inplane(s_old, a)
        ep_old = rotate_inplane(ep_old, a, 0.5)
    # A cell with no strain increment in this sweep keeps its state: the
    # model is not called. (A remapped stress can sit slightly outside the
    # yield surface; it is returned at the cell's next strain increment, as
    # in 1D. Without this, a sweep across a direction where nothing moves
    # would return it early, and a 1D problem on the 2D mesh would differ
    # from lr1d_plastic.py.)
    act = (deps[..., 0] != 0.0) | (gam != 0.0)
    s_new, ep_new, epq_new = s_old.copy(), ep_old.copy(), f["epq"].copy()
    if act.any():
        s_new[act], ep_new[act], epq_new[act], _ = model(
            deps[act], s_old[act], ep_old[act], f["epq"][act], params)
    made = (mr * (epq_new - f["epq"])).sum()
    rho_new = rho[:, NG:-NG] - rho[:, NG:-NG] * (h_new - hr) / h_new

    # --- remap (increments of the Lagrangian values, as lr2d) ---------------
    hp = pad(h_new, bc)
    rp = pad(rho_new, bc)
    mp = rp * hp
    mass_flux = d * swept_average(rp, hp, d)
    out = {}
    dh = h_new - h0
    dmf = mass_flux[:, :-1] - mass_flux[:, 1:]
    out["rho"] = rho_new + (rho_new * dh + dmf) / h0
    m_new = out["rho"] * h0

    def carry(val, odd=False, lim=limiter):
        flux = mass_flux_of(pad(val, bc, odd), mp, mass_flux, lim)
        if val.ndim == 3:
            return val + (flux[:, :-1] - flux[:, 1:] - val * dmf[..., None]) / m_new[..., None]
        return val + (flux[:, :-1] - flux[:, 1:] - val * dmf) / m_new

    if kin is None:
        out["vn"] = carry(vn_new, True)
        out["vt"] = carry(vt_new, True)
    else:
        out["vn"], out["vt"] = vn_new, vt_new
    out["s"] = carry(s_new)
    out["ep"] = carry(ep_new, lim=superbee)
    out["epq"] = carry(epq_new, lim=superbee)
    return out, made


# ----------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------

class PlasticBox2D:
    """Plane-strain box with a WP 1.4 model, split Lagrange plus remap."""

    def __init__(self, nx, ny, lx, ly, rho, model, params,
                 bc=("periodic", "periodic"), rotation="jaumann", split="symmetric"):
        self.nx, self.ny = nx, ny
        self.lx, self.ly = lx, ly
        self.dx, self.dy = lx / nx, ly / ny
        self.params = params
        self.model = batch.MODELS[model] if isinstance(model, str) else model
        self.K, self.G = models.bulk_shear(params)
        self.M = self.K + 4.0 * self.G / 3.0
        self.rho0 = rho
        self.cp = math.sqrt(self.M / rho)
        self.bc = bc
        self.rotation = rotation
        self.split = split
        self.s_impedance = "shear"   # "p": S system with the P impedance (study)
        self.prescribed = None
        self.xc = (np.arange(nx) + 0.5) * self.dx
        self.yc = (np.arange(ny) + 0.5) * self.dy
        self.rho = np.full((ny, nx), rho)
        self.vx = np.zeros((ny, nx))
        self.vy = np.zeros((ny, nx))
        self.s = np.zeros((ny, nx, 6))
        self.ep = np.zeros((ny, nx, 6))
        self.epq = np.zeros((ny, nx))
        self.t = 0.0
        self.nstep = 0
        self.plastic_made = 0.0

    KEYS = ("rho", "vx", "vy", "s", "ep", "epq")

    def state(self):
        return {k: getattr(self, k).copy() for k in self.KEYS}

    def set_state(self, st):
        for k in self.KEYS:
            setattr(self, k, st[k])

    def set_prescribed(self, f):
        """Kinematic mode: f(x, y) -> (vx, vy) on arrays."""
        self.prescribed = f
        X, Y = np.meshgrid(self.xc, self.yc)
        self.vx, self.vy = f(X, Y)

    def _kin(self, axis):
        f = self.prescribed
        if f is None:
            return None
        if axis == 0:
            X, Y = np.meshgrid(np.arange(self.nx + 1) * self.dx, self.yc)
            vx, vy = f(X, Y)
            return vx, vy
        X, Y = np.meshgrid(self.xc, np.arange(self.ny + 1) * self.dy)
        vx, vy = f(X, Y)
        return vy.T, vx.T

    def sweep(self, st, dt, axis):
        mat = (self.M, self.G, self.M if self.s_impedance == "p" else self.G)
        kin = self._kin(axis)
        if axis == 0:
            fl = {"rho": st["rho"], "vn": st["vx"], "vt": st["vy"], "s": st["s"],
                  "ep": st["ep"], "epq": st["epq"]}
            g, made = sweep_1d(fl, self.dx, dt, mat, self.bc[0], self.model,
                               self.params, self.rotation, kin)
            made *= self.dy
            return {"rho": g["rho"], "vx": g["vn"], "vy": g["vt"], "s": g["s"],
                    "ep": g["ep"], "epq": g["epq"]}, made
        tr = lambda a: np.swapaxes(a, 0, 1)                      # noqa: E731
        fl = {"rho": tr(st["rho"]), "vn": tr(st["vy"]), "vt": tr(st["vx"]),
              "s": tr(st["s"][..., PERM_Y]), "ep": tr(st["ep"][..., PERM_Y]),
              "epq": tr(st["epq"])}
        g, made = sweep_1d(fl, self.dy, dt, mat, self.bc[1], self.model,
                           self.params, self.rotation, kin)
        made *= self.dx
        return {"rho": tr(g["rho"]), "vx": tr(g["vt"]), "vy": tr(g["vn"]),
                "s": tr(g["s"])[..., PERM_Y], "ep": tr(g["ep"])[..., PERM_Y],
                "epq": tr(g["epq"])}, made

    def cell_mass(self):
        return self.rho * self.dx * self.dy

    def total_mass(self):
        return self.cell_mass().sum()

    def total_momentum(self):
        m = self.cell_mass()
        return (m * self.vx).sum(), (m * self.vy).sum()

    def stable_dt(self, cfl):
        vmax = max(np.abs(self.vx).max(), np.abs(self.vy).max())
        if self.prescribed is not None:
            return cfl * min(self.dx, self.dy) / vmax
        return cfl * min(self.dx, self.dy) / (self.cp + vmax)

    def step(self, dt):
        st = self.state()
        if self.split == "additive":
            self._step_additive(st, dt)
            return
        a, ma = self.sweep(st, dt, 0)
        a, mb = self.sweep(a, dt, 1)
        b, mc = self.sweep(st, dt, 1)
        b, md = self.sweep(b, dt, 0)
        rho = 0.5 * (a["rho"] + b["rho"])
        new = {"rho": rho}
        dr = a["rho"] - b["rho"]
        # Mass-weighted mean, exact when a == b (as lr2d.Box2D.step)
        for k in ("vx", "vy", "epq"):
            new[k] = 0.5 * (a[k] + b[k]) + 0.25 * dr * (a[k] - b[k]) / rho
        for k in ("s", "ep"):
            new[k] = 0.5 * (a[k] + b[k]) + 0.25 * (dr / rho)[..., None] * (a[k] - b[k])
        self.set_state(new)
        if self.prescribed is not None:
            X, Y = np.meshgrid(self.xc, self.yc)
            self.vx, self.vy = self.prescribed(X, Y)
        self.plastic_made += 0.5 * (ma + mb + mc + md)
        self.t += dt
        self.nstep += 1

    def _step_additive(self, st, dt):
        """Unsplit proxy (study only): both sweeps start from the same state
        and their conservative increments are added, the way a face-based
        unsplit finite-volume update uses all faces of a cell at once.
        Mass and m*q increments add, so it stays conservative."""
        a, ma = self.sweep(st, dt, 0)
        b, mb = self.sweep(st, dt, 1)
        rho = a["rho"] + b["rho"] - st["rho"]
        new = {"rho": rho}
        for k in ("vx", "vy", "epq", "s", "ep"):
            r0, ra, rb = st["rho"], a["rho"], b["rho"]
            if st[k].ndim == 3:
                r0, ra, rb, rr = r0[..., None], ra[..., None], rb[..., None], rho[..., None]
            else:
                rr = rho
            new[k] = (ra * a[k] + rb * b[k] - r0 * st[k]) / rr
        self.set_state(new)
        if self.prescribed is not None:
            X, Y = np.meshgrid(self.xc, self.yc)
            self.vx, self.vy = self.prescribed(X, Y)
        self.plastic_made += ma + mb
        self.t += dt
        self.nstep += 1

    def run(self, t_end, cfl=0.5, callback=None):
        nsteps = int(np.ceil(t_end / self.stable_dt(cfl)))
        dt = t_end / nsteps
        for k in range(nsteps):
            self.step(dt)
            if callback is not None:
                callback(self, k + 1)
        return nsteps


def principal_2d(s):
    """In-plane principal stresses (larger, smaller) and the angle (rad) of
    the larger one from x, for Voigt arrays (..., 6)."""
    sxx, syy, sxy = s[..., 0], s[..., 1], s[..., 3]
    c = 0.5 * (sxx + syy)
    r = np.sqrt((0.5 * (sxx - syy)) ** 2 + sxy ** 2)
    ang = 0.5 * np.arctan2(2.0 * sxy, sxx - syy)
    return c + r, c - r, ang


# ----------------------------------------------------------------------------
# Problems
# ----------------------------------------------------------------------------

RHO = 2000.0
VM = dict(E=20.0e6, nu=0.3, sigma_y=10.0e3, H=0.0)       # as lr1d_plastic.VM
ELASTIC = dict(E=20.0e6, nu=0.3)
DP = dict(E=20.0e6, nu=0.3, phi=30.0, c=5.0e3, psi=0.0, fit="compression")


def frame(box, nvec):
    """Normal and tangential velocity and stress for a wave normal nvec."""
    nx, ny = nvec
    tx, ty = -ny, nx
    sxx, syy, sxy = box.s[..., 0], box.s[..., 1], box.s[..., 3]
    vn = box.vx * nx + box.vy * ny
    vt = box.vx * tx + box.vy * ty
    snn = sxx * nx * nx + 2 * sxy * nx * ny + syy * ny * ny
    snt = sxx * nx * tx + sxy * (nx * ty + ny * tx) + syy * ny * ty
    stt = sxx * tx * tx + 2 * sxy * tx * ty + syy * ty * ty
    return {"vn": vn, "vt": vt, "snn": snn, "stt": stt, "szz": box.s[..., 2],
            "snt": snt, "tau": 0.5 * (snn - stt), "mean": 0.5 * (snn + stt)}


FIELDS = ("vn", "vt", "snn", "stt", "szz", "snt", "tau", "mean")


def impact(angle, n, model, params, v_frac=3.0, v0_frac=0.25, split="symmetric",
           s_impedance="shear", cfl=0.5, travel=0.22):
    """Plane impact along the normal at angle 0 or 45 degrees (periodic box).

    Material with s < 0.5 (s the distance along the normal, periodic) moves
    at +V along the normal, the rest at -V, plus a diagonal streaming
    velocity of v0_frac c_e: a compression structure runs out of s = 0.5 and
    a tension structure out of s = 0. V is v_frac times the elastic
    precursor's particle velocity of the WP 1.2 von Mises material (final
    stress about 2.6 HEL for von Mises; the same V for the other models).
    The running minimum and maximum of each field in the wave frame (normal
    and tangential velocity, s_nn, s_tt, s_zz, s_nt) are tracked over the
    start-up (until the elastic wave has run 10 cells along the normal) and
    after it. Returns a dict."""
    K, G = models.bulk_shear(params)
    M = K + 4.0 * G / 3.0
    ce = math.sqrt(M / RHO)
    Kv, Gv = models.bulk_shear(VM)
    Mv = Kv + 4.0 * Gv / 3.0
    hel = VM["sigma_y"] * Mv / (2.0 * Gv)
    V = v_frac * hel / (RHO * math.sqrt(Mv / RHO))
    if angle == 45:
        box = PlasticBox2D(n, n, 1.0, 1.0, RHO, model, params, split=split)
        nvec = (math.sqrt(0.5), math.sqrt(0.5))
        X, Y = np.meshgrid(box.xc, box.yc)
        period = math.sqrt(0.5)                # along the normal
        # Distance along the normal from integer indices: computed from the
        # cell centres in floating point, cells on one diagonal fall on
        # different sides of the threshold and the front is jagged (that was
        # most of the 45 degree overshoot reported in decision 0012).
        J, I = np.indices((n, n))
        s = ((I + J + 1) % n) * box.dx * period
        hn = box.dx / math.sqrt(2.0)           # cell size along the normal
    else:
        box = PlasticBox2D(n, 4, 1.0, 4.0 / n, RHO, model, params, split=split)
        nvec = (1.0, 0.0)
        X, Y = np.meshgrid(box.xc, box.yc)
        s = X
        hn = box.dx
        period = 1.0
    box.s_impedance = s_impedance
    sign = np.where(s < 0.5 * period, 1.0, -1.0)
    vd = v0_frac * ce / math.sqrt(2.0)
    box.vx = vd + V * sign * nvec[0]
    box.vy = vd + V * sign * nvec[1]
    m0 = box.total_mass()
    p0 = box.total_momentum()
    pscale = (box.cell_mass() * np.hypot(box.vx, box.vy)).sum()
    rng = {"start": {}, "main": {}}

    vs = (vd * (nvec[0] + nvec[1]), vd * (nvec[0] - nvec[1]))   # streaming (n, t)

    def rel(b):
        f = frame(b, nvec)
        f["vn"] = f["vn"] - vs[0]
        f["vt"] = f["vt"] - vs[1]
        return f

    def track(b, k):
        vals = rel(b)
        w = rng["start" if ce * b.t < 10 * hn else "main"]
        for key in FIELDS:
            lo, hi = w.get(key, (np.inf, -np.inf))
            w[key] = (min(lo, vals[key].min()), max(hi, vals[key].max()))

    t_end = travel * period / ce
    t0 = time.perf_counter()
    box.run(t_end, cfl=cfl, callback=track)
    cost = time.perf_counter() - t0
    final = rel(box)
    dm = abs(box.total_mass() - m0) / m0
    px, py = box.total_momentum()
    dp = max(abs(px - p0[0]), abs(py - p0[1])) / pscale
    return dict(box=box, V=V, ce=ce, rng=rng, final=final, dm=dm, dp=dp,
                cost=cost, steps=box.nstep)


def extra_extrema(r45, r0):
    """Largest amount by which a field of the 45 degree run goes beyond the
    range the same field reaches in the along-x run (1D-like, no 2D
    effects), as a fraction of the along-x range of s_nn (stresses) or of
    2V (velocities). Per window (start-up, main) and at the end."""
    out = {"start": 0.0, "main": 0.0, "end": 0.0, "field": {}}
    lo0 = {k: min(r0["rng"]["start"][k][0], r0["rng"]["main"][k][0]) for k in FIELDS}
    hi0 = {k: max(r0["rng"]["start"][k][1], r0["rng"]["main"][k][1]) for k in FIELDS}
    sscale = hi0["snn"] - lo0["snn"]
    for k in FIELDS:
        scale = 2.0 * r0["V"] if k.startswith("v") else sscale
        for win in ("start", "main"):
            lo, hi = r45["rng"][win][k]
            e = max(hi - hi0[k], lo0[k] - lo, 0.0) / scale
            out[win] = max(out[win], e)
            out["field"][(win, k)] = e
        q = r45["final"][k]
        out["end"] = max(out["end"], max(q.max() - hi0[k], lo0[k] - q.min(), 0.0) / scale)
    return out


# ----------------------------------------------------------------------------
# Acceptance checks for WP 1.3b
# ----------------------------------------------------------------------------

def rotated(s0, angle):
    """Voigt stress s0 (6,) rotated in plane by angle (rad)."""
    c, sn = math.cos(angle), math.sin(angle)
    Q = np.array([[c, -sn], [sn, c]])
    T = np.array([[s0[0], s0[3]], [s0[3], s0[1]]])
    R = Q @ T @ Q.T
    return np.array([R[0, 0], R[1, 1], s0[2], R[0, 1], 0.0, 0.0])


def disc_case(cells=100, margin=8, band=6, turns=1.0, rotation="jaumann", cfl=0.5,
              quiet=False):
    """Pre-stressed elastic disc in rigid rotation (prescribed velocity).

    Disc radius 0.5 m, `cells` cells across, centred in a periodic box with
    `margin` cells of stress-free material around it; the whole box rotates
    at omega = 1 rad/s about the disc centre. Initial stress in the disc:
    principal -40 kPa at 30 degrees from x, -100 kPa across it, s_zz -60 kPa
    (cell averages of the disc by 8 x 8 subsampling). At each quarter turn
    the stress in the interior (more than `band` cells inside the edge) is
    compared with the initial stress rotated by the turned angle, and the
    stress outside (more than `band` cells outside) with zero."""
    R = 0.5
    dx = 2.0 * R / cells
    n = cells + 2 * margin
    L = n * dx
    box = PlasticBox2D(n, n, L, L, RHO, "linearElastic", ELASTIC, rotation=rotation)
    xc = yc = 0.5 * L
    w = 1.0
    box.set_prescribed(lambda X, Y: (-w * (Y - yc), w * (X - xc)))
    sub = 8
    frac = np.zeros((n, n))
    for a in range(sub):
        for b in range(sub):
            X, Y = np.meshgrid((np.arange(n) + (a + 0.5) / sub) * dx,
                               (np.arange(n) + (b + 0.5) / sub) * dx)
            frac += (X - xc) ** 2 + (Y - yc) ** 2 < R * R
    frac /= sub * sub
    s1, s2, th = -40.0e3, -100.0e3, math.radians(30.0)
    s0 = rotated(np.array([s1, s2, -60.0e3, 0.0, 0.0, 0.0]), th)
    box.s = frac[..., None] * s0
    X, Y = np.meshgrid(box.xc, box.yc)
    r = np.hypot(X - xc, Y - yc)
    inner = r < R - band * dx
    outer = r > R + band * dx
    scale = abs(s2)
    res = []
    t0 = time.perf_counter()
    for q in range(1, int(round(4 * turns)) + 1):
        box.run(0.25 * 2.0 * math.pi / w, cfl=cfl)
        ex = rotated(s0, w * box.t)
        e_tensor = np.abs(box.s[inner] - ex).max() / scale
        p1, p2, ang = principal_2d(box.s[inner])
        e_princ = max(np.abs(p1 - s1).max(), np.abs(p2 - s2).max()) / scale
        da = (ang - (th + w * box.t) + 0.5 * math.pi) % math.pi - 0.5 * math.pi
        e_ang = math.degrees(np.abs(da).max())
        e_out = np.abs(box.s[outer]).max() / scale
        # Edge width: radial distance over which |s_zz| goes from 90 % to 10 %
        # of its disc value, along the x axis through the centre
        prof = np.abs(box.s[n // 2, :, 2]) / 60.0e3
        width = int(np.sum((prof > 0.1) & (prof < 0.9)))
        res.append(dict(turn=q / 4, steps=box.nstep, e_tensor=e_tensor, e_princ=e_princ,
                        e_ang=e_ang, e_out=e_out, width=width))
        if not quiet:
            print(f"  {q / 4:.2f} turn, {box.nstep:4d} steps: interior tensor error "
                  f"{100 * e_tensor:.4f} %, principal stresses {100 * e_princ:.4f} %, "
                  f"principal direction {e_ang:.4f} deg; outside {100 * e_out:.4f} %; "
                  f"edge cells (10-90 %, both sides) {width}")
    return res, time.perf_counter() - t0


def shear_case(model, params, gamma_end=1.0, nsteps=2000, p0=0.0):
    """Simple shear v_x = gdot (y - L/2), v_y = 0 on an 8 x 8 periodic block
    (prescribed velocity) to shear strain gamma_end. Returns gamma and the
    stress history of one cell, and the spread over all cells."""
    n = 8
    box = PlasticBox2D(n, n, 1.0, 1.0, RHO, model, params)
    gd = 1.0
    box.set_prescribed(lambda X, Y: (gd * (Y - 0.5), 0.0 * X))
    box.s[..., :3] = p0
    dt = gamma_end / gd / nsteps
    hist = np.empty((nsteps, 6))
    spread = 0.0
    for k in range(nsteps):
        box.step(dt)
        hist[k] = box.s[3, 5]
        spread = max(spread, np.abs(box.s - box.s[3, 5]).max())
    gam = gd * dt * np.arange(1, nsteps + 1)
    return gam, hist, spread, box


def largest_drop(q):
    q = np.asarray(q)
    return float(np.max(np.maximum.accumulate(q) - q))


def check_simple_shear():
    ok = True
    G = models.bulk_shear(ELASTIC)[1]
    gam, h, spread, _ = shear_case("linearElastic", ELASTIC)
    e = max(np.abs(h[:, 3] - G * np.sin(gam)).max(),
            np.abs(h[:, 0] - G * (1 - np.cos(gam))).max(),
            np.abs(h[:, 1] + G * (1 - np.cos(gam))).max()) / G
    good = e <= 1e-3
    ok &= good
    print(f"  elastic (Jaumann closed form s_xy = G sin(gamma), s_xx = -s_yy = "
          f"G (1 - cos(gamma))): largest error {e:.2e} G at gamma up to 1 "
          f"(<= 1e-3: {'PASS' if good else 'FAIL'})")
    print(f"    for the record: at gamma 1 the Jaumann elastic s_xy is "
          f"{h[-1, 3] / G:.4f} G (small-strain G gamma = 1 G)")

    k = VM["sigma_y"] / math.sqrt(3.0)
    gam, h, spread, box = shear_case("vonMises", VM)
    gy = k / G
    after = gam > 2.0 * gy
    tau = h[:, 3]
    dev = np.abs(tau[after] - k).max() / k
    drop = largest_drop(tau) / k
    good = dev <= 0.01 and drop <= 1e-3
    ok &= good
    print(f"  von Mises (Y {VM['sigma_y']:.0f} Pa, k = Y/sqrt(3) = {k:.2f} Pa): after "
          f"yield s_xy within {100 * dev:.2e} % of k (<= 1 %), largest drop "
          f"{100 * drop:.2e} % of k (no oscillation, <= 0.1 %): "
          f"{'PASS' if good else 'FAIL'}")
    print(f"    s_xy at gamma 1: {tau[-1]:.4f} Pa; normal stresses s_xx {h[-1, 0]:.3f}, "
          f"s_yy {h[-1, 1]:.3f} Pa; cells differ by at most {spread:.1e} Pa")

    p0 = -50.0e3
    eta, xi, _ = models.dp_coefficients(DP)
    lim = xi * DP["c"] - eta * p0
    gam, h, spread, _ = shear_case("DruckerPrager", DP, p0=p0)
    sqj = np.sqrt(batch._j2(h))
    yielded = gam > 2.0 * lim / G
    dev = np.abs(sqj[yielded] - lim).max() / lim
    pdev = np.abs(h[:, :3].mean(axis=1) - p0).max() / abs(p0)
    drop = largest_drop(h[:, 3]) / lim
    good = dev <= 0.01 and drop <= 1e-3 and pdev <= 1e-6
    ok &= good
    print(f"  Drucker-Prager (phi 30, c 5 kPa, psi 0) under p = -50 kPa: after yield "
          f"sqrt(J2) within {100 * dev:.2e} % of xi c - eta p = {lim:.1f} Pa, mean "
          f"stress change {pdev:.1e}, s_xy largest drop {100 * drop:.2e} %: "
          f"{'PASS' if good else 'FAIL'}")
    return ok


def check_plastic_step_1d(n=1000):
    """The WP 1.2 von Mises impact on an n x 2 mesh against lr1d_plastic."""
    import lr1d_plastic as lp
    K, G = models.bulk_shear(lp.VM)
    M = K + 4.0 * G / 3.0
    ce = math.sqrt(M / lp.RHO)
    hel = lp.VM["sigma_y"] * M / (2.0 * G)
    V = 3.0 * hel / (lp.RHO * ce)
    v0 = 0.25 * ce
    bar = lp.impact_bar(n, lp.VM, "vonMises", V, v0)
    box = PlasticBox2D(n, 2, lp.LENGTH, 2.0 * lp.LENGTH / n, lp.RHO, "vonMises", lp.VM)
    box.vx = np.tile(bar.v, (2, 1))
    m0, (px0, py0) = box.total_mass(), box.total_momentum()
    pscale = (box.cell_mass() * np.abs(box.vx)).sum()
    t_end = 0.22 * lp.LENGTH / ce          # WP 1.2's end time
    nsteps = int(np.ceil(t_end / bar.stable_dt(0.5)))
    dt = t_end / nsteps
    for _ in range(nsteps):
        bar.step(dt)
        box.step(dt)
    jump = np.abs(bar.sig6).max()
    ds = np.abs(box.s - bar.sig6[None]).max() / jump
    dv = np.abs(box.vx - bar.v).max() / V
    de = np.abs(box.epq - bar.eps_p_eq).max() / bar.eps_p_eq.max()
    dr = np.abs(box.rho - bar.rho).max() / lp.RHO
    side = max(np.abs(box.vy).max(), np.abs(box.s[..., 3]).max())
    dm = abs(box.total_mass() - m0) / m0
    px, py = box.total_momentum()
    dp = max(abs(px - px0), abs(py - py0)) / pscale
    worst = max(ds, dv, de, dr)
    good = worst <= 1e-9 and side == 0.0 and dm <= 1e-12 and dp <= 1e-12
    print(f"  {n} cells, {nsteps} steps to WP 1.2's end time; largest difference from "
          f"lr1d_plastic: stress {ds:.1e} (of the largest stress), velocity {dv:.1e} "
          f"(of V), eps_p_eq {de:.1e}, density {dr:.1e}; v_y and s_xy stay {side:.0e}")
    print(f"  plastic cells (eps_p_eq above 1 % of its peak) "
          f"{int((box.epq[0] > 0.01 * box.epq.max()).sum())}; mass change {dm:.1e}, "
          f"momentum change {dp:.1e}")
    print(f"  reproduces WP 1.2 to 1e-9, conserves to 1e-12: {'PASS' if good else 'FAIL'}")
    return good


def check_oblique(n=100):
    """Plastic impact at 45 degrees against the same impact along x."""
    ok = True
    for name, params in (("vonMises", VM), ("DruckerPrager", DP), ("linearElastic", ELASTIC)):
        r0 = impact(0, n, name, params)
        r = impact(45, n, name, params)
        e = extra_extrema(r, r0)
        good = e["main"] <= 0.01 and r["dm"] <= 1e-12 and r["dp"] <= 1e-12
        ok &= good
        print(f"  {name:13s}: beyond the along-x range after the start-up "
              f"{100 * e['main']:.3f} % (<= 1 %), at the end {100 * e['end']:.3f} %; "
              f"start-up {100 * e['start']:.3f} % (reported); mass {r['dm']:.0e}, "
              f"momentum {r['dp']:.0e}: {'PASS' if good else 'FAIL'}")
    return ok


def verify():
    ok = True
    print("[1] Vectorised WP 1.4 models (constitutive/batch.py) against models.py,")
    print("    200 random paths of 30 increments crossing yield")
    for name, d in batch.check_batch().items():
        good = d <= 1e-12
        ok &= good
        print(f"  {name:13s}: largest difference {d:.1e} (<= 1e-12: "
              f"{'PASS' if good else 'FAIL'})")

    print("\n[2] Pre-stressed elastic disc in rigid rotation, one full turn, 100 cells")
    print("    across; interior = more than 6 cells inside the edge, outside = more")
    print("    than 6 cells outside")
    res, cost = disc_case()
    e_p = max(x["e_princ"] for x in res)
    e_t = max(x["e_tensor"] for x in res)
    e_a = max(x["e_ang"] for x in res)
    e_o = max(x["e_out"] for x in res)
    good = e_p <= 0.02 and e_t <= 0.02 and e_a <= 1.0 and e_o <= 0.01
    ok &= good
    print(f"  principal stresses within {100 * e_p:.2e} % (<= 2 %), tensor within "
          f"{100 * e_t:.2e} % (<= 2 %), direction within {e_a:.3f} deg (<= 1), "
          f"outside below {100 * e_o:.2e} % (<= 1 %): {'PASS' if good else 'FAIL'} "
          f"({cost:.0f} s)")

    print("\n[3] Simple shear to shear strain 1 (prescribed velocity, 2000 steps)")
    ok &= check_simple_shear()

    print("\n[4] Plastic step wave along x on the 2D mesh vs the WP 1.2 bar")
    ok &= check_plastic_step_1d()

    print("\n[5] Impact at 45 degrees vs along x, 100 cells per metre, streaming")
    print("    0.25 c_e diagonally (decision 0013, item 6)")
    ok &= check_oblique()

    print("\nWP 1.3b checks:", "ALL PASS" if ok else "FAILED")
    return ok


# ----------------------------------------------------------------------------
# Formulation-decision measurements (not pass criteria; decision 0013)
# ----------------------------------------------------------------------------

def study():
    print("[A] Disc without the rotation correction (rotation='none'), 50 cells")
    res, _ = disc_case(cells=50, margin=5, band=4, rotation="none", quiet=True)
    for x in res:
        print(f"  {x['turn']:.2f} turn: principal direction error {x['e_ang']:.1f} deg, "
              f"tensor error {100 * x['e_tensor']:.1f} %")
    print("\n[B] Disc edge after one turn against resolution (Jaumann)")
    for cells, margin in ((50, 5), (100, 8)):
        res, _ = disc_case(cells=cells, margin=margin, band=6, quiet=True)
        print(f"  {cells} cells across: edge cells (10-90 %) {res[-1]['width']}, "
              f"{res[-1]['steps']} steps")
    print("\n[C] Impact at 45 degrees, extra extrema over the along-x run:")
    print("    symmetric split (prototype) and an unsplit proxy (both sweeps from")
    print("    the same state, increments added, as a face-based unsplit update)")
    for name, params in (("linearElastic", ELASTIC), ("vonMises", VM)):
        for split, cfl in (("symmetric", 0.5), ("additive", 0.5), ("additive", 0.4),
                           ("additive", 0.3)):
            r0 = impact(0, 100, name, params, split=split, cfl=cfl, travel=0.5)
            r = impact(45, 100, name, params, split=split, cfl=cfl, travel=0.5)
            e = extra_extrema(r, r0)
            print(f"  {name:13s} {split:9s} CFL {cfl}: start-up {100 * e['start']:.2f} %, "
                  f"after {100 * e['main']:.2f} %, end {100 * e['end']:.2f} %")
    print("\n[D] S system solved with the P impedance (more dissipation), symmetric")
    for name, params in (("linearElastic", ELASTIC), ("vonMises", VM)):
        r0 = impact(0, 100, name, params, s_impedance="p")
        r = impact(45, 100, name, params, s_impedance="p")
        e = extra_extrema(r, r0)
        print(f"  {name:13s}: start-up {100 * e['start']:.2f} %, after "
              f"{100 * e['main']:.2f} %")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--verify", action="store_true", help="run the WP 1.3b checks")
    p.add_argument("--study", action="store_true",
                   help="formulation-decision measurements (decision 0013)")
    args = p.parse_args(argv)
    if args.verify:
        return 0 if verify() else 1
    if args.study:
        return study()
    p.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
