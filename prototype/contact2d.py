"""WP 4.0: 2D experiments comparing two ways a rigid body pushes on Eulerian material.

A deliberately small, first-order, unsplit cell-centred scheme (plane strain, linear
elastic, material plus void, gravity), the same family as prototype/void2d.py but
written so the two couplings differ only where they touch the body:

  "penalty"   (LS-DYNA CLIS / S-ALE penalty, Abaqus CEL general contact): the body is
              invisible to the face fluxes; material that overlaps the body (cell fill
              a above the open fraction phi) gets a spring force k d along the body
              normal (d = overlap depth), a damper, and Coulomb friction.
  "points"    (LS-DYNA-style coupling points, our reconstruction, not verified against
              the manual): points every h / nquad on the body surface track the
              penetration of the material (integrated relative normal velocity),
              pressure k d or a curve, damper, friction; optional leak control
              (ileak) adds the "penalty" overlap spring.
  "riemann"   (cut cell / embedded boundary, Natarajan et al. 2021): faces are cut to
              their open length, the body surface inside a cut cell is one more face
              where the material meets a moving wall (one-sided acoustic Riemann
              solution with the body velocity; Coulomb cap on the shear, no tension).
              No mass crosses a cut. Cut cells with phi < 0.5 are mixed with their
              neighbour away from the body; a cell the body covers hands its content
              to that neighbour (mass and momentum kept).

Bodies: axis-aligned rectangles moving at constant velocity (exact cell and face
fractions) and fixed half-planes (exact by polygon clipping).

Run `python3 contact2d.py` for the three experiments and a summary table.
"""
import math
import sys

import numpy as np

FULL = 0.999      # fill above which a cell carries stress (void-first closure)
EPS = 1e-9


# --- geometry ------------------------------------------------------------------------

def clip_poly(poly, a, b, c):
    """Part of polygon poly (list of (x, y)) where a x + b y <= c."""
    out = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        fp = a * p[0] + b * p[1] - c
        fq = a * q[0] + b * q[1] - c
        if fp <= 0:
            out.append(p)
        if fp * fq < 0:
            t = fp / (fp - fq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def poly_area(poly):
    if len(poly) < 3:
        return 0.0
    s = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        s += x0 * y1 - x1 * y0
    return 0.5 * abs(s)


def seg_inside(p, q, a, b, c):
    """Length fraction of segment pq with a x + b y <= c."""
    fp = a * p[0] + b * p[1] - c
    fq = a * q[0] + b * q[1] - c
    if fp <= 0 and fq <= 0:
        return 1.0
    if fp > 0 and fq > 0:
        return 0.0
    t = fp / (fp - fq)
    return t if fp <= 0 else 1.0 - t


def overlap_1d(a0, a1, b0, b1):
    return np.clip(np.minimum(a1, b1) - np.maximum(a0, b0), 0.0, None)


class Rect:
    """Axis-aligned rigid rectangle [x0, x1] x [y0, y1] moving at velocity v."""

    def __init__(self, x0, x1, y0, y1, v=(0.0, 0.0)):
        self.box0 = np.array([x0, x1, y0, y1], float)
        self.v = np.array(v, float)

    def box(self, t):
        b = self.box0.copy()
        b[:2] += self.v[0] * t
        b[2:] += self.v[1] * t
        return b

    def fractions(self, g, t):
        """Body fraction per cell b, covered fraction of x faces bx (ny, nx+1) and
        y faces by (ny+1, nx), outward body normal per cell."""
        x0, x1, y0, y1 = self.box(t)
        h = g.h
        ox = overlap_1d(g.xf[:-1], g.xf[1:], x0, x1) / h
        oy = overlap_1d(g.yf[:-1], g.yf[1:], y0, y1) / h
        b = oy[:, None] * ox[None, :]
        inx = ((g.xf >= x0) & (g.xf <= x1)).astype(float)
        iny = ((g.yf >= y0) & (g.yf <= y1)).astype(float)
        # a face exactly on the body edge counts as open (material side)
        inx = np.where((np.abs(g.xf - x0) < 1e-12) | (np.abs(g.xf - x1) < 1e-12), 0.0, inx)
        iny = np.where((np.abs(g.yf - y0) < 1e-12) | (np.abs(g.yf - y1) < 1e-12), 0.0, iny)
        bx = oy[:, None] * inx[None, :]
        by = iny[:, None] * ox[None, :]
        # outward normal: from the body centre towards the cell, snapped to the
        # nearest side (a rectangle's faces)
        X, Y = np.meshgrid(g.xc, g.yc)
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        dx = (X - cx) / max(x1 - x0, EPS)
        dy = (Y - cy) / max(y1 - y0, EPS)
        nx = np.where(np.abs(dx) >= np.abs(dy), np.sign(dx), 0.0)
        ny = np.where(np.abs(dx) < np.abs(dy), np.sign(dy), 0.0)
        return b, bx, by, np.stack([nx, ny], -1)

    def points(self, g, t, spacing):
        """Coupling points (LS-DYNA style) at time t, fixed in the body: the
        same order every call, so a penetration history can be kept per point."""
        P, N, A = rect_points(self.box0, spacing)
        return P + self.v * t, N, A

    def depth(self, g, t, n):
        """Distance from each cell centre to the body surface along n (inside > 0)."""
        x0, x1, y0, y1 = self.box(t)
        X, Y = np.meshgrid(g.xc, g.yc)
        d = np.where(n[..., 0] > 0, x1 - X, 0.0) + np.where(n[..., 0] < 0, X - x0, 0.0) \
            + np.where(n[..., 1] > 0, y1 - Y, 0.0) + np.where(n[..., 1] < 0, Y - y0, 0.0)
        return d


def rect_points(box, spacing):
    """Coupling points on the four sides of a rectangle box = (x0, x1, y0, y1):
    positions, outward normals, and the side length each stands for."""
    x0, x1, y0, y1 = box
    P, N, A = [], [], []
    for (a, b), nrm in ((((x0, y0), (x1, y0)), (0.0, -1.0)), (((x1, y0), (x1, y1)), (1.0, 0.0)),
                        (((x1, y1), (x0, y1)), (0.0, 1.0)), (((x0, y1), (x0, y0)), (-1.0, 0.0))):
        L = math.hypot(b[0] - a[0], b[1] - a[1])
        n = max(1, int(round(L / spacing)))
        for k in range(n):
            f = (k + 0.5) / n
            P.append((a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])))
            N.append(nrm)
            A.append(L / n)
    return np.array(P), np.array(N), np.array(A)


class HalfPlane:
    """Fixed rigid half-plane n . x <= c (n the outward unit normal)."""

    def __init__(self, n, c):
        self.n = np.array(n, float) / np.linalg.norm(n)
        self.c = c
        self.v = np.zeros(2)
        self._cache = None

    def fractions(self, g, t):
        if self._cache is None:
            a, bb, c = self.n[0], self.n[1], self.c
            b = np.zeros((g.ny, g.nx))
            bx = np.zeros((g.ny, g.nx + 1))
            by = np.zeros((g.ny + 1, g.nx))
            for j in range(g.ny):
                for i in range(g.nx):
                    x0, x1, y0, y1 = g.xf[i], g.xf[i + 1], g.yf[j], g.yf[j + 1]
                    sq = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
                    b[j, i] = poly_area(clip_poly(sq, a, bb, c)) / g.h ** 2
            for j in range(g.ny):
                for i in range(g.nx + 1):
                    bx[j, i] = seg_inside((g.xf[i], g.yf[j]), (g.xf[i], g.yf[j + 1]), a, bb, c)
            for j in range(g.ny + 1):
                for i in range(g.nx):
                    by[j, i] = seg_inside((g.xf[i], g.yf[j]), (g.xf[i + 1], g.yf[j]), a, bb, c)
            nrm = np.broadcast_to(self.n, (g.ny, g.nx, 2)).copy()
            self._cache = (b, bx, by, nrm)
        return self._cache

    def depth(self, g, t, n):
        X, Y = np.meshgrid(g.xc, g.yc)
        return self.c - (self.n[0] * X + self.n[1] * Y)

    def points(self, g, t, spacing):
        """Coupling points along the part of the wall line inside the box."""
        tv = np.array([-self.n[1], self.n[0]])
        p0 = self.c * self.n
        W, H = g.xf[-1], g.yf[-1]
        # parameter range of the line inside [0, W] x [0, H]
        lo, hi = -1e9, 1e9
        for k, (a, b) in enumerate(((0.0, W), (0.0, H))):
            if abs(tv[k]) > 1e-14:
                s1, s2 = (a - p0[k]) / tv[k], (b - p0[k]) / tv[k]
                lo, hi = max(lo, min(s1, s2)), min(hi, max(s1, s2))
        L = hi - lo
        n = max(1, int(round(L / spacing)))
        sk = lo + (np.arange(n) + 0.5) * L / n
        P = p0 + sk[:, None] * tv
        return P, np.broadcast_to(self.n, P.shape).copy(), np.full(n, L / n)


# --- the solver ----------------------------------------------------------------------

class Grid:
    def __init__(self, nx, ny, h):
        self.nx, self.ny, self.h = nx, ny, h
        self.xf = np.arange(nx + 1) * h
        self.yf = np.arange(ny + 1) * h
        self.xc = 0.5 * (self.xf[:-1] + self.xf[1:])
        self.yc = 0.5 * (self.yf[:-1] + self.yf[1:])


def riemann(uL, sL, zL, uR, sR, zR):
    """Acoustic Riemann solution (u*, s*) between two homogenised states; a side
    with zero impedance is a free surface (s = 0)."""
    zs = zL + zR
    ok = zs > 0
    zsafe = np.where(ok, zs, 1.0)
    u = np.where(ok, (zL * uL + zR * uR + sR - sL) / zsafe, 0.5 * (uL + uR))
    s = np.where(ok, (zR * sL + zL * sR + zL * zR * (uR - uL)) / zsafe, 0.0)
    return u, s


class Box:
    """Material plus void on a Cartesian grid with one rigid body.

    Walls of the box are frictionless (mirror) except as set in self.walls_mu.
    State: mass m, momentum via v (vx, vy), mass-weighted stress s (xx, yy, xy)."""

    def __init__(self, grid, E, nu, rho0, body, scheme, mu_body=0.0, gravity=(0.0, 0.0),
                 pfac=1.0, cfl=0.4, nquad=2, damp=0.0, curve=None, ileak=False):
        self.g = grid
        self.E, self.nu, self.rho0 = E, nu, rho0
        self.lam = E * nu / ((1 + nu) * (1 - 2 * nu))
        self.G = E / (2 * (1 + nu))
        self.M = self.lam + 2 * self.G
        self.zp = math.sqrt(self.M * rho0)
        self.zs = math.sqrt(self.G * rho0)
        self.c = math.sqrt(self.M / rho0)
        self.body, self.scheme = body, scheme
        self.mu_b = mu_body
        self.grav = np.array(gravity, float)
        self.kpen = pfac * self.M / grid.h          # penalty stiffness [Pa/m]
        # coupling points (scheme "points", LS-DYNA CLIS/S-ALE style): nquad points
        # per cell length of body surface, damping as a fraction of the impedance
        # (critical for a wave leaving the surface), optional pressure-penetration
        # curve (d, p) pairs replacing k d
        self.nquad, self.damp, self.curve = nquad, damp, curve
        self.dpen = None
        self.ileak = ileak       # points plus the cell-overlap penalty (leak control)
        self.dt = cfl * grid.h / self.c
        ny, nx = grid.ny, grid.nx
        self.m = np.zeros((ny, nx))
        self.v = np.zeros((ny, nx, 2))
        self.s = np.zeros((ny, nx, 3))
        self.t = 0.0
        self.force = np.zeros(2)                     # force of the material on the body

    def fill(self, inside):
        X, Y = np.meshgrid(self.g.xc, self.g.yc)
        b = self.body.fractions(self.g, 0.0)[0]
        self.m = self.rho0 * self.g.h ** 2 * np.minimum(inside(X, Y), 1.0 - b)

    # geometry seen by the face fluxes
    def geom(self, t):
        g = self.g
        b, bx, by, nb = self.body.fractions(g, t)
        if self.scheme in ("penalty", "points"):
            z = np.zeros_like
            return z(b), z(bx), z(by), b, nb
        return b, bx, by, b, nb

    def step(self):
        g, dt, h = self.g, self.dt, self.g.h
        ny, nx = g.ny, g.nx
        tm = self.t + 0.5 * dt
        b, bx, by, btrue, nb = self.geom(tm)
        phi = 1.0 - b
        V = h * h
        vol = np.maximum(phi * V, 1e-30)
        a = self.m / (self.rho0 * V)                 # fill of the cell
        fill = np.where(phi > EPS, a / np.maximum(phi, EPS), 0.0)
        frac = np.clip(fill, 0.0, 1.0)
        has = self.m > 0.0
        msafe = np.where(has, self.m, 1.0)
        sig = np.where(has[..., None], self.s / msafe[..., None], 0.0)
        sig = np.where((fill >= FULL)[..., None], sig, 0.0)    # void-first
        vel = self.v
        zp = frac * self.zp
        zs = frac * self.zs

        F = np.zeros((ny, nx, 2))                    # force on each cell
        L = np.zeros((ny, nx, 2, 2))                 # sum A u* (x) n
        flux_x = np.zeros((ny, nx + 1))              # u*_n A on faces
        flux_y = np.zeros((ny + 1, nx))

        # x faces (interior), n = +x: normal component x, tangential y, shear sxy
        ap = (1.0 - bx[:, 1:-1]) * h
        un, sn = riemann(vel[:, :-1, 0], sig[:, :-1, 0], zp[:, :-1],
                         vel[:, 1:, 0], sig[:, 1:, 0], zp[:, 1:])
        ut, st = riemann(vel[:, :-1, 1], sig[:, :-1, 2], zs[:, :-1],
                         vel[:, 1:, 1], sig[:, 1:, 2], zs[:, 1:])
        fx, fy = ap * sn, ap * st
        F[:, :-1, 0] += fx
        F[:, :-1, 1] += fy
        F[:, 1:, 0] -= fx
        F[:, 1:, 1] -= fy
        L[:, :-1, 0, 0] += ap * un
        L[:, :-1, 1, 0] += ap * ut
        L[:, 1:, 0, 0] -= ap * un
        L[:, 1:, 1, 0] -= ap * ut
        flux_x[:, 1:-1] = ap * un
        # y faces, n = +y: normal y, tangential x, shear sxy
        ap = (1.0 - by[1:-1, :]) * h
        un, sn = riemann(vel[:-1, :, 1], sig[:-1, :, 1], zp[:-1],
                         vel[1:, :, 1], sig[1:, :, 1], zp[1:])
        ut, st = riemann(vel[:-1, :, 0], sig[:-1, :, 2], zs[:-1],
                         vel[1:, :, 0], sig[1:, :, 2], zs[1:])
        F[:-1, :, 1] += ap * sn
        F[:-1, :, 0] += ap * st
        F[1:, :, 1] -= ap * sn
        F[1:, :, 0] -= ap * st
        L[:-1, :, 1, 1] += ap * un
        L[:-1, :, 0, 1] += ap * ut
        L[1:, :, 1, 1] -= ap * un
        L[1:, :, 0, 1] -= ap * ut
        flux_y[1:-1, :] = ap * un
        # box walls: frictionless, fixed (mirror state: u* = 0 normal, shear 0)
        for side, (j, sgn) in {"left": (0, -1), "right": (-1, 1)}.items():
            ap = (1.0 - bx[:, j]) * h
            s_n = sig[:, j, 0] - sgn * zp[:, j] * vel[:, j, 0]
            F[:, j, 0] += sgn * ap * s_n
            L[:, j, 1, 0] += sgn * ap * vel[:, j, 1]
        for side, (j, sgn) in {"bottom": (0, -1), "top": (-1, 1)}.items():
            ap = (1.0 - by[j, :]) * h
            s_n = sig[j, :, 1] - sgn * zp[j] * vel[j, :, 1]
            F[j, :, 1] += sgn * ap * s_n
            L[j, :, 0, 1] += sgn * ap * vel[j, :, 0]

        self.force[:] = 0.0
        vb = self.body.v
        if self.scheme == "riemann":
            # body face per cut cell: A_b n_b = -(sum of open face areas times n)
            Ax = (1.0 - bx[:, 1:]) * h - (1.0 - bx[:, :-1]) * h
            Ay = (1.0 - by[1:, :]) * h - (1.0 - by[:-1, :]) * h
            An = -np.stack([Ax, Ay], -1)              # A_b n_b, n_b into the body
            A = np.linalg.norm(An, axis=-1)
            cut = (A > 1e-12 * h) & (phi > EPS)
            n = np.where(cut[..., None], An / np.where(cut, A, 1.0)[..., None], 0.0)
            tt = np.stack([-n[..., 1], n[..., 0]], -1)
            # stress traction on the body face, normal and shear components
            tr = np.stack([sig[..., 0] * n[..., 0] + sig[..., 2] * n[..., 1],
                           sig[..., 2] * n[..., 0] + sig[..., 1] * n[..., 1]], -1)
            s_nn = (tr * n).sum(-1)
            s_nt = (tr * tt).sum(-1)
            u_n = (vel * n).sum(-1)
            u_t = (vel * tt).sum(-1)
            w_n = (vb * n).sum(-1)
            w_t = (vb * tt).sum(-1)
            # wall on the + side of the face: s* = sL + zL (uw - uL)
            sn_star = s_nn + zp * (w_n - u_n)
            sep = sn_star > 0.0                        # no tension: separate
            sn_star = np.where(sep, 0.0, sn_star)
            un_star = np.where(sep, u_n + np.where(zp > 0, -s_nn / np.where(zp > 0, zp, 1), 0),
                               w_n)
            st_star = s_nt + zs * (w_t - u_t)
            cap = self.mu_b * np.abs(sn_star)
            st_star = np.clip(st_star, -cap, cap)
            ut_star = u_t + np.where(zs > 0, (st_star - s_nt) / np.where(zs > 0, zs, 1), 0)
            fb = (sn_star[..., None] * n + st_star[..., None] * tt) * A[..., None]
            fb = np.where(cut[..., None], fb, 0.0)
            F += fb
            self.force -= fb.sum((0, 1))
            ub = un_star[..., None] * n + ut_star[..., None] * tt
            L += np.where(cut[..., None, None],
                          ub[..., :, None] * An[..., None, :], 0.0)
        if self.scheme == "points":
            self._couple_points(F, tm, vel, sig)
        if self.scheme == "penalty" or (self.scheme == "points" and self.ileak):
            # penalty on the overlap of material and body (cells with b > 0)
            # material volume (mass over the density its stress implies), so that
            # elastic compression of a full cell is not taken for overlap
            ev = (sig[..., 0] + sig[..., 1]) / (2.0 * (self.lam + self.G))
            av = a * (1.0 + ev)
            ov = np.where(btrue > 0.0, np.clip(av - (1.0 - btrue), 0.0, None), 0.0)
            d = ov * h                                            # overlap depth
            on = ov > 0.0
            # stiffness and damping scale with the cell's material mass (as LS-DYNA
            # scales PFAC on the critical stiffness of the coupled mass)
            am = np.clip(a, 0.0, 1.0)
            fn = self.kpen * am * d * h                           # force per unit depth
            rel = vel - vb
            rn = (rel * nb).sum(-1)
            damp = np.where(on & (rn < 0), -self.zp * h * am * rn, 0.0)
            fn = np.where(on, fn + damp, 0.0)
            tvec = rel - rn[..., None] * nb
            tm_ = np.linalg.norm(tvec, axis=-1)
            ft = np.minimum(self.mu_b * fn, self.zs * h * tm_)
            tdir = np.where(tm_[..., None] > 0, tvec / np.where(tm_ > 0, tm_, 1)[..., None], 0)
            fb = fn[..., None] * nb - ft[..., None] * tdir
            F += fb
            self.force -= fb.sum((0, 1))

        # gravity, momentum
        F += self.m[..., None] * self.grav
        mom = self.m[..., None] * vel + dt * F
        # stress (elastic, plane strain, full cells only)
        Lg = L / vol[..., None, None]
        D = 0.5 * (Lg + np.swapaxes(Lg, -1, -2))
        tr_ = D[..., 0, 0] + D[..., 1, 1]
        ds = np.stack([self.lam * tr_ + 2 * self.G * D[..., 0, 0],
                       self.lam * tr_ + 2 * self.G * D[..., 1, 1],
                       2 * self.G * D[..., 0, 1]], -1)
        full = (fill >= FULL)[..., None]
        sig_new = np.where(full, sig + dt * ds, 0.0)
        s_m = self.m[..., None] * sig_new                     # mass-weighted

        # remap: donor cell with the face velocities (fluxes are u*_n A)
        m = self.m.copy()
        dens = np.where(phi > EPS, self.m / vol, 0.0)         # homogenised density
        vnew = np.where(has[..., None], mom / msafe[..., None], 0.0)
        q = {"m": m, "p": self.m[..., None] * vnew, "s": s_m}
        newq = {k: v.copy() for k, v in q.items()}
        per_m = {"p": vnew, "s": sig_new}
        mfs = []
        out = np.zeros_like(m)
        for axis, fl in ((1, flux_x[:, 1:-1]), (0, flux_y[1:-1, :])):
            V_sw = fl * dt
            if axis == 1:
                dL, dR = dens[:, :-1], dens[:, 1:]
                sl = (slice(None), slice(None, -1)), (slice(None), slice(1, None))
            else:
                dL, dR = dens[:-1], dens[1:]
                sl = (slice(None, -1), slice(None)), (slice(1, None), slice(None))
            mf = np.where(V_sw > 0, V_sw * dL, V_sw * dR)
            out[sl[0]] += np.clip(mf, 0, None)
            out[sl[1]] += np.clip(-mf, 0, None)
            mfs.append((sl, mf))
        # a donor never gives more than it has (small cut cells)
        r = np.where(out > m, m / np.where(out > 0, out, 1.0), 1.0)
        for sl, mf in mfs:
            mf = np.where(mf > 0, mf * r[sl[0]], mf * r[sl[1]])
            newq["m"][sl[0]] -= mf
            newq["m"][sl[1]] += mf
            for k in ("p", "s"):
                qL, qR = per_m[k][sl[0]], per_m[k][sl[1]]
                qf = np.where((mf > 0)[..., None], qL, qR) * mf[..., None]
                newq[k][sl[0]] -= qf
                newq[k][sl[1]] += qf
        m = np.maximum(newq["m"], 0.0)
        p, s_m = newq["p"], newq["s"]

        # riemann scheme: small and covered cells (body pose at the end of step)
        self.t += dt
        if self.scheme == "riemann":
            bE, _, _, nE = self.body.fractions(g, self.t)
            phiE = 1.0 - bE
            small = (phiE < 0.5) & (m > 0)
            for j, i in zip(*np.nonzero(small)):
                n = nE[j, i]
                if abs(n[0]) >= abs(n[1]):
                    jj, ii = j, i + int(np.sign(n[0]))
                else:
                    jj, ii = j + int(np.sign(n[1])), i
                if not (0 <= jj < ny and 0 <= ii < nx):
                    continue
                if phiE[j, i] <= EPS:                 # covered: hand everything over
                    m[jj, ii] += m[j, i]
                    p[jj, ii] += p[j, i]
                    s_m[jj, ii] += s_m[j, i]
                    m[j, i] = 0.0
                    p[j, i] = 0.0
                    s_m[j, i] = 0.0
                else:                                 # mix velocity and stress
                    mt = m[j, i] + m[jj, ii]
                    if mt > 0:
                        vm = (p[j, i] + p[jj, ii]) / mt
                        sm = (s_m[j, i] + s_m[jj, ii]) / mt
                        p[j, i], p[jj, ii] = m[j, i] * vm, m[jj, ii] * vm
                        s_m[j, i], s_m[jj, ii] = m[j, i] * sm, m[jj, ii] * sm
        has = m > 1e-14 * self.rho0 * V
        ms = np.where(has, m, 1.0)
        self.m = np.where(has, m, 0.0)
        self.v = np.where(has[..., None], p / ms[..., None], 0.0)
        self.s = np.where(has[..., None], s_m, 0.0)

    def _couple_points(self, F, t, vel, sig):
        """LS-DYNA-style coupling: points on the body surface track how far the
        material has moved into the body (integrated relative normal velocity of
        the material interpolated at the point); a penalty pressure p(d) plus a
        damper acts on the material there and the reaction on the body."""
        g, h, dt = self.g, self.g.h, self.dt
        P, N, A = self.body.points(g, t, h / self.nquad)
        if self.dpen is None:
            self.dpen = np.zeros(len(P))
        # bilinear weights over the four nearest cell centres holding material
        fx = P[:, 0] / h - 0.5
        fy = P[:, 1] / h - 0.5
        i0 = np.floor(fx).astype(int)
        j0 = np.floor(fy).astype(int)
        tx, ty = fx - i0, fy - j0
        idx, wts = [], []
        for di, dj, w in ((0, 0, (1 - tx) * (1 - ty)), (1, 0, tx * (1 - ty)),
                          (0, 1, (1 - tx) * ty), (1, 1, tx * ty)):
            ii = np.clip(i0 + di, 0, g.nx - 1)
            jj = np.clip(j0 + dj, 0, g.ny - 1)
            w = w * (self.m[jj, ii] > 0)
            idx.append((jj, ii))
            wts.append(w)
        wsum = sum(wts)
        inbox = (P[:, 0] > 0) & (P[:, 0] < g.xf[-1]) & (P[:, 1] > 0) & (P[:, 1] < g.yf[-1])
        on = (wsum > 1e-12) & inbox
        ws = np.where(on, wsum, 1.0)
        um = sum(w[:, None] * vel[jj, ii] for (jj, ii), w in zip(idx, wts)) / ws[:, None]
        rel = um - self.body.v
        rin = -(rel * N).sum(-1)                         # into the body > 0
        self.dpen = np.where(on, np.maximum(self.dpen + dt * rin, 0.0), 0.0)
        d = self.dpen
        if self.curve is None:
            p = self.kpen * d
        else:
            cd, cp = np.array(self.curve).T
            slope = (cp[-1] - cp[-2]) / (cd[-1] - cd[-2])
            p = np.where(d <= cd[-1], np.interp(d, cd, cp), cp[-1] + slope * (d - cd[-1]))
        p = np.where(d > 0, np.maximum(p + self.damp * self.zp * rin, 0.0), 0.0)
        fn = p * A
        tvec = rel + rin[:, None] * N
        tmag = np.linalg.norm(tvec, axis=-1)
        ft = np.minimum(self.mu_b * fn, self.zs * A * tmag)
        tdir = np.where(tmag[:, None] > 0, tvec / np.where(tmag > 0, tmag, 1)[:, None], 0)
        fp = fn[:, None] * N - ft[:, None] * tdir     # on the material
        # the force goes to the cells in proportion to weight times mass, so
        # nearly empty cells are not flung (LS-DYNA spreads it to the nodes)
        wm = [w * self.m[jj, ii] for (jj, ii), w in zip(idx, wts)]
        wms = sum(wm)
        wms = np.where(wms > 0, wms, 1.0)
        for (jj, ii), w in zip(idx, wm):
            np.add.at(F, (jj, ii), (w / wms)[:, None] * fp)
        self.force -= fp.sum(0)

    # diagnostics
    def leak(self):
        """Material inside the body: mass beyond the open volume (overfill above 2 %
        elastic compression is counted), over the material mass."""
        b = self.body.fractions(self.g, self.t)[0]
        V = self.g.h ** 2
        cap = 1.02 * self.rho0 * V * (1.0 - b)
        return float(np.clip(self.m - cap, 0, None)[b > 0].sum() / self.m.sum())

    def mass(self):
        return float(self.m.sum())


# --- experiments -----------------------------------------------------------------------

MAT = dict(E=10e6, nu=0.3, rho0=1800.0)


def piston(scheme, h=0.005, offset=0.0, pfac=1.0, v=2.0, t_end=None, **kw):
    """Full-width rigid piston pushed down into an elastic layer at speed v. Exact
    contact stress until the wave returns from the floor: z v (1D impedance)."""
    W, H = 0.1, 0.3
    nx, ny = int(round(W / h)), int(round(H / h))
    g = Grid(nx, ny, h)
    top = 0.2
    y0 = top + offset * h
    body = Rect(-1.0, W + 1.0, y0, H + 1.0, v=(0.0, -v))
    box = Box(g, MAT["E"], MAT["nu"], MAT["rho0"], body, scheme, pfac=pfac, **kw)
    box.fill(lambda X, Y: np.clip((y0 - (Y - 0.5 * h)) / h, 0.0, 1.0))
    t_end = t_end or 1.6 * top / box.c       # before the floor reflection returns
    m0 = box.mass()
    F, T, Lk = [], [], []
    while box.t < t_end:
        box.step()
        F.append(box.force[1] / W)               # contact stress (compression > 0)
        T.append(box.t)
        Lk.append(box.leak())
    F, T = np.array(F), np.array(T)
    exact = box.zp * v
    win = T > 0.3 * top / box.c
    mean = F[win].mean() / exact
    noise = F[win].std() / exact
    return dict(mean=mean, noise=noise, leak=max(Lk), mass=abs(box.mass() / m0 - 1))


def block(scheme, h=0.005, offset=0.0, pfac=1.0, v=2.0, cells=6, **kw):
    """Rigid block (4 x 4 cells) translating along x through material at rest for
    `cells` cells; returns leakage, drag noise and mass error."""
    W, H = 0.3, 0.15
    nx, ny = int(round(W / h)), int(round(H / h))
    g = Grid(nx, ny, h)
    x0 = 0.08 + offset * h
    y0 = 0.055 + offset * h
    body = Rect(x0, x0 + 4 * h, y0, y0 + 4 * h, v=(v, 0.0))
    box = Box(g, MAT["E"], MAT["nu"], MAT["rho0"], body, scheme, pfac=pfac, **kw)
    box.fill(lambda X, Y: np.ones_like(X))
    m0 = box.mass()
    t_end = cells * h / v
    F, Lk = [], []
    while box.t < t_end:
        box.step()
        F.append(-box.force[0])                  # drag (against the motion > 0)
        Lk.append(box.leak())
    F = np.array(F)
    n2 = len(F) // 2
    late = F[n2:]
    ref = box.zp * v * 4 * h                    # impedance force on the front face
    # noise: step-to-step jumps of the drag (high-frequency part), over z v A
    noise = np.abs(np.diff(late)).mean() / ref
    return dict(drag=late.mean() / ref, noise=noise, leak=max(Lk),
                leak_end=Lk[-1], mass=abs(box.mass() / m0 - 1))


def incline(scheme, h=0.005, pfac=1.0, theta=30.0, delta=20.0, t_end=0.12, **kw):
    """Elastic block on a fixed frictional wall inclined at theta (not aligned with the
    mesh), under gravity. Exact sliding acceleration g (sin theta - tan delta cos
    theta); measured from the block's velocity along the slope over the second half."""
    W, H = 0.4, 0.3
    nx, ny = int(round(W / h)), int(round(H / h))
    g = Grid(nx, ny, h)
    th = math.radians(theta)
    n = (-math.sin(th), math.cos(th))           # outward normal of the wall (up-left)
    # wall: n . x <= c, surface through (0, 0.05)
    c = n[1] * 0.05
    body = HalfPlane(n, c)
    grav = (0.0, -9.81)
    box = Box(g, MAT["E"], MAT["nu"], MAT["rho0"], body, scheme,
              mu_body=math.tan(math.radians(delta)), gravity=grav, pfac=pfac, **kw)
    # block: 0.12 along the slope, 0.05 thick, centred at s = 0.18 along the wall
    tv = np.array([math.cos(th), math.sin(th)])
    nv = np.array(n)
    origin = np.array([0.0, 0.05])

    def inside(X, Y):
        P = np.stack([X - origin[0], Y - origin[1]], -1)
        s = P @ tv
        d = P @ nv
        return ((s > 0.12) & (s < 0.24) & (d > 0.0) & (d < 0.05)).astype(float)
    box.fill(inside)
    m0 = box.mass()
    T, U, Lk = [], [], []
    while box.t < t_end:
        box.step()
        mt = box.m.sum()
        U.append(-((box.m[..., None] * box.v).sum((0, 1)) / mt @ tv))   # down the slope
        T.append(box.t)
        Lk.append(box.leak())
    T, U = np.array(T), np.array(U)
    half = T > 0.5 * t_end
    acc = np.polyfit(T[half], U[half], 1)[0]
    exact = 9.81 * (math.sin(th) - math.tan(math.radians(delta)) * math.cos(th))
    return dict(acc=acc / exact, leak=max(Lk), mass=abs(box.mass() / m0 - 1))


def main(argv=None):
    rows = []
    for scheme, pf in (("riemann", None), ("penalty", 0.1), ("penalty", 1.0),
                       ("penalty", 10.0)):
        label = scheme if pf is None else f"penalty k={pf} M/h"
        pf = pf or 1.0
        p = [piston(scheme, offset=o, pfac=pf) for o in (0.0, 0.25, 0.5)]
        bl = [block(scheme, offset=o, pfac=pf) for o in (0.0, 0.5)]
        inc = incline(scheme, pfac=pf)
        rows.append((label, p, bl, inc))
        print(f"{label}:")
        for o, r in zip((0.0, 0.25, 0.5), p):
            print(f"  piston offset {o}: stress/zv {r['mean']:.4f} noise {r['noise']:.4f}"
                  f" leak {r['leak']:.2e} mass {r['mass']:.1e}")
        for o, r in zip((0.0, 0.5), bl):
            print(f"  block offset {o}: drag/zvA {r['drag']:.3f} noise {r['noise']:.4f}"
                  f" leak max {r['leak']:.2e} end {r['leak_end']:.2e} mass {r['mass']:.1e}")
        print(f"  incline 30 deg, delta 20: acc/exact {inc['acc']:.3f} leak {inc['leak']:.2e}"
              f" mass {inc['mass']:.1e}")
        sys.stdout.flush()
    return rows


if __name__ == "__main__":
    main(sys.argv[1:])
