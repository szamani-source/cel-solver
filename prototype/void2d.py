#!/usr/bin/env python3
"""WP 3.2p: one material plus void in the 2D prototype (free surface).

Builds on ``silo2d.py`` (WP 3.1p: well-balanced gravity, wall boundaries,
settling) and adds a material volume fraction ``al`` per cell. The rules
(decision 0021):

1. **State.** ``rho`` is the bulk density (cell mass over cell volume), so the
   mass and every mass-weighted remap are exactly as before; ``al`` is the
   material volume fraction; the material density is rho / al. Velocity,
   stress and history are per unit material mass. A cell with al <= EPS is
   empty: no mass, zero velocity, zero stress.
2. **Where the material sits (1D SLIC per sweep).** In a sweep a partly
   filled cell holds its material as one slab against the side of its fuller
   neighbour (``lean`` = +1 right, -1 left); when both neighbours are equally
   full the material is spread evenly (a layer parallel to the sweep). Walls
   mirror the cell next to them; free boundaries count as empty.
3. **Faces.** A face is a contact face when both sides have material that
   touches it (non-empty, not leaning away). In the Lagrange step a partly
   filled cell is homogenised: traction al s, impedance al z_material, mass
   al rho_m h, so its wave speed is the material's and the explicit update
   is stable for any al (with the material's own stress and impedance a cell
   of al = 1e-3 blows up: its stress acts on a thousandth of the mass).
   Contact faces get the acoustic Riemann solution of these. Every other face is a free surface: zero traction (s_hat = -G in
   the gravity sweep), and the velocity from the one-sided relation of the
   side that touches it (material impedance); if neither touches it, from the
   side(s) holding material (so a slab leaning away is not sheared by the
   face), zero between empty cells.
   Reconstruction is first order in and next to partly filled cells.
4. **Void closure, (a) of the brief: volume changes go to the void first.**
   Compression of a partly filled cell closes its void before it strains the
   material. Expansion first unloads the material elastically to zero normal
   stress, and the rest opens void, so cohesionless material separates
   instead of pulling, also inside the domain. Shear goes to the material.
5. **Remap.** Volume flux through a face = the material part of the swept
   region of the donor cell, from the slab position of rule 2 (exact for a
   slab moving into void). Mass flux = volume flux times the donor's
   Lagrangian material density (MUSCL as before where the donor and its
   neighbours are full). Velocity, stress and history are carried with the
   mass flux as before (MC limiter, superbee for history), first order next
   to partly filled cells. A cell overfilled by round-off or by inconsistent
   slab positions is clipped to al = 1 (its mass is kept).

Rules 6 (round-off dust) and 7 (gravity in first-order cells) and the
reasons for all of them are in docs/decisions/0021-void-treatment.md.

8. **Shear cap (WP 3.2q, decision 0023).** The jump term of the shear Riemann
   solution may not push the face stress outside the Mohr-Coulomb surface;
   a capped face slips (each side keeps its own characteristic velocity).
   Without it the jump term acts as a mesh-independent viscosity in one-cell
   shear bands and slow granular flow depends on E.
9. **Outlet (WP 3.3p, decision 0024).** ``outlet`` marks rows whose lo end is
   an outflow face: zero traction, empty ghosts with free-end parities, and
   the material swept across the face leaves the domain (counted in
   ``VoidBox.discharged``).
10. **Low-Mach jump terms (WP 3.3p, decision 0024; off by default).** Both
   Riemann jump terms z (vR - vL) / 2 are multiplied by min(1, |v| / c), so
   the numerical viscosity is rho |v| dv rather than rho c dv.

``python3 void2d.py --verify`` runs the WP 3.2p checks
(tests/test_proto_void.sh); ``--collapse`` the column collapse checks (WP 3.2q,
slope open in WP 3.2r) (tests/manual/test_proto_collapse.sh). SI units, stress positive in tension, thickness
1 m, arrays (ny, nx).
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

import silo2d as sd  # noqa: E402
from lr1d import limiter  # noqa: E402
from lr1d_plastic import superbee  # noqa: E402
from lr2d_plastic import NG, PERM_Y, gradients, swept_average  # noqa: E402
from lr2d_plastic import rotate_inplane  # noqa: E402
from silo2d import mirror, pad_ends, parity, weight_profile, boundary_face  # noqa: E402

EPS = 1e-12          # al at or below: empty
FULL = 1.0 - 1e-12   # al at or above: full
DISS = 1.0           # study knob: scale of the Riemann jump terms
SHEAR_CAP = True     # rule 8 (decision 0023): face shear traction capped at yield
DUST = 1e-9          # al below which a cell with mass is round-off dust


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

def pad_alpha(al, lo, hi):
    """Volume fraction with ghosts: mirrored at walls (a wall neither
    attracts nor repels the slab), empty beyond free ends."""
    p = mirror(al)
    if lo == "free":
        p[:, :NG] = 0.0
    if hi == "free":
        p[:, -NG:] = 0.0
    return p


def lean_of(alp, n):
    """Slab side of each cell (+1 right, -1 left, 0 spread) from the padded
    volume fraction: towards the fuller neighbour."""
    left = alp[:, NG - 1:NG - 1 + n]
    right = alp[:, NG + 1:NG + 1 + n]
    dif = right - left
    return np.where(np.abs(dif) > 1e-9, np.sign(dif), 0.0)


def smooth_mask(alp):
    """1 where a cell and both neighbours are full (second order allowed),
    for the padded cells 1 .. -2 (the cells that get slopes)."""
    full = alp >= FULL
    return (full[:, 1:-1] & full[:, :-2] & full[:, 2:]).astype(float)


def out_overlap(A, H, d, lean, side):
    """Material volume in the swept region of length |d| at one end of a cell
    of Lagrangian length H holding material length A as a slab (lean)."""
    ad = np.abs(d)
    spread = np.where(H > 0.0, A / np.where(H > 0.0, H, 1.0), 0.0) * ad
    near = np.minimum(ad, A)                       # slab against this end
    far = np.maximum(0.0, ad - (H - A))            # slab against the other end
    toward = lean > 0 if side == "right" else lean < 0
    away = lean < 0 if side == "right" else lean > 0
    return np.where(toward, near, np.where(away, far, spread))


def shear_cap(fn, ft, seff, zs, vt_star, snt_star, Gf, params):
    """Rule 8 (decision 0023): the jump term of the shear Riemann solution
    (z_s (vR - vL) / 2, a numerical viscosity) may not push the face stress
    outside the Mohr-Coulomb surface. The face stress is the Riemann normal
    and shear tractions with the mean of the two cells' tt and zz stresses;
    the shear jump is scaled back (bisection) until the yield function is
    no larger than without it (or zero). A capped face slips: each side
    takes the velocity of its own characteristic at the capped traction, as
    at a frictional wall. Returns the new (vt_star, snt_star) and the
    (left side, right side) tangential face velocities, all per face."""
    sphi = math.sin(math.radians(params["phi"]))
    ccos = params.get("c", 0.0) * math.cos(math.radians(params["phi"]))
    stt = mirror(seff[..., 1])
    szz = mirror(seff[..., 2])
    ftt = 0.5 * (stt[:, 1:-2] + stt[:, 2:-1])
    fzz = 0.5 * (szz[:, 1:-2] + szz[:, 2:-1])
    snn_face = fn[4] + fn[5] + Gf                  # actual normal face stress
    tb, tj = ft[4], ft[5]

    def yf(snt):
        cen = 0.5 * (snn_face + ftt)
        rad = np.sqrt((0.5 * (snn_face - ftt)) ** 2 + snt ** 2)
        s1 = np.maximum(cen + rad, fzz)
        s3 = np.minimum(cen - rad, fzz)
        return 0.5 * (s1 - s3) + 0.5 * (s1 + s3) * sphi - ccos

    lim = np.maximum(yf(tb), 0.0)
    bad = yf(tb + tj) > lim
    if not bad.any():
        return vt_star, snt_star, None
    lo = np.zeros_like(tb)
    hi = np.ones_like(tb)
    for _ in range(20):
        mid = 0.5 * (lo + hi)
        ok = yf(tb + mid * tj) <= lim
        lo = np.where(ok, mid, lo)
        hi = np.where(ok, hi, mid)
    t_new = np.where(bad, tb + lo * tj, snt_star)
    zl, zr = zs[:, 1:-2], zs[:, 2:-1]
    va = ft[1] + (t_new - ft[0]) / np.where(zl > 0.0, zl, 1.0)
    vb = ft[3] - (t_new - ft[2]) / np.where(zr > 0.0, zr, 1.0)
    return vt_star, t_new, (np.where(bad, va, vt_star), np.where(bad, vb, vt_star))


# ----------------------------------------------------------------------------
# One sweep with void
# ----------------------------------------------------------------------------

def void_sweep(f, h0, dt, mat, ends, model, params, grav=0.0, mu=0.0,
               rotation="jaumann", no_tension=True, outlet=None, low_mach=False):
    """silo2d.silo_sweep with a material volume fraction f["al"] (rules 1 to 5
    in the module docstring). Returns (new state, eps_p_eq made, boundary
    record).

    outlet: optional bool array per row (WP 3.3p, rule 9): those rows' lo end
    is an outflow face instead of the lo wall: zero traction, empty ghosts,
    and material swept across it leaves the domain. The mass that left (per
    unit face area) is returned in the record as rec["out"].
    low_mach: rule 10 (WP 3.3p), the Riemann jump terms scaled by the face
    Mach number."""
    M, G, rho0 = mat
    c_mat = math.sqrt(M / rho0)
    rows, n = f["rho"].shape
    lo, hi = ends
    rho, al = f["rho"], f["al"]
    m = rho * h0                                     # mass per unit face area
    has = al > EPS
    full = al >= FULL
    msafe = np.where(m > 0.0, m, 1.0)
    # homogenised cell: traction al * s, impedance al * z_material, mass
    # al * rho_m * h, so the wave speed is the material's for any al and the
    # explicit update is stable in nearly empty cells. rho_m is bounded so
    # that round-off in al cannot make the impedance meaningless.
    rho_m = np.where(has, rho / np.where(has, al, 1.0), rho0)
    rho_m = np.clip(rho_m, 0.5 * rho0, 2.0 * rho0)
    zeff = np.where(has, al * np.sqrt(M * rho_m), 0.0)
    zseff = np.where(has, al * np.sqrt(G * rho_m), 0.0)
    seff = f["s"] * al[..., None]
    Gf, Gc = weight_profile(rho, h0, grav)

    if outlet is not None and not outlet.any():
        outlet = None
    alp = pad_alpha(al, lo, hi)
    if outlet is not None:
        alp[outlet, :NG] = 0.0                       # beyond an outlet: empty
    lean = lean_of(alp, n)
    lean = np.where(full, 0.0, lean)
    mask = smooth_mask(alp)                          # (rows, n + 2)

    # Material touching each face: right end of cell j-1, left end of cell j
    touch_r = has & ~((~full) & (lean < 0))         # touches its right face
    touch_l = has & ~((~full) & (lean > 0))         # touches its left face

    h = np.full((rows, n + 2 * NG), h0)
    psn_lo, pvn_lo = parity(lo, "n")
    psn_hi, pvn_hi = parity(hi, "n")
    pst_lo, pvt_lo = parity(lo, "t")
    pst_hi, pvt_hi = parity(hi, "t")
    snn = pad_ends(seff[..., 0] - Gc, psn_lo, psn_hi, -Gf[:, 0], -Gf[:, -1])
    vn = pad_ends(f["vn"], pvn_lo, pvn_hi)
    snt = pad_ends(seff[..., 3], pst_lo, pst_hi)
    vt = pad_ends(f["vt"], pvt_lo, pvt_hi)
    if outlet is not None:
        # outlet rows: ghosts of a free end (rule 9)
        o_rows = outlet[:, None]
        fsn, fvn = parity("free", "n")
        fst, fvt = parity("free", "t")
        for arr, q, par, c in ((snn, seff[..., 0] - Gc, fsn, -Gf[:, 0]),
                               (vn, f["vn"], fvn, 0.0), (snt, seff[..., 3], fst, 0.0),
                               (vt, f["vt"], fvt, 0.0)):
            g = pad_ends(q, par, 1, c)
            arr[:, :NG] = np.where(o_rows, g[:, :NG], arr[:, :NG])
    zp = mirror(zeff)
    zs = mirror(zseff)
    mp = mirror(m)

    def characteristic(s, v, z):
        zc, hc = z[:, 1:-1], h[:, 1:-1]
        zsafe = np.where(zc > 0.0, zc, 1.0)
        gs_l, gs_r = gradients(s, h)
        gv_l, gv_r = gradients(v, h)
        wr = limiter(gs_l - zc * gv_l, gs_r - zc * gv_r)
        wl = limiter(gs_l + zc * gv_l, gs_r + zc * gv_r)
        return 0.5 * (wr + wl) * hc * mask, 0.5 * (wl - wr) / zsafe * hc * mask

    def riemann(s, v, z, ds, dv, mod):
        mc, hc, zc = mp[:, 1:-1], h[:, 1:-1], z[:, 1:-1]
        mcs = np.where(mc > 0.0, mc, 1.0)
        s_p = s[:, 1:-1] + 0.5 * dt * mod * dv / hc
        v_p = v[:, 1:-1] + 0.5 * dt * ds / mcs
        sL = (s_p + 0.5 * ds)[:, :-1]
        vL = (v_p + 0.5 * dv)[:, :-1]
        zL = zc[:, :-1]
        sR = (s_p - 0.5 * ds)[:, 1:]
        vR = (v_p - 0.5 * dv)[:, 1:]
        zR = zc[:, 1:]
        zsum = zL + zR
        zsum = np.where(zsum > 0.0, zsum, 1.0)
        v_star = (zL * vL + zR * vR + DISS * (sR - sL)) / zsum
        s_base = (zR * sL + zL * sR) / zsum
        s_jump = DISS * zL * zR * (vR - vL) / zsum
        if low_mach:
            # rule 10 (WP 3.3p, decision 0024): jump term scaled by the face
            # Mach number of the material speed, so its stress is of order
            # rho |v| dv instead of rho c dv (Thornber et al. 2008)
            ma = np.maximum(np.abs(vL), np.abs(vR)) / c_mat
            s_jump = s_jump * np.minimum(ma, 1.0)
        return v_star, s_base + s_jump, (sL, vL, sR, vR, s_base, s_jump)

    ds_n, dv_n = characteristic(snn, vn, zp)
    if grav > 0.0:
        # First-order cells in the gravity sweep: the s_hat difference across
        # the cell is -w phi (w the cell's own weight). phi = 0 keeps s_hat
        # uniform (exact for a supported cell in equilibrium, whose centre
        # stress carries at least half its own weight); phi = 1 keeps the
        # actual stress uniform (exact in free fall, stress 0). Between,
        # phi = 1 - |s_nn| / (w / 2), clipped to [0, 1] (rule 7).
        wp = mirror(grav * rho * h0)[:, 1:-1]
        sp = mirror(seff[..., 0])[:, 1:-1]
        phi = np.where(wp > 0.0, np.clip(1.0 - np.abs(sp) / np.where(wp > 0.0, 0.5 * wp, 1.0),
                                         0.0, 1.0), 0.0)
        ds_n = ds_n + (1.0 - mask) * (-wp * phi)
    vn_star, snn_star, fn = riemann(snn, vn, zp, ds_n, dv_n, M)
    ds_t, dv_t = characteristic(snt, vt, zs)
    vt_star, snt_star, ft = riemann(snt, vt, zs, ds_t, dv_t, G)
    slip = None
    if SHEAR_CAP and "phi" in params:
        vt_star, snt_star, slip = shear_cap(fn, ft, seff, zs, vt_star, snt_star, Gf, params)

    # --- free-surface faces (interior faces 1 .. n-1) -------------------------
    pL = touch_r[:, :-1]                     # left cell touches the face
    pR = touch_l[:, 1:]
    free = ~(pL & pR)
    if no_tension:
        # material without tensile strength separates where the contact
        # solution would pull (actual normal traction positive)
        free |= snn_star[:, 1:-1] + Gf[:, 1:-1] > 0.0
    if free.any():
        zmL, zmR = zeff[:, :-1], zeff[:, 1:]
        zsL, zsR = zseff[:, :-1], zseff[:, 1:]
        sfree = -Gf[:, 1:-1]

        def one_sided(sstar, F, zl, zr):
            sL, vL, sR, vR = (a[:, 1:-1] for a in F[:4])
            vl = vL + (sstar - sL) / np.where(zl > 0.0, zl, 1.0)
            vr = vR - (sstar - sR) / np.where(zr > 0.0, zr, 1.0)
            hL, hR = has[:, :-1], has[:, 1:]
            other = np.where(hL & hR, 0.5 * (vl + vr), np.where(hL, vl, np.where(hR, vr, 0.0)))
            return np.where(pL, vl, np.where(pR, vr, other))

        vfn = one_sided(sfree, fn, zmL, zmR)
        vft = one_sided(0.0 * sfree, ft, zsL, zsR)
        inner = (slice(None), slice(1, -1))
        snn_star[inner] = np.where(free, sfree, snn_star[inner])
        vn_star[inner] = np.where(free, vfn, vn_star[inner])
        snt_star[inner] = np.where(free, 0.0, snt_star[inner])
        vt_star[inner] = np.where(free, vft, vt_star[inner])

    # capped faces slip: the cell on each side sees its own tangential velocity
    vtA = vt_star.copy()          # seen by the cell on the left of the face
    vtB = vt_star.copy()          # seen by the cell on the right
    if slip is not None:
        inner = (slice(None), slice(1, -1))
        vtA[inner] = np.where(free, vtA[inner], slip[0][:, 1:-1])
        vtB[inner] = np.where(free, vtB[inner], slip[1][:, 1:-1])

    # --- boundary faces --------------------------------------------------------
    rec = {}
    for side, kind, j, sgn, k, cell, pres in (
            ("lo", lo, 0, -1, 0, 0, touch_l[:, 0]),
            ("hi", hi, -1, 1, -1, -1, touch_r[:, -1])):
        o = 2 if sgn < 0 else 0
        zI = np.where(pres, zeff[:, cell], 1.0)
        zIs = np.where(pres, zseff[:, cell], 1.0)
        s_b, v_b = boundary_face(kind, sgn, fn[o][:, k], fn[o + 1][:, k], zI, Gf[:, j])
        o_here = outlet is not None and side == "lo"
        if o_here:
            s_o, v_o = boundary_face("free", sgn, fn[o][:, k], fn[o + 1][:, k], zI, Gf[:, j])
            s_b, v_b = np.where(outlet, s_o, s_b), np.where(outlet, v_o, v_b)
        s_b = np.where(pres, s_b, -Gf[:, j])
        v_b = np.where(pres, v_b, 0.0)
        snn_star[:, j], vn_star[:, j] = s_b, v_b
        normal = s_b + Gf[:, j]
        t_b, w_b = boundary_face(kind, sgn, ft[o][:, k], ft[o + 1][:, k], zIs,
                                 None, normal=normal, mu=mu)
        if o_here:
            t_o, w_o = boundary_face("free", sgn, ft[o][:, k], ft[o + 1][:, k], zIs,
                                     None, normal=normal, mu=mu)
            t_b, w_b = np.where(outlet, t_o, t_b), np.where(outlet, w_o, w_b)
        t_b = np.where(pres, t_b, 0.0)
        w_b = np.where(pres, w_b, 0.0)
        snt_star[:, j], vt_star[:, j] = t_b, w_b
        vtA[:, j], vtB[:, j] = w_b, w_b
        rec[side] = (normal, t_b.copy(), v_b.copy(), w_b.copy())

    # --- Lagrange update ----------------------------------------------------
    vn_new = np.where(has, f["vn"] + dt * (snn_star[:, 1:] - snn_star[:, :-1]) / msafe, 0.0)
    vt_new = np.where(has, f["vt"] + dt * (snt_star[:, 1:] - snt_star[:, :-1]) / msafe, 0.0)

    d = dt * vn_star
    h_new = h0 + d[:, 1:] - d[:, :-1]
    dh = h_new - h0
    A0 = np.where(full, h0, np.where(has, al * h0, 0.0))
    vd = h0 - A0
    s_nn_old = f["s"][..., 0]
    cap = A0 * np.maximum(0.0, -s_nn_old) / M
    if not no_tension:
        cap = np.where(full, np.inf, cap)    # a full cell of a solid stretches
    dmat = np.where(dh < 0.0, dh - np.maximum(dh, -vd), np.minimum(dh, cap))
    dmat = np.where(has, dmat, 0.0)
    A_new = A0 + dmat
    amid = A0 + 0.5 * dmat
    hmid = 0.5 * (h0 + h_new)
    gam = np.where(has, dt * (vtA[:, 1:] - vtB[:, :-1]) / hmid, 0.0)
    deps = np.zeros((rows, n, 6))
    deps[..., 0] = np.where(has, dmat / np.where(amid > 0.0, amid, 1.0), 0.0)
    deps[..., 3] = gam
    s_old, ep_old = f["s"], f["ep"]
    if rotation == "jaumann":
        a = 0.5 * gam
        s_old = rotate_inplane(s_old, a)
        ep_old = rotate_inplane(ep_old, a, 0.5)
    act = has & ((deps[..., 0] != 0.0) | (gam != 0.0))
    s_new, ep_new, epq_new = s_old.copy(), ep_old.copy(), f["epq"].copy()
    if act.any():
        s_new[act], ep_new[act], epq_new[act], _ = model(
            deps[act], s_old[act], ep_old[act], f["epq"][act], params)
    made = (m * (epq_new - f["epq"])).sum()

    # --- remap: slab volume fluxes, no flux through the boundary faces --------
    dR = d[:, 1:]                      # right face of each cell
    dL = d[:, :-1]
    outR = np.where(dR > 0.0, out_overlap(A_new, h_new, dR, lean, "right"), 0.0)
    outL = np.where(dL < 0.0, out_overlap(A_new, h_new, dL, lean, "left"), 0.0)
    # a cell that would keep less than 1e-9 of a cell of material gives all
    # of it (no round-off dust left behind a moving edge)
    tot_out = outR + outL
    drain = (A_new > 0.0) & (tot_out > 0.0) & (A_new - tot_out < 1e-9 * h0)
    scale = np.where(drain, A_new / np.where(tot_out > 0.0, tot_out, 1.0), 1.0)
    outR, outL = outR * scale, outL * scale
    vf = np.zeros((rows, n + 1))
    dI = d[:, 1:-1]
    vf[:, 1:-1] = np.where(dI > 0.0, outR[:, :-1], -outL[:, 1:])
    if outlet is not None:
        # rule 9: material swept across an outlet face leaves the domain
        vf[:, 0] = np.where(outlet & (d[:, 0] < 0.0), -outL[:, 0], 0.0)
    # round-off face motion (for example 1e-19 m from a 1e-17 m/s velocity)
    # must not seed dust in empty cells: no flux below 1e-9 of a cell
    vf = np.where(np.abs(vf) > 1e-9 * h0, vf, 0.0)
    rhoL = np.where(A_new > 0.0, m / np.where(A_new > 0.0, A_new, 1.0), 0.0)
    # MUSCL density in fully filled regions (identical to silo_sweep there)
    hp = mirror(h_new)
    rp = mirror(rhoL)
    rho_muscl = swept_average(rp, hp, d)
    donor_smooth = np.zeros((rows, n + 1), dtype=bool)
    ms = mask[:, 1:-1] > 0.0                           # (rows, n)
    donor_smooth[:, 1:-1] = np.where(dI > 0.0, ms[:, :-1], ms[:, 1:])
    rho_donor = np.zeros((rows, n + 1))
    rho_donor[:, 1:-1] = np.where(dI > 0.0, rhoL[:, :-1], rhoL[:, 1:])
    rho_donor[:, 0] = rhoL[:, 0]
    mf = vf * np.where(donor_smooth, rho_muscl, rho_donor)
    if outlet is None:
        mf[:, 0] = 0.0
    else:
        mf[:, 0] = np.where(outlet, mf[:, 0], 0.0)
    mf[:, -1] = 0.0
    rec["out"] = -mf[:, 0]

    dmf = mf[:, :-1] - mf[:, 1:]
    m_new = m + dmf
    m_new = np.where(m_new > 0.0, m_new, 0.0)
    out = {"rho": m_new / h0}
    a_tot = A_new + vf[:, :-1] - vf[:, 1:]
    al_new = np.clip(a_tot / h0, 0.0, 1.0)
    al_new = np.where(al_new >= FULL, 1.0, al_new)
    # dust: round-off mass left in a cell that has (almost) no material
    # volume. Its volume is set from its mass at the reference density, and it
    # carries no stress and no velocity of its own (rule 6, decision 0021).
    dust = (m_new > 0.0) & (al_new < DUST)
    al_new = np.where(dust, m_new / (rho0 * h0), al_new)
    al_new = np.where(m_new > 0.0, al_new, 0.0)
    out["al"] = al_new
    mnew_safe = np.where(m_new > 0.0, m_new, 1.0)

    def carry(val, par=(1, 1), lim=limiter):
        if val.ndim == 3:
            pv = mirror(val)
            mk = mask[..., None]
            mpp, mff = mp[..., None], mf[..., None]
        else:
            pv = pad_ends(val, *par)
            mk, mpp, mff = mask, mp, mf
        mps = np.where(mpp > 0.0, mpp, 1.0)
        gl, gr = gradients(pv, mps[..., 0] if val.ndim == 3 else mps)
        sl = lim(gl, gr) * mk
        qc, mc = pv[:, 1:-1], mpp[:, 1:-1]
        am = np.abs(mff)
        from_left = qc[:, :-1] + sl[:, :-1] * 0.5 * (mc[:, :-1] - am)
        from_right = qc[:, 1:] - sl[:, 1:] * 0.5 * (mc[:, 1:] - am)
        flux = mff * np.where(mff > 0.0, from_left, from_right)
        mm = m[..., None] if val.ndim == 3 else m
        mn = mnew_safe[..., None] if val.ndim == 3 else mnew_safe
        keep = (m_new > 0.0)[..., None] if val.ndim == 3 else m_new > 0.0
        return np.where(keep, (mm * val + flux[:, :-1] - flux[:, 1:]) / mn, 0.0)

    out["vn"] = np.where(dust, 0.0, carry(vn_new, (pvn_lo, pvn_hi)))
    out["vt"] = np.where(dust, 0.0, carry(vt_new, (pvt_lo, pvt_hi)))
    out["s"] = np.where(dust[..., None], 0.0, carry(s_new))
    out["ep"] = carry(ep_new, lim=superbee)
    out["epq"] = carry(epq_new, lim=superbee)
    return out, made, rec


# ----------------------------------------------------------------------------
# Box with void
# ----------------------------------------------------------------------------

class VoidBox(sd.SiloBox):
    """SiloBox with a material volume fraction ``al`` (1 full, 0 empty)."""

    KEYS = ("rho", "al", "vx", "vy", "s", "ep", "epq")

    def __init__(self, *args, **kw):
        super().__init__(*args, **kw)
        self.al = np.ones((self.ny, self.nx))
        # no tensile strength: cohesionless, or a tension cut-off at zero
        p = self.params
        self.no_tension = p.get("c", None) == 0.0 or p.get("sigma_t", None) == 0.0
        self.low_mach = False     # rule 10 (decision 0024)
        self.outlet = None        # bool per column: bottom face is an outflow (rule 9)
        self.discharged = 0.0     # mass that left through the outlet (kg per m depth)
        self._out = []

    def fill(self, inside, K0=None):
        """Material where inside(X, Y, h) gives the filled fraction per cell
        (0 to 1), void elsewhere; stress K0-lithostatic if K0 is given."""
        X, Y = np.meshgrid(self.xc, self.yc)
        al = np.clip(inside(X, Y), 0.0, 1.0)
        self.al = np.where(al >= FULL, 1.0, np.where(al > EPS, al, 0.0))
        self.rho = self.rho0 * self.al
        self.vx = np.zeros_like(self.al)
        self.vy = np.zeros_like(self.al)
        self.s = np.zeros((self.ny, self.nx, 6))
        self.ep = np.zeros((self.ny, self.nx, 6))
        self.epq = np.zeros((self.ny, self.nx))
        if K0 is not None:
            Gf, Gc = weight_profile(self.rho.T, self.dy, self.g)
            alT = self.al.T
            syy = np.where(alT > EPS, (Gc - Gf[:, -1:]) / np.where(alT > EPS, alT, 1.0),
                           0.0).T                  # material stress = cell stress / al
            self.s[..., 1] = syy
            self.s[..., 0] = K0 * syy
            self.s[..., 2] = K0 * syy

    def sweep(self, st, dt, axis):
        mat = (self.M, self.G, self.rho0)
        if axis == 0:
            fl = {"rho": st["rho"], "al": st["al"], "vn": st["vx"], "vt": st["vy"],
                  "s": st["s"], "ep": st["ep"], "epq": st["epq"]}
            g, made, rec = void_sweep(fl, self.dx, dt, mat,
                                      (self.walls["left"], self.walls["right"]),
                                      self.model, self.params, 0.0, self.mu, self.rotation,
                                      self.no_tension, low_mach=self.low_mach)
            self._rec.append(("left", rec["lo"]))
            self._rec.append(("right", rec["hi"]))
            return {"rho": g["rho"], "al": g["al"], "vx": g["vn"], "vy": g["vt"],
                    "s": g["s"], "ep": g["ep"], "epq": g["epq"]}, made * self.dy
        tr = lambda a: np.swapaxes(a, 0, 1)                      # noqa: E731
        fl = {"rho": tr(st["rho"]), "al": tr(st["al"]), "vn": tr(st["vy"]),
              "vt": tr(st["vx"]), "s": tr(st["s"][..., PERM_Y]),
              "ep": tr(st["ep"][..., PERM_Y]), "epq": tr(st["epq"])}
        g, made, rec = void_sweep(fl, self.dy, dt, mat,
                                  (self.walls["bottom"], self.walls["top"]),
                                  self.model, self.params, self._g_now, self.mu,
                                  self.rotation, self.no_tension, outlet=self.outlet,
                                  low_mach=self.low_mach)
        self._out.append(rec["out"].sum() * self.dx)
        self._rec.append(("bottom", rec["lo"]))
        self._rec.append(("top", rec["hi"]))
        return {"rho": tr(g["rho"]), "al": tr(g["al"]), "vx": tr(g["vt"]),
                "vy": tr(g["vn"]), "s": tr(g["s"])[..., PERM_Y],
                "ep": tr(g["ep"])[..., PERM_Y], "epq": tr(g["epq"])}, made * self.dx

    def step(self, dt):
        """Symmetric split (x-y and y-x averaged, mass-weighted as
        PlasticBox2D.step); al is averaged by volume. Damping and the
        boundary record as SiloBox.step."""
        self._g_now = self.g_at(self.t + 0.5 * dt)
        self._rec = []
        self._out = []
        st = self.state()
        a, ma = self.sweep(st, dt, 0)
        a, mb = self.sweep(a, dt, 1)
        b, mc = self.sweep(st, dt, 1)
        b, md = self.sweep(b, dt, 0)
        ra, rb = a["rho"], b["rho"]
        rho = 0.5 * (ra + rb)
        tot = ra + rb
        ts = np.where(tot > 0.0, tot, 1.0)
        new = {"rho": rho, "al": 0.5 * (a["al"] + b["al"])}
        for k in ("vx", "vy", "epq", "s", "ep"):
            wa, wb, t2 = ra, rb, ts
            if a[k].ndim == 3:
                wa, wb, t2 = ra[..., None], rb[..., None], ts[..., None]
            new[k] = (wa * a[k] + wb * b[k]) / t2
        self.set_state(new)
        self.plastic_made += 0.5 * (ma + mb + mc + md)
        self.discharged += 0.5 * sum(self._out)    # the two orderings are averaged
        self.t += dt
        self.nstep += 1
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

    # --- diagnostics -----------------------------------------------------------

    def centroid(self):
        m = self.cell_mass()
        X, Y = np.meshgrid(self.xc, self.yc)
        M = m.sum()
        return (m * X).sum() / M, (m * Y).sum() / M

    def column_height(self):
        """Material height per column (m): sum of al dy."""
        return self.al.sum(axis=0) * self.dy

    def void_stress(self):
        """Largest |stress| in cells that are empty."""
        e = self.al <= EPS
        return np.abs(self.s[e]).max() if e.any() else 0.0


def make_box(nx, ny, lx, ly, model="MohrCoulomb", params=None, walls=None, delta_w=0.0,
             g=sd.G_ACC, rho=sd.RHO):
    walls = walls or {"left": "smooth", "right": "smooth", "bottom": "smooth",
                      "top": "smooth"}
    return VoidBox(nx, ny, lx, ly, rho, model, params or sd.T0_MC, walls=walls,
                   delta_w=delta_w, g=g)


def rect(x0, x1, y0, y1, dx, dy):
    """Cell fill fraction of the rectangle [x0, x1] x [y0, y1]."""
    def inside(X, Y):
        fx = np.clip((np.minimum(X + 0.5 * dx, x1) - np.maximum(X - 0.5 * dx, x0)) / dx, 0, 1)
        fy = np.clip((np.minimum(Y + 0.5 * dy, y1) - np.maximum(Y - 0.5 * dy, y0)) / dy, 0, 1)
        return fx * fy
    return inside


# ----------------------------------------------------------------------------
# Cases
# ----------------------------------------------------------------------------

def lake_at_rest(nx=10, ny=60, W=0.1, H=0.6, fill=0.455, nsteps=1000):
    """K0 column filled to `fill` (a partly filled top cell) between smooth
    walls in a closed box, full gravity, no damping."""
    box = make_box(nx, ny, W, H)
    box.fill(rect(-1.0, 2.0 * W, -1.0, fill, box.dx, box.dy), K0=sd.k0_elastic(0.3))
    m0 = box.total_mass()
    dt = box.stable_dt(0.5)
    vmax = 0.0
    for _ in range(nsteps):
        box.step(dt)
        vmax = max(vmax, np.abs(box.vx).max(), np.abs(box.vy).max())
    top = box.al > EPS
    return dict(vmax=vmax, void_stress=box.void_stress(),
                dm=abs(box.total_mass() - m0) / m0,
                partial=float(box.al[(box.al > EPS) & (box.al < FULL)].mean()),
                top_row=int(np.where(top.any(axis=1))[0].max()), box=box)


def translating_block(n_cells=50, v=5.0, cell=0.01, model="linearElastic"):
    """10 x 10 cell block in void, no gravity, velocity (v, 0), travels
    n_cells cells."""
    nx, ny = 80, 20
    box = make_box(nx, ny, nx * cell, ny * cell, model=model, g=0.0)
    x0, y0 = 5 * cell, 5 * cell
    box.fill(rect(x0, x0 + 10 * cell, y0, y0 + 10 * cell, cell, cell))
    box.vx = np.where(box.al > EPS, v, 0.0)
    m0 = box.total_mass()
    cx0, cy0 = box.centroid()
    t_end = n_cells * cell / v
    nsteps = int(math.ceil(t_end / box.stable_dt(0.5)))
    dt = t_end / nsteps
    smax = 0.0
    for _ in range(nsteps):
        box.step(dt)
        smax = max(smax, np.abs(box.s).max())
    cx, cy = box.centroid()
    # exact block after the travel: [x0 + v t, + 10 cells] x [y0, y0 + 10 cells]
    xa, xb = x0 + v * t_end, x0 + v * t_end + 10 * cell
    X, Y = np.meshgrid(box.xc, box.yc)
    dist_x = np.maximum(np.maximum(xa - X, X - xb), 0.0)
    dist_y = np.maximum(np.maximum(y0 - Y, Y - (y0 + 10 * cell)), 0.0)
    dist = np.maximum(dist_x, dist_y) / cell
    has = box.al > EPS
    spread = dist[has].max() + 0.5              # cells beyond the exact edges
    # edge sharpness: partly filled cells per row through the block
    partial = ((box.al > 1e-6) & (box.al < 1.0 - 1e-6)).sum(axis=1).max()
    return dict(dm=abs(box.total_mass() - m0) / m0, smax=smax,
                dcx=(cx - cx0) / (v * t_end) - 1.0, dcy=abs(cy - cy0),
                spread=spread, partial=int(partial), steps=nsteps,
                vrange=(box.vx[has].min(), box.vx[has].max()), box=box)


def falling_block(cell=0.01, model="linearElastic", t_after=0.03):
    """10 x 10 cell elastic block released at rest 0.1 m above a rigid
    smooth floor, full gravity. Free fall: centroid acceleration from a
    quadratic fit; landing: no mass loss, no material below the floor."""
    nx, ny = 20, 30
    box = make_box(nx, ny, nx * cell, ny * cell, model=model)
    y0 = 0.10
    box.fill(rect(5 * cell, 15 * cell, y0, y0 + 10 * cell, cell, cell))
    m0 = box.total_mass()
    dt = box.stable_dt(0.5)
    t_land = math.sqrt(2.0 * y0 / sd.G_ACC)
    ts, ys = [], []
    smax_fall = 0.0
    while box.t < 0.9 * t_land:
        box.step(dt)
        ts.append(box.t)
        ys.append(box.centroid()[1])
        smax_fall = max(smax_fall, np.abs(box.s).max())
    c = np.polyfit(np.array(ts), np.array(ys), 2)
    acc = -2.0 * c[0]
    ymin = 1e9
    vy_min = 0.0
    fb = []
    while box.t < t_land + t_after:
        box.step(dt)
        ymin = min(ymin, box.centroid()[1])
        m = box.cell_mass()
        vy_c = (m * box.vy).sum() / m.sum()
        vy_min = min(vy_min, vy_c)
        fb.append(-(box.face["bottom"][0]).sum() * box.dx)
    m = box.cell_mass()
    return dict(acc=acc, e_acc=abs(acc / sd.G_ACC - 1.0), smax_fall=smax_fall,
                dm=abs(box.total_mass() - m0) / m0, ymin=ymin,
                vy_end=(m * box.vy).sum() / m.sum(),
                v_impact=math.sqrt(2.0 * sd.G_ACC * y0), fb_max=max(fb), box=box,
                weight=m0 * sd.G_ACC)


def separation(cell=0.01, v=0.5, t_end=0.02):
    """Cohesionless block (T0 Mohr-Coulomb, zero stress, no gravity) pulled
    apart: left half moving at -v, right half at +v. It must separate (void
    opens inside the domain) without carrying tension."""
    nx, ny = 60, 10
    box = make_box(nx, ny, nx * cell, ny * cell, g=0.0)
    box.fill(rect(20 * cell, 40 * cell, -1.0, 1.0, cell, cell))
    X, _ = np.meshgrid(box.xc, box.yc)
    box.vx = np.where(box.al > EPS, np.where(X < 0.3, -v, v), 0.0)
    m0 = box.total_mass()
    p0 = box.total_momentum()[0]
    dt = box.stable_dt(0.5)
    smax = 0.0
    while box.t < t_end:
        box.step(dt)
        smax = max(smax, box.s[..., :3].max())
    mid = box.al[:, 28:32].min()
    gap_cells = int(((box.al < 1e-6).all(axis=0) & (box.xc > 0.2) & (box.xc < 0.4)).sum())
    inner = (box.xc > 0.25) & (box.xc < 0.35)
    void_cells = float((1.0 - box.al[:, inner]).sum(axis=1).mean())
    return dict(dm=abs(box.total_mass() - m0) / m0, smax=smax, gap=gap_cells,
                dp=abs(box.total_momentum()[0] - p0) / (m0 * v), mid=mid,
                expect_gap=2.0 * v * box.t / cell, void_cells=void_cells, box=box)


# ----------------------------------------------------------------------------
# Granular column collapse
# ----------------------------------------------------------------------------

COLLAPSE_E = 1.0e6    # softer than T0 (10 MPa): the run-out is rate-independent
COLLAPSE_DELTA_B = 35.0


def column_collapse(a, cells_short=20, cell=None, delta_b=COLLAPSE_DELTA_B,
                    E=COLLAPSE_E, phi=35.0, t_max=2.0, ke_stop=1e-3, log=None):
    """Column L_i x H_i = L_i x a L_i against the smooth left wall, released
    from K0 rest on a frictional floor (delta_b), runs until the kinetic
    energy is below ke_stop times its peak."""
    if a <= 1.0:
        H_i = 0.1
        L_i = H_i / a
    else:
        L_i = 0.1
        H_i = a * L_i
    cell = cell or min(H_i, L_i) / cells_short
    run = 1.2 * a * 1.6                       # room for 1.5 times the fit
    lx = L_i * (1.0 + run) + 4 * cell
    ly = H_i + 3 * cell
    nx = int(math.ceil(lx / cell))
    ny = int(math.ceil(ly / cell))
    params = dict(sd.T0_MC, E=E, phi=phi)
    walls = {"left": "smooth", "right": "smooth", "bottom": "friction", "top": "smooth"}
    box = make_box(nx, ny, nx * cell, ny * cell, params=params, walls=walls,
                   delta_w=delta_b)
    box.fill(rect(-1.0, L_i, -1.0, H_i, cell, cell), K0=sd.k0_elastic(params["nu"]))
    m0 = box.total_mass()
    t0 = time.perf_counter()
    ke_peak = 0.0
    hist = []
    dt = box.stable_dt(0.5)
    while box.t < t_max:
        dt = box.stable_dt(0.5)
        box.step(dt)
        if box.nstep % 20 == 0:
            ke = box.kinetic_energy()
            ke_peak = max(ke_peak, ke)
            hist.append((box.t, ke))
            if log:
                log(box, ke)
            if box.t > 0.05 and ke <= ke_stop * ke_peak:
                break
    sec = time.perf_counter() - t0
    hc = box.column_height()
    mcol = box.cell_mass().sum(axis=0) / (sd.RHO * cell * cell)   # cells' worth of mass
    front = np.where(mcol >= 0.5)[0]
    L_f = (front.max() + 1) * cell if front.size else 0.0
    H_f = hc.max()
    dep = hc >= cell
    slopes = np.degrees(np.arctan(np.abs(np.diff(hc)) / cell))
    sl_mask = dep[:-1] & dep[1:]
    # slope of the surface where the material is at least a cell deep, and
    # the same over two-cell baselines
    slopes2 = np.degrees(np.arctan(np.abs(hc[2:] - hc[:-2]) / (2 * cell)))
    sl2_mask = dep[:-2] & dep[2:]
    return dict(a=a, L_i=L_i, H_i=H_i, cell=cell, nx=nx, ny=ny, L_f=L_f, H_f=H_f,
                runout=(L_f - L_i) / L_i, dm=abs(box.total_mass() - m0) / m0,
                slope_max=slopes[sl_mask].max() if sl_mask.any() else 0.0,
                slope2_max=slopes2[sl2_mask].max() if sl2_mask.any() else 0.0,
                steps=box.nstep, t=box.t, seconds=sec, hc=hc,
                ke_end=hist[-1][1] / max(ke_peak, 1e-300), box=box,
                void_stress=box.void_stress())


# ----------------------------------------------------------------------------
# Acceptance checks for WP 3.2p
# ----------------------------------------------------------------------------

def verify():
    ok = True
    t_all = time.perf_counter()
    P = lambda good: "PASS" if good else "FAIL"            # noqa: E731

    print("[1] Lake at rest with a free surface inside the domain: K0 column between")
    print("    smooth walls in a closed 0.1 m x 0.6 m box, filled to 0.455 m (top cell")
    print("    half full), 10 mm cells, full gravity, no damping, 1000 steps")
    r = lake_at_rest()
    good = r["vmax"] <= 1e-8 and r["void_stress"] == 0.0 and r["dm"] <= 1e-10
    ok &= good
    print(f"  largest speed {r['vmax']:.1e} m/s (<= 1e-8); largest stress in empty cells "
          f"{r['void_stress']:.1e} Pa (0); partly filled cells al = {r['partial']:.3f}; "
          f"mass change {r['dm']:.0e} (<= 1e-10): {P(good)}")

    print("\n[2] Block of 10 x 10 cells translating through void at 5 m/s for 50 cells")
    print("    (elastic, no gravity, 10 mm cells)")
    r = translating_block()
    good = (r["dm"] <= 1e-10 and r["smax"] <= 1e-6 and r["spread"] <= 3.0
            and abs(r["dcx"]) <= 1e-3)
    ok &= good
    print(f"  {r['steps']} steps; mass change {r['dm']:.0e} (<= 1e-10); largest |stress| "
          f"{r['smax']:.1e} Pa (<= 1e-6); material at most {r['spread']:.2f} cells beyond "
          f"the exact block (<= 3); partly filled cells per row {r['partial']}; centroid "
          f"travel error {100 * r['dcx']:.3f} % (<= 0.1 %), vx {r['vrange'][0]:.6f} to "
          f"{r['vrange'][1]:.6f} m/s: {P(good)}")

    print("\n[3] Elastic block (10 x 10 cells) falling 0.1 m in void onto a rigid floor")
    r = falling_block()
    good_fall = r["e_acc"] <= 0.005 and r["smax_fall"] <= 1e-6
    good_land = r["dm"] <= 1e-10 and r["ymin"] >= 0.05 - 0.01 and r["fb_max"] > r["weight"]
    ok &= good_fall and good_land
    print(f"  free fall: centroid acceleration {r['acc']:.5f} m/s^2, {100 * r['e_acc']:.3f} % "
          f"from g (<= 0.5 %); largest |stress| while falling {r['smax_fall']:.1e} Pa: "
          f"{P(good_fall)}")
    print(f"  landing (impact at {r['v_impact']:.2f} m/s): lowest centroid {r['ymin']:.4f} m "
          f"(>= 0.04 m: half the block less one cell; the floor is at 0); peak floor force "
          f"{r['fb_max']:.0f} N/m (weight {r['weight']:.0f} N/m); centroid vy after "
          f"{r['box'].t:.3f} s {r['vy_end']:+.3f} m/s; mass change {r['dm']:.0e} "
          f"(<= 1e-10): {P(good_land)}")

    print("\n[4] Separation: cohesionless block (T0) pulled apart at +-0.5 m/s, no gravity")
    r = separation()
    e_gap = abs(r["void_cells"] / r["expect_gap"] - 1.0)
    good = (r["dm"] <= 1e-10 and r["smax"] <= 1e-6 and r["gap"] >= 1 and r["dp"] <= 1e-10
            and e_gap <= 0.05)
    ok &= good
    print(f"  void opened inside the block {r['void_cells']:.3f} cells per row (2 v t = "
          f"{r['expect_gap']:.3f}, {100 * e_gap:.1f} % off, <= 5 %), {r['gap']} fully empty "
          f"column(s); largest normal stress {r['smax']:.1e} Pa "
          f"(no tension, <= 1e-6); momentum change {r['dp']:.0e}; mass change "
          f"{r['dm']:.0e}: {P(good)}")

    print(f"\nWP 3.2p checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return ok


def collapse_check(cells_short=10):
    """Column collapse against decision 0019 (WP 3.2q done-when, decision 0021)."""
    ok = True
    t_all = time.perf_counter()
    P = lambda good: "PASS" if good else "FAIL"            # noqa: E731
    print("[5] Granular column collapse, T0 Mohr-Coulomb (c = 0, phi 35, psi 0; E 1 MPa),")
    print(f"    smooth left wall, floor friction {COLLAPSE_DELTA_B:g} deg, from K0 rest")
    for a, lo_hi in ((0.5, (0.30, 0.90)), (2.0, (1.2, 3.6))):
        r = column_collapse(a, cells_short=cells_short)
        fit = 1.2 * a
        g_run = lo_hi[0] <= r["runout"] <= lo_hi[1]
        g_top = abs(r["H_f"] - r["H_i"]) <= r["cell"] if a < 1.0 else True
        g_slope = r["slope_max"] <= 37.0
        good = g_run and g_top and g_slope and r["dm"] <= 1e-10
        ok &= good
        print(f"  a = {a}: L_i {r['L_i']:.3f} m, H_i {r['H_i']:.3f} m, {r['nx']} x {r['ny']} "
              f"cells of {1000 * r['cell']:.1f} mm; {r['steps']} steps to t = {r['t']:.3f} s "
              f"({r['seconds']:.0f} s), KE {r['ke_end']:.0e} of peak")
        print(f"    run-out (L_f - L_i)/L_i = {r['runout']:.3f} (fit 1.2 a = {fit:.2f}, band "
              f"{lo_hi[0]} to {lo_hi[1]}); H_f = {r['H_f']:.4f} m against H_i {r['H_i']:.3f} "
              f"m{' (within 1 cell required)' if a < 1.0 else ''}; steepest deposit slope "
              f"{r['slope_max']:.1f} deg (<= 37; over two cells {r['slope2_max']:.1f}); "
              f"mass change {r['dm']:.0e}: {P(good)}")

    print(f"\nColumn collapse checks: {'ALL PASS' if ok else 'FAILED'} "
          f"({time.perf_counter() - t_all:.0f} s)")
    return ok


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--verify", action="store_true", help="run the WP 3.2p checks")
    p.add_argument("--collapse", action="store_true",
                   help="run the column collapse checks (WP 3.2q, tests/manual)")
    p.add_argument("--case", choices=("lake", "block", "fall", "sep", "c05", "c2"))
    args = p.parse_args(argv)
    if args.verify:
        return 0 if verify() else 1
    if args.collapse:
        return 0 if collapse_check() else 1
    if args.case:
        f = {"lake": lake_at_rest, "block": translating_block, "fall": falling_block,
             "sep": separation, "c05": lambda: column_collapse(0.5),
             "c2": lambda: column_collapse(2.0)}[args.case]
        r = f()
        print({k: v for k, v in r.items() if k not in ("box", "hc")})
        return 0
    p.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
