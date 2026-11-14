# 0012 2D prototype: dimensional splitting and the oblique-step overshoot (WP 1.3a)

Date: 2026-11-14. Status: accepted for the prototype; the 45 degree overshoot goes to
the formulation decision in WP 1.3b.

Numbers 0009 to 0011 are taken by WP 2.5a, 2.1 and 2.5b.

## Context

WP 1.3a extends the cell-centred 1D scheme (decision 0005) to 2D plane strain on a
fixed Cartesian mesh (`prototype/lr2d.py`). The brief suggested direction splitting with
the 1D Lagrange step plus remap in each sweep.

## Decision

1. **Dimensional splitting, symmetric.** Each step runs x-then-y and y-then-x from the
   same state and averages density, momentum and stress (velocity as the mass-weighted
   mean). This costs two sweeps more per step than alternating the order, but a problem
   mirrored about the diagonal gives the mirrored result bit-for-bit (alternating does
   not). `split="alternate"` stays as an option.
2. **Per-sweep physics.** In a sweep along n: P system (v_n, s_nn) with `rho c_p`,
   S system (v_t, s_nt) with `rho c_s`, each with MC-limited characteristic slopes and an
   acoustic Riemann solver; s_tt and s_zz change only with the normal strain
   (`lambda` times it). Shear strain increment `dt (v_t*_R - v_t*_L)/h_mid`.
   Stresses remapped volume-weighted, velocities mass-weighted, as in 1D.
3. **Remap and Lagrangian density written as increments.** `lr1d.py` now computes
   `rho - rho dh/h_new` instead of `m/h_new`, and the remap as increments of the
   Lagrangian values (algebraically identical). With that, a sweep across a direction
   where nothing varies is an exact identity, and a wave along x in 2D matches `lr1d.py`
   bit-for-bit. Before the change it matched only to 2e-9: the strain `(h_new - h)/h`
   cancels badly when the material streams (h_new - h is 1e-6 of h), so a one-ulp
   difference anywhere grows to 1e-9. The 1D tests give the same orders, speeds and
   widths as before; only round-off-level numbers move (step extremum 1.81e-6 to
   1.78e-6 Pa, mass and momentum change 0 to 8e-16).
4. **45 degree step P-wave: reported, not a pass criterion.** The step P-wave along x
   (material streaming diagonally) has no new extrema (0.000 percent). The same step at
   45 degrees overshoots: up to about 14 percent of the jump in s_xy and 8 percent in
   s_xx during the first 20 steps, settling to 3 to 5 percent. It does not depend on the
   time step (CFL 0.5, 0.25 and 0.1 give the same), on the limiter (first order also
   overshoots) or on the streaming (it is there at rest), and it is the same in cells at
   100 and 200 cells per metre. So it is not a splitting error: it is the
   direction-by-direction Riemann solver. Seen along x, a 45 degree P-wave carries a
   large left-going S part (`s_xy + z_s v_y` jumps by -0.32 G eps), which only cancels
   against the y-direction fluxes in the limit; a 1D Riemann solver per face cannot
   see that. The brief's step criterion ("material streaming diagonally") is read as
   the wave along x; the 45 degree number is printed by the test for the record.

## Consequences

- Smooth waves are fine at every angle: orders 1.77 to 1.86, speeds within 0.002
  percent, and the 45 degree errors are smaller than the axis-aligned ones (the pulse
  covers more cells per width along the diagonal).
- A shock or contact front oblique to the mesh will ring by a few percent with this
  scheme. In CEL most fronts are oblique. WP 1.3b must weigh this in the formulation
  decision. Known fixes: a multidimensional (nodal) Riemann solver as in cell-centred
  Lagrangian schemes (GLACE, EUCCLHYD; Maire's elastic-plastic version), or the
  staggered layout of Abaqus/LS-DYNA (decision 0005 measured its own ringing).
- Changing `lr1d.py` arithmetic changed no WP 1.1 or WP 1.2 result beyond round-off.
