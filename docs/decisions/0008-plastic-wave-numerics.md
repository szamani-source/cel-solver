# 0008 Numerics for elastic-plastic waves and history remap (WP 1.2)

Date: 2026-11-14. Status: accepted. Sia chose option 1 of item 3 (check after the
start-up) on 2026-11-14.

Numbers 0006 and 0007 are taken by WP 2.2 and WP 0.4 branches running in parallel.

## Context

WP 1.2 adds plasticity to the WP 1.1 bar (`prototype/lr1d_plastic.py`). Three numerical
questions came up while making `tests/test_proto_1d_plastic.sh` pass.

## Decision

1. **Godunov step unchanged; elastic wave speeds.** The acoustic Riemann solver, the
   characteristic slopes and the time step use the elastic impedance
   `rho c_e`, `c_e = sqrt((K + 4G/3)/rho)`, as the brief asked. No separate two-wave
   Riemann solver. Wave speeds and the HEL come out within 0.2 percent at 1000 cells.
   The price: across a plastic wave the true impedance is lower (`rho c_p`), and an
   acoustic solver with too large an impedance puts the face stress outside the two cell
   values by `(z_e/z_p - 1)/2` of their difference (13 percent here). Inside a moving wave
   the limited slopes hide this (0.1 percent extrema); at an impulsive start they do not
   (item 3).
2. **History remap uses the superbee limiter; stress keeps MC.** All stress components
   and history are remapped mass-weighted (as velocity in WP 1.1), so `sum(m q)` is exact.
   With MC, a plastic-strain edge carried at 0.25 c kept widening: 10-90 percent width 5
   cells after 60 cells of travel, 6 after 120, failing the 4-cell criterion. Superbee
   holds it at 2 cells (1-99 percent: 5 cells) and does not widen. Superbee squares off
   smooth peaks, so velocity, density and stress stay on MC. Consequence: stress and
   history of a cell are remapped with different limiters. Harmless for perfect plasticity;
   with hardening the yield stress and the stress could drift slightly apart at an edge,
   which the next stress update returns to the surface.
3. **Impact start-up is reported, not checked.** A sudden velocity jump at a face
   (the step-loading test) makes a start-up overshoot at the impact plane of 2.7 percent
   of the jump (von Mises), 2.5 percent (Drucker-Prager), and 1.7 percent for the purely
   elastic bar, so it is a property of the WP 1.1 scheme under impulsive loading, made
   worse by item 1. It is the same in cells at every mesh size, and decays to
   0.8 percent after the elastic wave has run 5 cells and 0.1 percent after 10. The test
   therefore checks the 1 percent extrema criterion after the first 10 cells of wave
   travel and prints the start-up value. Options given to Sia: this; an elastic-plastic
   Riemann solver to remove it (one or two more sessions); or ramping the impact over a
   few cells (about 1.0 percent, borderline).

## Consequences

- Known issue in `docs/STATUS.md`: impulsive start-up overshoot. Contact (WP 4.x) starts
  impulsively, so the 2D/3D solver may need either the elastic-plastic Riemann solver
  or a little artificial viscosity at contact. Revisit in WP 1.3.
- The C++ remap (WP 2.4) needs a per-field limiter choice.
- The uniaxial strain increment per step is `(h_new - h)/h_mid` (as WP 1.1); with no
  rotation in 1D no objective rate is needed. 2D needs one (WP 1.3).
