# 0005 Cell-centred Godunov Lagrange step in the 1D prototype

Date: 2026-11-14. Status: accepted (Sia approved the switch on 2026-11-14);
revisit in WP 1.3.

## Context

The WP 1.1 brief asked for a staggered grid (velocity at nodes, stress at cells)
with a central-difference Lagrange step, and for a step wave with no new extrema.
That version was built first (commit "WP 1.1: staggered-grid prototype", kept in git
history). It met the smooth-wave order (1.84, 1.88) and conservation, but the step wave
always rang:

| Artificial viscosity tried | Step overshoot | Smooth order |
| --- | --- | --- |
| Linear, unlimited (alpha 1 to 2) | below 1e-5 | about 0.7 (first-order damping) |
| Christensen-type limited (MC, minmod, van Leer) | 4 to 7 percent | 1.4 to 1.96 |
| Smoothness sensor (JST-type), tuned | 1e-4 to 3e-3 | 1.2 to 1.9, never both |

No setting met both criteria. A central-difference step is dispersive, and a viscosity
strong enough to stop the ringing also damps smooth waves at first order.

## Decision

Keep every field at the cell centre and use a second-order Godunov Lagrange step:
MUSCL-Hancock with limited slopes in characteristic variables, and an acoustic Riemann
solver at each face for the face velocity and stress. The remap stays a conservative
MUSCL swept-region remap. Both use the MC (monotonized central) limiter, because the
harmonic van Leer limiter clipped smooth peaks enough to keep the order at 1.73 to 1.78.

Result (`tests/test_proto_1d_elastic.sh`, material streaming at 0.25 c): orders 1.88 and
1.83; largest new extremum 2e-7 of the jump, which is round-off and the bar's small
nonlinearity, not ringing; front speed error 0.04 percent; mass and momentum change
below 1e-15.

## Consequences

- The 1D prototype matches the layout of the OpenFOAM solver (cell-centred FV,
  CLAUDE.md), so WP 2.4 can port it more directly.
- No artificial viscosity is needed for elastic waves. Shocks in later packages may
  still need a check.
- Abaqus CEL uses a node-based (staggered) layout. WP 1.3 decides the 2D/3D
  formulation and can revisit this choice with the numbers above.
- The smooth-wave order is set by the limiter clipping the wave peak. It stays at about
  1.8 and drops slowly on finer meshes and at faster streaming (1.78 at 0.5 c). It
  passes at the test setting only. Recorded as a known issue in `docs/STATUS.md`.
