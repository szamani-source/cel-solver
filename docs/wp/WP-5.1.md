# WP 5.1: Hex8 element with hourglass control

State: done (2026-11-18). Needs the constitutive library (WP 1.x).
Branch: the thread's `claude/project-thread-*` branch

## Goal

An explicit deformable Lagrangian body of 8-node hexahedra (`hex8Body` in celLagrangian),
using the same constitutive models as the Eulerian side, ready to be coupled to the
material in WP 5.3.

## Read first

- `docs/decisions/0070`; `src/celLagrangian/hex8Body/hex8Body.H`

## Build

- `src/celLagrangian/hex8Body`: one-point Flanagan-Belytschko hexahedron with exact mean
  gradients, stiffness hourglass control, central difference, Hughes-Winget rotation,
  lumped mass, velocity constraints per node and direction, external nodal forces,
  gravity, plane strain option, energy and work accounting.
- `tests/hex8/hex8Check` (`celHex8Check`) and `tests/test_hex8.sh`.

## Done when

`tests/test_hex8.sh` passes: MacNeal-Harder patch test (seven distorted hexahedra):
stress the same in every element to 1e-12, inner node forces zero to 1e-10; 1D strain
wave in a steel bar: stress -rho c v within 2 %, front at c t within 2 elements; free
elastic cube with hourglass modes over 2000 steps: total energy (kinetic, strain,
hourglass) within 1 %, momentum to 1e-12.

## Out of scope

- Reading a body mesh from a case, coupling to the Eulerian material (5.3), contact
  between bodies, the Taylor bar (5.2).

## Log

- 2026-11-18: Done in one session (decision 0070). Patch test to 1e-15, wave -45.98
  against -45.83 MPa, energy conserved to 1.4e-7 with 21 % of it in hourglass modes
  (the start is a checkerboard velocity). Learned: the check of energy must use the
  central-difference invariant 0.5 v(n-1/2) M v(n+1/2) + U(x_n); the plain sum of
  kinetic energy and internal work is an identity and checks nothing. Next: WP 5.2.
