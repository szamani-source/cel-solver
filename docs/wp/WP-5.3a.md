# WP 5.3a: Two-way coupling of a deformable body

State: done (2026-11-18). Needs WP 4.3b, 5.1. Split from WP 5.3 (soft bodies to 5.3b).
Branch: the thread's `claude/project-thread-*` branch

## Goal

A deformable hex8 body (WP 5.1) as a cut-cell wall of the material: its surface moves
with its elements, the material pushes back on its nodes.

## Read first

- `docs/decisions/0053`, `0070`, `0072`; `src/celLagrangian/rigidBody/rigidBody.H`
  (deformable entry); `src/celCoupling/immersedWalls/immersedWalls.C` (moveBodies)

## Build

- `deformable { min; max; cells; density; material; planeStrain; damping; drive; }` in a
  body entry of `constant/rigidBodies`: a box of hex8 elements, surface = its boundary.
- celCoupling `moveBodies`: wall-cell forces to the nearest surface triangle's nodes,
  summed over processors; the elements run one coupling step ahead (subcycled) so the
  surface velocity the material sees is exactly the surface's displacement over the
  step; bodyCells rebuilds the search of a deformable surface each step.
- hex8Body `damping` (viscous stress).

## Done when

`tests/test_deformable_push.sh`: the footing push of WP 4.3b rigid and as a steel hex8
body driven at its top: mass to 1e-10, nothing lost, density within 2 %, steel mean
force over the push within 5 % of rigid, both separate on pull-out.

## Out of scope

- Soft bodies (WP 5.3b), body meshes other than a box, contact between bodies.

## Log

- 2026-11-18: Done (decision 0072). Steel mean force 39.2 N against rigid 40.4 N
  (-3 %), density 1.013 rho0. Learned: an undamped steel body rings (node velocities
  0.7 m/s around a 0.25 m/s drive) and the material, seeing that as wall velocity, packs
  to 1.6 rho0; the elements must also run a step ahead, or the surface moves other than
  the material saw (1.18 rho0). Soft body (E 20 MPa): shortening 1.44 times the elastic
  column and density 1.62 rho0 under the footing: WP 5.3b.
