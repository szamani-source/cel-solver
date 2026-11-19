# WP 5.3b: Soft deformable bodies

State: done (2026-11-18). Needs WP 5.3a. Split from WP 5.3.
Branch: the thread's `claude/project-thread-*` branch

## Goal

A soft hex8 body (stiffness of the order of the material's) coupled both ways without
packing the material.

## Read first

- `docs/decisions/0072`, `0073`; `src/celLagrangian/rigidBody/rigidBody.C`
  (addSurfaceForce, nearestFace, faceSamples); `src/celCoupling/immersedWalls`
  (moveBodies); `src/celCoupling/bodyCells/bodyCells.C` (wall velocity)

## Done when

`tests/test_soft_push.sh`: a soft footing (E 20 MPa, nu 0.3) pushed into T0 shortens as
the elastic column within 20 % and leaves the material density within 2 %.

## Log

- 2026-11-18: Done (decision 0073). Cause: material loads on one triangle of each
  element face stirred the hourglass modes; the bottom zigzagged and the material
  packed (1.62 rho0). Loads now go to the element face's 4 nodes, spread over the
  cell's wall face, and the material sees the mean surface velocity there; soft bodies
  use `hourglass 0.3`. Density 1.0155 rho0, shortening ratio 0.90. Steel case improved
  to -0.5 % of rigid. Tried and dropped: impedance dashpot, heavier body damping.
