# WP 4.3d: Free bodies

State: done (2026-11-18). Needs WP 4.3b. Split from WP 4.3b (free bodies were out of scope).
Branch: the thread's `claude/project-thread-*` branch

## Goal

A rigid STL body whose motion follows from the force of the material on it and gravity
(a falling block, an anchor), instead of a prescribed path (WP 4.3b).

## Read first

- `docs/decisions/0053`; `src/celLagrangian/rigidBody/rigidBody.H`;
  `src/celCoupling/immersedWalls/immersedWalls.H`

## Build

- `mass` (and optional `inertia`) in a body's entry makes it free; `velocity` and `omega`
  are then initial values.
- The solid keeps the force through each cell's wall face, so the moment about the
  body's centre can be summed (`bodyMoment`).
- celCoupling `moveBodies(solid, bodies)` replaces `bodies.move()` in the drivers: moves
  the bodies, then changes each free body's velocity by the step's force, moment and
  gravity (central difference, velocities staggered).

## Done when

`tests/test_free_body.sh` passes: a block thrown in the void follows the central-difference
parabola to 1e-9 with constant omega; a steel block dropped onto a T0 layer comes to rest
(speed below 1 % of the impact speed), the material carries its weight within 2 %, the
momentum of block and material balances gravity and the mesh walls within 1 %, mass kept
to 1e-10 with nothing lost, material density within 2 %.

## Out of scope

- Restart of a free body from a later time; contact between bodies; bodies hitting mesh
  walls; mixing across processors (as WP 4.3b).

## Log

- 2026-11-18: Done in one session (decision 0058). Test passes in about 1 min: the block
  sinks 2.6 mm and carries its weight to 0.01 %, momentum balance 1e-5 of the weight
  impulse. Found and fixed: the Jacobi eigen solver divided by a negligible off-diagonal
  (zx 1e-20 in a K0 cell) and overflowed under FOAM_SIGFPE. Learned: with rotation free
  the dropped block is kicked to 0.5 rad/s at impact (the material's reaction is a few
  % asymmetric) and its corner then squeezes a cut cell to 1.038 rho0 (as decision 0053);
  the drop test runs without inertia. Next: rough wall.
