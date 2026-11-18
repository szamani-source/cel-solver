# 0058 Free bodies: central difference with the step's wall force (WP 4.3d)

Date: 2026-11-18. Status: accepted.

## Context

WP 4.3b moves bodies on a prescribed path and reports the force of the material on them
(decision 0053). Many cases need the motion to follow from that force: a block dropped
onto a layer, an anchor pulled by a given force, an insert that settles.

## Decision

1. **Input.** A body with `mass` is free. `velocity` and `omega` are its initial values.
   `inertia` (symmTensor about the centre, in the time-zero pose) is optional; without it
   the body does not rotate. The mass is that of the part the mesh holds (2D: one cell
   deep), as the forces are. Gravity is the solid's `g` (with its ramp).
2. **Time integration.** Central difference with staggered velocities, as for the
   material: each step the Lagrange step uses the body's current velocity at its wall
   faces; `moveBodies` then moves the body with that same velocity (so the remap target
   of decision 0053 matches the body's real sweep) and adds the step's impulse,
   `v += dt (F/m + g)`, `w += dt J^-1 (M - w x J w)` with `J = R J0 R^T`.
3. **Force and moment.** F is minus the body force of the solid (summed over
   processors). The solid now keeps the force of each wall face (`wallForce_`); the
   moment about the body's centre uses the cut cells' centres as the points of action.
4. **Drivers.** celFoam and celVoidCheck call `moveBodies(solid, bodies)` in place of
   `bodies.move()`; prescribed bodies just move. Free bodies count as moving, so they
   need `freeSurface yes` to be coupled (as 0053); without it they fly uncoupled.
5. **Eigen solver.** Cyclic Jacobi now drops an off-diagonal that is negligible against
   its two diagonal entries (1e-22) instead of rotating by it: a K0 cell with zx 1e-20
   made theta overflow, which traps under FOAM_SIGFPE.

## Consequences

- `tests/test_free_body.sh` (about 1 min): free flight matches the scheme to 1e-14; a
  dropped steel block comes to rest on a T0 layer carrying its weight (force / weight
  1.0000), momentum of block and material balanced to 1e-5 of the weight impulse.
- With rotation free, the dropped block is kicked to about 0.5 rad/s at impact: the
  material's reaction under a symmetric block is a few % asymmetric (directional
  splitting of the remap is the likely cause; not isolated). The tilted corner then
  squeezes a cut cell to 1.038 rho0 (the corner-cell family of 0053). The drop test
  runs without inertia; a rotating free body in contact is usable but this is a known
  issue.
- The moment uses cell centres, not the cut faces' centroids: first-order in the cell
  size for the arm.
- Not done: restart of free bodies (stops with an error), body-body contact, bodies
  hitting mesh walls.
