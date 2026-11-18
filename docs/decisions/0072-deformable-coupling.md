# 0072 Deformable bodies as cut-cell walls: elements one step ahead (WP 5.3a)

Date: 2026-11-18. Status: accepted.

## Context

WP 5.3 couples the hex8 body (decision 0070) to the material both ways. The cut-cell
machinery of rigid bodies (decisions 0047, 0051, 0053) needs, per step, the surface the
material sees, a wall velocity per cut cell, and the remap target Vo + dt wallU.A. A
rigid body moves exactly with the velocity the material saw (0058). A deformable body is
stepped by its own, much smaller stable step (steel 1e-6 s against 4e-5 s for the
material) and its nodes move with many velocities within one coupling step.

## Decision

1. **Where it lives.** A body entry with `deformable { ... }` holds a hex8Body (a box of
   elements for now) inside rigidBody, so bodyCells, setImmersedWalls,
   moveImmersedWalls and the force output work unchanged. Its frame is fixed (c = c0,
   R = I); its surface (the elements' boundary, two triangles per face) moves, and
   bodyCells rebuilds its search each update.
2. **Forces.** The force of the material on each wall cell (the solid now keeps it per
   cell) goes to the nodes of the surface triangle nearest the cell centre, with the
   barycentric weights; summed over processors.
3. **Elements one coupling step ahead.** In moveBodies after the Lagrange step of step
   n: the surface moves to t(n+1) with the velocity the material saw; the elements,
   already at t(n+1), are subcycled to t(n+1) + dt(n) with step n's forces; the surface
   velocity for step n+1 is (element surface - committed surface) / dt(n). So the
   surface moves exactly as the material saw it move, whatever the elements do inside
   the step; a change of dt corrects itself the next step.
4. **Damping.** hex8Body `damping` adds a viscous stress 2 damping rho c Le D (Le = V /
   largest face area). Without it a stiff body rings for ever after each impact and the
   material, seeing the ringing as wall velocity, packs. The tests use 0.2.
5. **Drive.** Nodes inside `drive (min) (max)` move with the body's `velocity`
   Function1; without a drive box the body is free (gravity, material).

## Consequences

- `tests/test_deformable_push.sh` (about 3.5 min): steel footing as hex8 against the
  rigid STL footing: mean force -3 %, density 1.013 rho0, mass 1e-14, separates.
- Without item 3 the steel run packed the material to 1.18 rho0; without item 4 to
  1.6 rho0 (mean force -31 % at 0.1 s).
- Known issue (WP 5.3b): a soft body (E 20 MPa, the order of the material's own
  stiffness) packs the material under it to 1.62 rho0 and shortens 1.44 times the
  elastic column. Likely the remap target Vo + dt wallU.A with wallU from the cell
  centre is not the volume a bending surface sweeps; not isolated.
- The force lag of one coupling step is in the elements' dynamics, not in the geometry.
