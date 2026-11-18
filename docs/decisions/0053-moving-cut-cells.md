# 0053 Moving bodies: remap onto the volume the body leaves, then re-cut (WP 4.3b)

Date: 2026-11-18. Status: accepted.

## Context

Decisions 0042, 0047 and 0051 couple fixed STL bodies to the material as cut-cell walls.
A body with prescribed motion (WP 4.1) must push the material, and the cut cells must
follow it, without losing mass or letting material into the body. The wall face of the
Lagrange step already moves with the body velocity (0047), so the question is how the
cell geometry changes between steps.

## Decision

1. **Order of a step.** Lagrange step at the old pose; bodies moved; remap; cut cells
   recomputed at the new pose (`bodyCells::update`) and the material carried onto them
   (`eulerianSolid::moveImmersedWalls`, called through celCoupling
   `moveImmersedWalls`). Drivers call one function after the remap.
2. **Remap target.** With moving bodies the remap maps each wall cell onto the open
   volume the body leaves, `Vt = Vo + dt wallU.A` (A into the body; at least 1e-6 Vo),
   instead of `Vo`. The Lagrangian volume already holds the wall's swept volume (its
   volume flux), so a full cell stays full and no material is squeezed into the part
   the body takes. Mixing waits for the new geometry.
3. **Carry-over.** Per cell, mass, material volume, momentum and the mass-weighted
   stress and history are kept and spread over the new open volume. A cell the body now
   covers hands its content to the cell it was mixed with, else to the face neighbour
   with the largest open volume (stepping through covered cells, up to 3 steps). Members
   of a mixing group may hold a small negative amount after the remap (overdrawn, 0051);
   it is handed over all the same (dropping it created mass). Group members are not
   clipped before mixing. A cell the body uncovers starts empty (void), so moving
   bodies need `freeSurface yes` (FatalError otherwise). Mass a covered cell cannot hand
   over is counted (`lostMass`); none was lost in the tests.
4. **Reaction output.** celFoam writes the force of the material on each body every step
   to `postProcessing/bodyForces.dat` and prints it at write times with the hand-over
   count and lost mass.
5. **Wall velocity of an open cell whose face lies on a body** (no covering body of its
   own) is taken from that neighbour; before, it was zero, which is harmless only for
   fixed bodies.

## Consequences

- `tests/test_stl_push.sh` (about 1 min): footing pushed 4.25 cm into a T0 layer at
  0.25 m/s and pulled out. Mass to 1e-14, nothing lost, no gap under the footing,
  largest material density 1.008 rho0, momentum balance of the material against body,
  wall and gravity forces to 1e-6 of the footing impulse, force 0.01 % of its peak after
  separating.
- The test uses the T0 rule set (mu(I), flowImpedance, normalCap). With plain
  Mohr-Coulomb (no flowImpedance, so no nearly-full snap, 0044) the footing's corner
  squeezes a nearly full surface cell (alpha 0.9995) to 1.39 rho0: compression there
  closes void without stress. Same family as the 1.07 rho0 of the WP 4.3c falling block
  with fixed walls. Use the T0 rule set with moving bodies.
- Body speed does not limit the time step; at 0.25 m/s it moves 0.001 cell per step.
  A fast body (more than about 0.3 cell per step) would need a limit.
- Not done: free bodies (motion from the force), mixing and hand-over across processors,
  moving bodies without void.
