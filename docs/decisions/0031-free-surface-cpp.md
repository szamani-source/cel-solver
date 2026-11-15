# 0031 Material plus void in the C++ cycle (WP 3.2a)

Date: 2026-11-15. Status: accepted (WP 3.2a). Follows decision 0021 (prototype rules 1
to 7) face by face in the unsplit C++ step.

## Context

The C++ silo track (decision 0029) needs void before the column collapse (3.2b), the
outlet (3.3) and the discharge (3.4) can run in celFoam. The prototype applies its void
rules per 1D sweep; the C++ Lagrange step and remap are unsplit and face based
(decisions 0017, 0018). Decision 0021 suggested isoAdvector for the material volume
flux, but rule 6 ties the flux to the Lagrangian cell lengths and the slab position,
which isoAdvector does not know.

## Decision

1. **Switch.** `materialProperties` `freeSurface yes;` turns the void on; field `alpha`
   (read if present, default 1) is the material volume fraction. Without the switch the
   code takes the old branches, bit for bit (every earlier test unchanged).
2. **Cartesian hex meshes only.** Every face must be normal to an axis (checked at the
   start, fatal otherwise). The prototype's per-sweep slab rule then becomes a per-axis
   rule on each face: no isoAdvector. Revisit only if a non-Cartesian mesh is needed.
3. **Slab side per axis** (rule 2): each cell leans along each axis towards the fuller
   of its two neighbours on that axis (difference above 1e-9), spread if equal or full.
   Walls, symmetry planes and prescribed-velocity patches mirror the cell, traction-free
   patches count as empty, coupled patches give the neighbour. A cell touches a face
   unless it is empty or partly filled and leans away. Rule 9 (top layers spread
   sideways) is not in: it belongs to the collapse (3.2b).
4. **Homogenised cells** (rule 3): face states and forces use alpha sigma, impedances
   alpha z of the material (density rho / alpha bounded to 0.5 to 2 rho0), mass rho V.
   A cell is second order only if it and every face neighbour are full; otherwise its
   face states are the cell values, and under gravity it keeps 1 - phi of its
   hydrostatic jump (rule 7, phi from the stress along g against half its own weight).
   No stress predictor in first-order cells; empty cells do not move.
5. **Faces** (rule 4): contact faces (both sides touch) get the Riemann solution, unless
   the material has no tensile strength (`c 0` or `sigmaT 0`) and the solution pulls;
   all others are traction-free with the one-sided velocity of the side that touches
   (owner first), else of the side(s) holding material, else zero. Boundary faces the
   material does not touch: zero velocity and traction.
6. **Void closure** (rule 5) per axis of the cell: compression closes void first;
   expansion unloads the material to zero normal stress, then opens void (a full cell of
   a material with tensile strength stretches). Shear goes to the material. A full cell
   that opens no void keeps the strain of the full scheme exactly.
7. **Remap** (rule 6): material volume flux through a face is the donor's material in the
   swept slice along the face axis, from its Lagrangian length, Lagrangian fraction and
   slab side; a donor left with under 1e-9 of its volume gives all of it; fluxes under
   1e-9 of a cell are dropped. Mass flux: MUSCL density times the swept volume where the
   donor is second order, else volume flux times the donor's Lagrangian material
   density. Velocity, stress and history mass-weighted, first order unless the donor is
   second order. Dust (alpha below 1e-9 with mass) takes its volume from its mass, with
   no velocity or stress. Material crossing a traction-free boundary leaves; nothing
   enters from outside.
8. **Time step**: the material's wave speed (density rho / alpha bounded as in 4).

## Consequences

- Measured (`tests/test_void.sh`, 17 s, 10 mm cells): lake at rest 2.5e-17 m/s over
  1000 steps (prototype 2.5e-16); block translating 50 cells with edges exact (0.5 cell,
  no partly filled cell), stress 4e-10 Pa, mass 6e-15; free fall at g to 1e-6 relative,
  stress 1e-10 Pa, landing at 0.049 m (no loss); cohesionless block cracks with void
  2.000 cells per row for 2 v t = 2.000, no tension; 2 processors match serial to 3e-16
  (translating across the split) and 2e-15 (crack on the split).
- The slab rule is per axis, so a tilted surface is a staircase (as in the prototype).
- Partly filled cells carry no stress in expansion and are frictionless (the prototype's
  rules); the collapse package 3.2b adds the face shear cap (decision 0023, rule 8) and
  rule 9 and measures the run-out at E = 0.1, 1, 10 MPa with neither rule 12 nor 13.
- The outlet (rule 10, 10b) is WP 3.3. A traction-free patch next to void already lets
  material leave.
