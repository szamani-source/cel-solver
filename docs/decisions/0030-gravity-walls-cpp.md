# 0030 Gravity and wall boundaries in the C++ Lagrange step (WP 3.1)

Date: 2026-11-15. Status: accepted (WP 3.1). Follows decision 0020 (prototype) in face
form, as 0020 suggested for an unstructured, decomposed mesh.

## Context

WP 3.2t (decision 0029) moved the question of rules 12 and 13 (normal jump cap,
dispersed void) to the C++ solver: run the discharge and the column collapse with one
code, starting with neither rule. The C++ cycle had no gravity, no walls, no void and no
outlet, so that comparison needs WP 3.1 (this), 3.2 (void) and 3.3 (outlet) first.
The prototype balances gravity with a per-column weight profile; a column sum does not
fit an unstructured mesh split across processors.

## Decision

1. **Gravity `constant/g`** (standard OpenFOAM file), ramped with (1 - cos) over
   `gravityRamp` (controlDict), evaluated at mid-step and used for both the face states
   and the momentum.
2. **Hydrostatic reconstruction per face.** For a face between cells P and Q with
   d = x_Q - x_P, each cell's hydrostatic jump is e = -rho (g . d) gHat gHat. The MUSCL
   jump becomes a - (e_P + e_Q)/2, each cell's gradient increment loses its own e, and
   the face state adds e/2 back (P side +e_P/2, Q side -e_Q/2). A column whose cell
   stresses differ by the mean e has equal face states and the face stresses balance
   each cell's weight exactly; slopes, predictor and Riemann solutions are then zero.
   First-order states get the same e/2. Non-coupled boundary faces mirror the
   deviation (ghost carries the cell's own e). Without gravity the code takes the old
   branch, bit for bit. Measured: K0 column 3e-17 m/s after 1000 steps.
3. **Walls** are named in `materialProperties` `walls { <patch> { friction <deg>; } }`
   and take precedence over the U patch type. Rigid wall at rest: v_n = 0 unless that
   would pull on the wall (then zero traction); sticking shear kept up to
   tan(friction) |t_n|, otherwise clipped to it with v_t from the interior
   characteristic (decision 0020 item 2, face form).
4. **The Hancock predictor sees the wall reaction.** At wall faces the predictor uses
   the wall Riemann solution of the face states instead of the cell value. Without it,
   the predicted velocity of a wall cell misses the friction force every step, the
   sticking shear comes out too large and is clipped, and the cells next to a
   frictional wall carried only about 55 percent of the wall shear (one-cell boundary
   layer, 0.4 x 1.2 m tank). Other boundary kinds are unchanged (their tests pass).
5. **Damping** `U <- U/(1 + damping dt)` after the momentum update (controlDict
   `damping`, `setDamping()` to switch it off once settled), as in the prototype.
6. **Wall ratio measure.** The test reads the wall shear ratio from cells extrapolated to
   the wall, as the prototype, but on the faces that slip in the settled state. The
   prototype's zone (slipped at least 10 percent of the largest slip during the run)
   also takes in faces near the bottom that slipped during the ramp and stick again;
   in C++ they sit at 0.92 to 1.0 of the limit (printed by the test).

## Consequences

- WP 3.1 results (10 mm, maxCo 0.4, dt 46 us): smooth column bottom stress 0.12
  percent from the weight, K0 exact; frictional tank Janssen within 1.7 percent with K
  0.424 at the walls (prototype 0.425, 1.7 percent), wall ratio 0.4 percent where it
  slips, force balance 0.05 percent; settles in 3580 steps, 38 s serial; 2 processors
  match serial to 7e-18 of rho g H.
- The free top of a full tank is still a traction-free boundary through which the remap
  lets material enter with the top cells' values as the column settles (mass +6e-4).
  The prototype had no flux there. WP 3.2 replaces it with void.
- Velocity, symmetry and traction boundaries still use the cell value in the predictor;
  revisit if a test shows the same one-cell layer there.
