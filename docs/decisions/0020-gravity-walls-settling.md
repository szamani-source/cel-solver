# 0020 Gravity, wall boundaries and settling in the Lagrange step

Date: 2026-11-14. Status: accepted (WP 3.1p, prototype; the C++ WP 3.1 follows it).

## Context

Stage 1 of T0 (`docs/spec/target-problems.md`) is a tank settling under gravity with
smooth or Coulomb-frictional walls. A cell-centred Godunov scheme with a plain gravity
source term is not well-balanced: a column in equilibrium has a linear stress profile
whose reconstruction, Hancock predictor and boundary ghosts each produce a small
spurious force, and the column drifts. Wall friction has to fit the acoustic Riemann
solver. Settling needs damping that is off again for the discharge.

## Decision

1. **Gravity by subtracting the discrete weight.** In the sweep along gravity each
   column gets the weight profile G (G at a face = g times the mass per unit area below
   it, from the cell densities at the start of the sweep). The Lagrange step works on
   `s_hat = s_nn - G` and has no separate source term, which is exact because G jumps by
   `rho g h` across a cell. Equilibrium means `s_hat` uniform, so slopes, predictor and
   Riemann solutions are uniform and nothing moves. Measured: a K0 column stays below
   3e-18 m/s for 1000 steps. The actual face stress is `s_hat* + G(face)`.
2. **Boundaries as one-sided Riemann solutions.** Ghost cells are parity mirrors (wall:
   stress even, normal velocity odd; free: stress odd about the traction, velocity even;
   smooth wall tangential: shear odd; frictional: the sticking ghost). The boundary
   face's (s*, v*) is then set from the interior state alone with
   `s* = s_I + z_I (v* - v_I)` (sign by side):
   - wall: v_n* = 0; if that would give tension, the material separates (zero
     traction);
   - friction: keep the sticking shear if it is at most `tan(delta_w)` times the wall
     pressure from the normal solution; otherwise clip it to that limit with the
     sticking sign and take v_t* from the characteristic, so friction opposes sliding;
   - free: zero traction.
3. **No mass flux through boundary faces** in the remap. Mass is exact. At the free top
   of a completely full tank, the top cells lose some density as the column settles (7
   percent in the top cell for 1.2 m); WP 3.2p replaces this with a real free surface.
4. **Settling = gravity ramp plus velocity damping, both switched off afterwards.**
   Gravity ramps from 0 with a (1 - cos) profile over 0.1 s, about two periods of the
   column's first mode; damping `v <- v / (1 + alpha dt)` with alpha = 250 1/s, about
   critical for that mode. Settled when kinetic energy is at most 1e-6 times the
   elastic strain energy.
5. **Mohr-Coulomb vectorised** in `prototype/constitutive/batch.py` (main plane, edge,
   apex on stacked eigen-decompositions); tension cut-off returns are handed to the
   point model, so it matches `models.py` in every case.

## Consequences

- The settling has to load the walls monotonically. With a faster ramp (0.05 s) and
  lighter damping (200 1/s) the column rebounds after the ramp, the top half of the wall
  unloads and the wall ratio ends 10 to 25 percent below `tan(delta_w)`. That is a real
  dynamic effect, not an error, but Janssen assumes full mobilisation, so the checks use
  the slower ramp. Overdamping (500 1/s) also works but takes twice as many steps.
- After settling, the state creeps very slowly (speeds about 3e-4 m/s; 1.5 percent of
  `rho g H` change over a further 10 percent of the run without damping) because the
  walls sit exactly at the friction limit. The discharge (3.3p) starts from it anyway.
- Costs (10 mm cells, CFL 0.5, dt 58 us): 2520 steps and about 110 s of Python for the
  40 x 120 frictional tank. At the spec's 5 mm cells this is about 4 times the cells and
  twice the steps. The C++ code (WP 3.1) should do the same thing: G per column needs a
  column sum, which on an unstructured, decomposed mesh is better replaced by the same
  idea in face form (subtract `rho g . (x_face - x_cell)` from the reconstructed face
  stress, hydrostatic reconstruction), to be checked against the same lake-at-rest test.
