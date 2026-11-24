# 0128 Effective stress and the gas carried by the grains (WP F.3)

Date: 2026-11-24. Status: accepted (follows decisions 0125 items 4, 5, 8, 9 and 0127).

## Context

F.2 solved the pore pressure with the solid frozen. F.3 couples the two both ways
(decision 0125 item 5): the constitutive models see the effective stress, the momentum
sees the total stress, and the solid's motion and volume change drive the gas. Sia asked
that every test use `ergunWenYu`, that 2 and 4 processors match serial to 1e-10 with gas
mass to 1e-12, that the cost with gas on be reported against gas off (F.2 measured up to
3.6x against the 1.5x F.9 target), and that F.3 bring gas validation pack test 3, the
pore-pressure decay in a suddenly loaded fine-powder column against Terzaghi.

## Decision

1. **sigma is the effective stress.** eulerianSolid's `sigma` (what the constitutive
   models update, the Riemann solver and the walls see) is sigma' once gas is on; nothing
   in the constitutive code changes. The gas pressure enters the momentum as a force on
   each cell, -sum(p_f S_f) over its open faces (linear face values, boundary values on
   patches, the immersed wall face with the cell value), so sigma' - p I acts in total and
   a uniform pressure exerts no force (a column at rest stays at rest to round-off). Gas
   inertia and gas weight are left out (Stage A). In a partly filled cell the force is
   scaled by alpha (F.4 will use it; F.3 runs without free surface). `eulerianSolid::setPorePressure(p)` switches it on;
   without gas nothing changes (v1 bit for bit). The pressure is the last solve's (one
   step behind, like any explicit force); its stiffness p_abs/e (0.25 MPa) is not added
   to the CFL step, which stays the v1 step.
2. **Gas carried with the grains.** d(phi_g rho_g)/dt + div(phi_g rho_g u_s) + div(rho_g q)
   = 0 is discretised to match the solid's remap exactly: the gas volume crossing a face
   in a step is the Lagrangian volume dt phiL less the grain volume dm/rhoP that the remap
   moved across it (`eulerianSolid::remapMass()`, kept only with gas), at the donor
   cell's gas density of the start of the step; the change of phi_g in each cell comes
   from the solid's new density. A compressed cell therefore squeezes its gas by exactly
   its Lagrangian volume change (undrained p0 dV/(phi_g V)), and the gas mass changes
   only by the boundary fluxes. Non-coupled boundaries carry no mixture volume in the
   full-cell remap; outlets, inlets and the void are F.4 and F.5.
3. **Surface loads** (needed for the test, useful later): `loads { <patch> { pressure
   <Pa>; startTime <s>; } }` in materialProperties puts a normal pressure on a traction
   patch, which moves as a free surface does under it.
4. **Two-way is the default**: `coupled yes` in gasProperties; `coupled no` keeps F.2's
   frozen solid (the F.1 and F.2 tests and pack test 1 now say so: F.1 checks that the
   fractions leave the solid byte for byte, and an unconfined bed with air blown through
   it would otherwise be lifted). Until F.4 (void at the boundary pressure), `coupled
   yes` with a `freeSurface` solid stops with an error: gas in void cells has no
   mobility, so material moving in would trap and compress it.
5. **Cost.** F.2's inner solve to 0.1 of the residual gained one digit per outer
   iteration (8 to 10 outer iterations per step, each with the drag law on every face).
   Now: inner relTol 1e-3 (default), the first guess extrapolated from the last two
   steps, and the iteration stops when the error left (the next correction, estimated from
   the ratio of the last two once they shrink tenfold) is below pTolerance. Sharp steps
   take 3 solves instead of 8 to 10; converged values are unchanged within pTolerance.

## Consequences

- `tests/test_consolidation.sh` (about 1 min) passes: column 1 m, 50 cells, fine quartz
  d 0.1 mm, e 0.40, E 1 MPa, drained top loaded with 2 kPa: U within +0.77 / -0.02 /
  -0.23 / -0.21 / -0.11 % of Terzaghi at T = 0.05 / 0.1 / 0.2 / 0.5 / 1 (limit 3 %),
  profiles within 0.011 p_i; the steps are those without gas (0.07 % more of them, as the
  column compresses later); a column at rest under gravity and the load stays at p0 to
  2e-17; 2 and 4 processors match serial to 5e-13 of p_i (p) and 4e-11 (U), gas mass
  1e-15 on every count.
- Validation pack test 3, `tests/validation/gas/terzaghiColumn` (`Allrun`, about 1.5 min):
  80-cell column to T = 1, U within 0.66 % of Terzaghi from T = 0.05 on 1, 2 and 4
  processors, which agree to 1e-12 of p_i, gas mass 1e-15; `results.md`, `times.csv`.
- Parallel caveat (solid, not gas): a column 4 cells wide cut along its length (2 x 2)
  parts from serial by 1e-9 of p_i, because the solid itself does: with gas off its
  stress parts by 4e-6 of the load there (round-off lateral motion, the serial/parallel
  issue of decision 0080). Cuts across the column match to 1e-12. The tests use those.
- Before T = 0.05 (0.36 s for this column) the elastic waves of the sudden load, which
  the quasi-static theory leaves out, are still being damped; the velocity damping
  (50 1/s) delays consolidation by about rho D c/M = 0.8 % (seen: 0.2 %).
- Cost per step, gas on against gas off: 1.4 to 1.5x on a 20 x 100 column through the sharp
  start (2.6 solves per step); 2.5x on the 50-cell 1D column, where the per-step overhead
  of the gas solve outweighs the tiny solid step; F.2's sealed blob 2.7 to 2.8x (was 3.6x). The
  1.5x target of F.9 holds for meshes of a few thousand cells, not for tiny ones.
- Not yet: the gas pressure on immersed bodies (`bodyForces` is the solid's part only),
  outlets, inlets, void and free surface (F.4, F.5); restart of the extrapolated first
  guess (a restarted run converges to the same p within pTolerance, not bit for bit:
  F.9).
