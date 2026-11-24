# 0127 Pore pressure in a rigid bed and the drag-law selector (WP F.2)

Date: 2026-11-24. Status: accepted (follows decisions 0125 and 0126).

## Context

F.2 adds the gas pressure to `celFluid` with the solid frozen (decision 0125 item 4):
the isothermal gas mass balance with a Darcy-Ergun flux, solved implicitly each step,
plus the run-time drag-law selector with `ergunWenYu` as the default and for every test
(item 7). Sia also asked that it be fully parallel (2 and 4 processors match serial to
1e-10, gas mass to 1e-12), that the cost against gas off be reported, and that it bring
gas validation pack test 1, the packed-bed pressure drop against Ergun (item 9).

## Decision

1. **Equation.** d(phi_g rho_g)/dt + div(rho_g q) = 0, q = -lambda grad p,
   rho_g = (p0 + p)/(R T). The mobility lambda = e^2/K(e, u) comes from the drag law with
   the slip u that satisfies |grad p| = K(e, u) u / e, so Ergun's inertial term is part
   of the flux law, not a correction. No gas weight (12 Pa/m, small against bed weights
   of 10 to 20 kPa/m) and no grain velocity (F.3).
2. **p is the gauge pressure** (absolute minus p0). With absolute pressure the solve's
   round-off (1e-16 of 1e5 Pa) was larger than the convergence asked for and the
   iteration stalled; the gas mass sums the p0 part and the gauge part apart for the same
   reason.
3. **Faces.** Each face's coefficient is the harmonic mean of the two cells' mobilities
   (each cell's gas fraction e_m and k, so a side without permeability shuts the face),
   both taken at the face pressure gradient (interpolated cell gradient with the normal
   part from the two cell values, the corrected snGrad), times the linear face gas
   density, times the open face fraction with STL bodies. The Laplacian is "Gauss linear
   corrected" by default, so non-orthogonal hex cells are treated as everywhere else
   (Sia, 2026-11-19). Boundaries from `0/p`: fixedValue (fixed pressure), zeroGradient
   (sealed), and fixedGradient patches named in gasProperties `inflow` get the gradient
   that carries a given superficial velocity at p0 (a fixed mass flux). Without `0/p`,
   p = 0 and every patch is sealed.
4. **Solution per step.** Newton on the correction: the residual of the mass balance with
   the flux law (secant mobility) at the current p, a matrix with the tangent mobility
   d|q|/d|grad p| (Picard matrix on the first iteration and whenever the last correction
   grew), PCG/DIC to 0.1 of its initial residual (inexact Newton), Anderson acceleration
   of depth 5, until the largest correction is below `pTolerance` (default 1e-14 p0, about
   1e-9 Pa). The step then ends with p = p_old - dt div(F)/S from the face fluxes F of
   the converged state, so the gas mass changes exactly by the boundary fluxes. Plain
   Picard needed 15 to 50 solves per step with strong inertia, Newton alone overshot from
   a sharp start, and the non-orthogonal correction (explicit in the Laplacian) converges
   only linearly; Anderson fixed both (17 to 50 solves down to 6 to 17 on the hard cases,
   2 to 4 on smooth ones).
5. **Residual by hand.** `fvMatrix::residual()` adds the processor-face terms on top of
   the interface update, which made the Newton correction diverge on 2 processors. The
   residual is assembled as `fvMatrix::solve` does (boundary diagonal, sources of the
   non-coupled patches, interfaces through `lduMatrix::residual`), and the correction is
   solved with `lduMatrix::solver` on the same interfaces.
6. **Drag selector.** `gasDragModel` (run-time selection, `drag <name>;`, coefficients in
   `<name>Coeffs`) gives K(e, u, d, mu, rho_g), the slip at a given gradient (default:
   bracketed Illinois to 1e-14), dK/du (default central difference), mobility and tangent
   mobility. `ergunWenYu` has closed forms on the Ergun branch (slip from the quadratic,
   dK/du exact) and writes the Wen-Yu branch with Cd Re, so u = 0 gives its Stokes limit.
   F.1's permeability now comes from the selected law (unchanged for ergunWenYu, 7e-16
   against `tools/gas_reference.py`). F.7 adds the other laws.

## Consequences

- `tests/test_pore_pressure.sh` (about 40 s) passes: permeameter (0.2 m, d 1 mm,
  e 0.40) 162.746 / 2487.215 Pa against compressible Ergun 162.746 / 2487.217 Pa at 0.05 /
  0.5 m/s (1e-6); Wen-Yu branch (e 0.85) 27.7652 Pa against the integrated law; step at a
  face against erfc 0.074 / 0.030 / 0.007 % of the step at 2.5 / 5 / 10 ms; sealed
  27-degree non-orthogonal box with a 20 kPa blob: gas mass to round-off, ends uniform at
  the mass-conserving pressure to 9e-13 of the blob; 2 and 4 processors match serial to
  3e-11 (p and qGas, every write), gas mass 7e-16 on every count.
- Validation pack test 1, `tests/validation/gas/ergunBed` (`Allrun`, about 1 min): sand
  d 0.5 mm, e 0.40, column 0.3 m (10 x 60), 0.02 to 0.2 m/s, all within 0.008 % of
  Ergun on 1, 2 and 4 processors, which agree to 5e-12; `times.csv`, `results.md`. At 600
  cells the 2 and 4 processor runs are no faster than serial (too small to scale).
- Cost per step, gas on against gas off (one core): 1.09 on the settled F.1 column (1,000
  cells, pressure uniform, 1 solve per step); 3.6 on a 100 x 100 sealed box with a 20 kPa
  blob (8.7 solves per step at the start, 1 at the end). The F.9 target (1.5) is met only
  when the pressure changes slowly per step; sharp starts cost several solves. Options
  if F.9 needs more: looser `pTolerance` (1e-12 p0), a better preconditioner (GAMG), or
  coefficients reused while they barely change.
- Parallel agreement depends on `pTolerance`: the worst case (3e-11 of a 160 Pa drop) has
  a margin of 3 to the 1e-10 asked for.
- Not yet: grain velocity and volume change (F.3), force on the solid (F.3), void at the
  boundary pressure and partly filled cells (F.4; void cells have no mobility now and
  keep their pressure), gas weight.
