# 0017 C++ remap: geometry from face fluxes, Barth-Jespersen form of MC and superbee

Date: 2026-11-14. Status: accepted (WP 2.4; Sia's standing rule: proceed and merge on
green CI). Numbers 0015 and 0016 are taken by WP 2.5c and the silo re-plan (WP 2.3 also drafted a 0016).

## Context

WP 2.4 fills `eulerianSolid::remap(dt)` against the contract of decision 0014: the
remap gets only the face volume flux `phiL` and the Lagrangian cell state. It must be
second order, conservative and monotone on any `fvMesh`, serial and parallel, reduce
to the 1D prototype remap on a uniform mesh, and take the limiters of decision 0013
item 5 per field. Brief: `docs/wp/WP-2.4.md`. Numbers below come from
`tests/test_remap.sh` (`tests/eulerian/remap_test.py`, `--study` for the Courant scan).

## Decision

1. **Swept-region geometry from the fluxes alone (contract 0014 unchanged).** A cell
   velocity `uC` is the least-squares fit of the face fluxes, `sum_f (uC.S_f -
   phiL_f)^2 -> min` (exact for a translation; the mean of the two face velocities in
   1D). The Lagrangian centroid is `C + dt uC`. The face displacement is the
   tangential part of the mean `uC` of the two cells plus `phiL/|S|` along the normal,
   times `dt`; the swept region's centroid is `Cf + d/2`. Using the full displacement
   vector (not only its normal part) gives the cross-derivative terms of a
   second-order unsplit scheme. No corner transport: the whole swept volume comes from
   the donor cell. A corner-transport remap would still need no change of the
   contract.
2. **Reconstruction.** Per field component, a least-squares gradient (weights
   `1/|d|^2`) on the Lagrangian centroids, evaluated at the swept centroid of the
   donor. Density by volume flux; velocity, stress and history per unit mass, fluxed
   with the mass flux, reconstructed in space (the 1D prototype reconstructs them in
   mass coordinates; the two agree exactly when density is uniform and differ at second
   order otherwise). Mirror ghost cells (zero difference) across wall, symmetry and
   empty-direction boundaries regularise the fit; they never limit.
3. **Limiter: one scalar factor per cell, Barth-Jespersen form.** For each face,
   `t = bound / (g.d_PN)`, where `bound` is the largest (`g.d > 0`) or smallest
   (`g.d < 0`) of `0` and `q_N - q` over **all** neighbours of the cell. Then MC
   `min(1, 2 t_min)`, superbee `max(min(2 t_min, t_max), min(t_min, 2 t_max))`, minmod
   `min(1, t_min)`, upwind `0`. On a uniform 1D mesh `t` is the one-sided slope over the
   central slope and these are exactly the limiters of `lr1d.py` and `lr1d_plastic.py`.
   A first version used each face's own difference `q_N - q` instead of the neighbour
   range; then round-off between identical rows (1e-13) made a cell a "local
   extremum" in the transverse direction and switched off its gradient along the
   motion, and the rows of the x-translation case drifted apart by 3 percent of the
   jump within 20 steps. With the neighbour range the rows agree to 4e-15. The factor
   is also capped so the value at every swept centroid the cell donates to stays in
   the neighbour range (never active in 1D).
4. **Limiters are chosen in `constant/materialProperties`**, optional subdictionary
   `remapLimiters { rho MC; U MC; sigma MC; epsP superbee; epsPEq superbee; }` (these
   are the defaults). `system/fvSchemes` was not used: these are not OpenFOAM schemes,
   and the material file is where the other remap-relevant choices live.
5. **Only internal and coupled faces carry flux.** `phiL` on wall, symmetry and empty
   faces is ignored. Processor and cyclic faces get neighbour values, gradients and
   cell velocities through temporary volFields, so each side computes the same flux to
   round-off.
6. **Update written as increments of the Lagrangian state**,
   `rho = rho_L + (rho_L (V_L - V) - sum dm)/V`, `q = q - sum dm (q_f - q)/m_new`, so a
   cell with no flux through any face keeps every value bit for bit.

## Measured (tests/test_remap.sh)

| Check | Result | Criterion |
| --- | --- | --- |
| x translation on 64x3x3, Courant 0.4, against `lr1d`/`lr1d_plastic` | 7e-14 | 1e-10 |
| Order, 2D diagonal 32/64/128 (rho, U, sigma), Courant 0.25 per direction | 2.75, 2.73, 2.73 | 1.7 |
| Order, 3D body diagonal 12/24/48 | 2.35, 2.34, 2.34 | 1.7 |
| Square, diagonal: new extrema in rho, U, sigma | 4.5e-15 of the range | 1e-10 |
| History edge width after 50 cells in x and y (71 cells of travel) | 2 cells | 4 |
| Slotted disc, one turn, 100 x 100: new extrema / L1 shape error | 0.00 % / 0.127 | 1 % / reported |
| Uniform state, swirl flux, 32 steps | 5.5e-16 | 1e-12 |
| Uniform U, sigma, history, compress then reverse | 0 | 1e-12 |
| Square diagonal on 2 / 4 processors vs serial | 1.8e-14 / 2.5e-14 | 1e-12 |
| Mass, momentum, sum(m q), every case | at most 1e-13 | 1e-12 |

Stability (square, diagonal, `--study`; material Courant number per direction):
2D clean up to 0.45 (sum 0.9), new extrema 0.6 percent at 0.5 (sum 1.0), negative
density at 0.6. 3D clean up to 0.30 (sum 0.9), 53 percent overshoot at 0.35 (sum 1.04),
negative density at 0.4. So the limit is a summed outflow Courant number of about 0.9,
as expected without corner transport. The acoustic time step (0.3 c) keeps material
Courant numbers far below that unless the material moves at a large fraction of the
wave speed.

Density under a compressing flux followed by its reverse does not come back to uniform
exactly (9.8e-4 relative at Courant 0.2 on 40 x 40): after the first step the density
really is non-uniform, and no conservative remap undoes that to round-off. The brief's
round-off requirement is therefore checked for velocity, stress and history (exact
there) and density is reported.

## Consequences

- WP 2.4b runs Lagrange step plus remap together; the remap needs nothing more from
  WP 2.3 than `phiL` on internal and coupled faces and `rho_` set to the Lagrangian
  density.
- Each remap computes 17 component gradients with two neighbour exchanges each; on the
  48^3 case (110 k cells) a step of the test program costs about 0.27 s serial,
  diagnostics included. Fine for now; batching the exchanges is
  the first optimisation if profiles show it.
- If a later test needs a larger Courant number (supersonic material speeds, contact),
  add corner transport (WP 2.4b or later), still from `phiL` and the geometry above.
