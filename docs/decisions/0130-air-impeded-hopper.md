# 0130 Air-impeded hopper discharge (WP F.5)

Date: 2026-11-24. Status: accepted (follows decisions 0125 items 7-9, 0128, 0129).

## Context

F.5 closes Stage A: with the pore gas on, coarse grains should discharge as in v1 and
fine grains should slow down as published. The brief names Crewdson, Ormond and
Nedderman (1977) as the reference. Sia's rules: ergunWenYu drag in every test; 2 and 4
processors match serial to 1e-10 with gas mass to 1e-12; report the cost with gas on
against off; validation pack test 4 on 1, 2 and 4 processors; aligned, body-fitted
meshes; coarse where possible; never ask Sia for data.

## Decision

1. **Reference.** Crewdson et al. (1977, Powder Technol. 16, 197) could not be read: the
   publisher, ResearchGate and the open repositories refused automated access, and
   Muite et al. (2004)'s PDF on the Caltech repository is blocked by the network
   policy. What could be read: an open thesis (Mamtani 2011, Univ. Florida) states
   Crewdson's model: the air pressure gradient at the orifice acts as an extra body
   force on top of the weight, with the gradient from Carman-Kozeny; the abstract of
   Hsiau, Hsu and Smid (2010, Phys. Fluids 22, 043306) says standard correlations
   over-predict below about 500 microns; Muite's abstract says the fluid matters when
   N = U_t / U_out is below 10. So the check is the published direction and size
   (coarse unchanged, fine slower, large below 0.5 mm) plus Crewdson's relation
   evaluated with the run's own outlet gradient, not a point-by-point comparison.
   Crewdson's measured rates are still owed (F.9 or whenever the paper is readable).
2. **Case.** Plane-flow (slot) hopper, full width, 2D: slot 60 mm, walls 30 deg from
   vertical, bin 0.2 m wide, fill 0.27 m; two blockMesh blocks, the hopper columns fan
   from the outlet edge to the wall top, rows graded 4:1 to the outlet; 544 cells,
   outlet cells 3.75 x 3.1 mm. Quartz grains, gas fraction 0.40, Mohr-Coulomb phi 30,
   c 0, psi 1, E 1 MPa, fullBlend; walls 20 deg; damping 50 1/s until the outlet opens
   at 0.15 s; rate by least squares from 0.3 to 0.5 s. Each run takes about a minute,
   so no long runs and no checkpoints were needed. The Couto hopper as a control was
   not rerun: its 4 mm wheat relaxes the gas in about 1 ms, and the 4 mm case of the
   sweep shows the gas changes the rate by 0.05 %. Muite's N for the Couto hopper
   (wheat 4 mm, 1380 kg/m^3, U_t 10.8 m/s, outlet speed 1.43 m/s at 118 kg/s) is 7.5:
   below Muite's 10, yet the sweep's 1 mm case (N 8.2) changes by only 2 %, so for air
   Muite's threshold is conservative.
3. **Mesh first tried.** A 20 mm slot with 2 mm outlet cells discharged at 1/15 of
   Beverloo in v1 even without gas, and more slowly on a finer mesh (2.2 instead of 7
   kg/(m s)); a 60 mm slot at 3.75 mm cells gives Beverloo. Recorded as a v1 issue,
   not pursued here.
4. **Rules.** The sweep uses v1's T0 rule set (mu(I) with the case's d, shear cap,
   rule 9, flow impedance, normal cap): dry rate 74.1 kg/(m s), Beverloo 72.6. The
   parallel checks use plain Mohr-Coulomb: with the T0 rules the solid alone (gas off)
   parts serial from 2 processors by 4 % of U by 0.25 s (decision 0080), with plain
   Mohr-Coulomb and gas on they agree to 4e-12.

## Results

`tests/test_air_impeded.sh` (about 8 min on 4 cores), T0 rules, gas on against off on
the same mesh, rate in kg/(m s):

| d [mm] | dry | gas | change | Muite N | Crewdson sqrt(1 + G/(rho_b g)) |
| --- | --- | --- | --- | --- | --- |
| 4 | 71.86 | 71.90 | +0.1 % | 20.2 | 0.999 |
| 1 | 74.13 | 72.47 | -2.2 % | 8.2 | 0.983 |
| 0.5 | 73.46 | 69.52 | -5.4 % | 4.9 | 0.933 |
| 0.2 | 74.40 | 54.58 | -26.6 % | 2.4 | 0.568 |
| 0.1 | 74.55 | 26.80 | -64.1 % | 1.9 | 0 |

(also 0.05 mm by hand: -87 %). 1 mm and up within 5 % of dry (limit 5 %); the rate
falls monotonically with d; the slowdown starts near 0.5 mm, as Hsiau et al. report,
and the gas matters where Muite's N is below 10. Crewdson's relation with the gradient
over the outlet row matches at 1 mm and 0.5 mm (0.983 against 0.978, 0.933 against
0.946) and over-predicts the slowdown at 0.2 and 0.1 mm (the half-cell gradient at the
orifice is not resolved: it exceeds the weight at 0.1 mm). Gas mass 5e-15; the outlet
holds p0.

Parallel (plain Mohr-Coulomb, 0.1 mm, gas on, to 0.25 s with material leaving, cuts
2 x 1 and 2 x 2 through the outlet, the bed and the free surface): p within 3.6e-13 of
the largest |p|, U within 4.1e-12, gas mass 1.4e-15.

Validation pack test 4 (`tests/validation/gas/fineHopper`, about 2 min), plain
Mohr-Coulomb, d 0.1 mm, to 0.5 s: gas 16.92 against dry 32.94 kg/(m s) (-49 %); dry on
1, 2 and 4 processors agree to 3e-12; gas on 2 processors to 7e-11. **Gas on 4
processors (2 x 2) fails the 1e-10 rule**: from 0.356 s one partly filled surface cell
by the left wall (the top of the draining bed) differs by up to 8e-10 of the largest
U; full cells stay within 6e-11, discharged mass within 2e-13, gas mass 2e-15. A restart
from the 0.35 s write shows the round-off of the pressure solve (1e-13 relative) is
amplified in nearly empty surface cells (alpha 0.04) as they pass; without gas the same
cells do not part. Fix in F.9 (which reruns the pack).

Cost per step, gas on against off, 0.1 mm grains (6 to 8 solves per step: stiff
diffusion at fine grains): 1.60x with the T0 rules, 2.2 to 2.4x plain on 1 processor,
2.9x on 2 and 3.2x on 4 (544 cells, too small to scale). Above the 1.5x target; the
fixes of decision 0129 (one mobility call per face, no non-orthogonal correction where
not needed, one assembly per iteration) are F.9's.

## Consequences

- `tests/gas/hopper_gas_cases.py` builds the hopper (also `draw_mesh` for the picture);
  `tests/test_air_impeded.sh` is F.5's done-when test; pack test 4 is
  `tests/validation/gas/fineHopper` (Allrun, `times.csv`, `results.md`).
- No solver change in F.5.
- Open: Crewdson's measured rates (paper not readable); pack test 4 on 4 processors
  (one surface cell, 8e-10); the narrow-slot v1 discharge (item 3); cost (F.9).
