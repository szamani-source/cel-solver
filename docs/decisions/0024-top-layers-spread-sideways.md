# 0024 Top layers spread sideways in the x sweep (deposit staircase)

Date: 2026-11-14. Status: accepted for the prototype (WP 3.2r). The column collapse
toe is not fixed and moves to WP 3.2s.

## Context

After WP 3.2q (decision 0023) the column collapse run-outs pass, but the deposit is a
staircase of whole cells: per-cell slopes 45 deg (a = 0.5) and 59 deg (a = 2) against
37 (decision 0019). Rule 2 of decision 0021 puts the material of a partly filled cell
as one slab against its fuller neighbour in each sweep. A thin top layer of a pile has
a fuller neighbour uphill, so in the x sweep it is modelled as a vertical slab with a
gap on the downhill side: it cannot move sideways until its cell is full, and column
heights snap to whole cells. WP 3.2q found that a 2D rule smoothed the surface without
the shear cap but blew up with it.

## What was found

- The blow-up is not the shear cap. Spreading partly filled cells in the **y sweep**
  (vertical films and pores, |d al/dx| > |d al/dy|) blew up with and without the cap
  (a = 2: kinetic energy 5 times the baseline, single cells at 3 m/s against the left
  wall). The x sweep alone is stable. Restricting the y rule to cells next to void and
  standing on material stopped the blow-up but left the toe as it was (a = 2, 20 cells:
  37.9 deg; a = 0.5, 20 cells: 58.6 deg).
- Spreading every top layer in the x sweep let films a few percent of a cell thick
  shoot along the floor: at 20 cells (a = 0.5) 1.2 cells of material reached the far
  wall and the run-out read 1.08 (band 0.30 to 0.90). Layers thinner than 0.1 of a
  cell therefore keep rule 2.
- The steep part left is the toe, the last one or two columns. Its cause is deeper:
  rule 4 of 0021 (compression closes void first) means a partly filled cell carries no
  stress until it is full. A thin surface layer therefore has no frictional strength,
  slides downhill below the friction angle and piles up behind a front held by floor
  friction (a = 0.5 still creeps at 1 to 2 cm/s at t = 1 s; it stops at 0.5 s with 20
  cells). Two ways of letting such cells strain were tried and rejected: all of the
  compression to the slab (unstable: a slab of length al h has a wave transit time
  al h / c, below the time step, and a 99.9 percent full cell kept a free face inside
  the bulk), and a share al of it in cells whose free side faces an empty cell (stable,
  but films sprayed over the floor, a = 2 run-out 3.3).
- A NaN path: the remap slopes divide by cell masses, which can be round-off small next
  to the surface; inf times the zero mask gave NaN in the flux (overflow warnings in
  every collapse run since rule 2 cells got thinner). Slopes are now set to zero where
  the mask is zero instead of multiplied by it.

## Decision

Rule 9 of the void prototype (`SLAB2D`, `SLAB2D_MIN` in `prototype/void2d.py`): in the
x sweep, a partly filled cell with al >= 0.1 whose al gradient (Youngs 1-2-1 stencil) is
larger across the sweep than along it is spread evenly (lean 0) instead of leaning on
its fuller neighbour. The y sweep keeps rule 2.

## Result (`tests/manual/test_proto_collapse.sh`, E = 1 MPa, steepest slope per cell,
two-cell slope in brackets, body = without the last pair of columns)

| Case | Before (0023) | This decision | Band |
| --- | --- | --- | --- |
| a = 0.5, 10 cells: run-out, top | 0.40, 0.06 cell low | 0.40, 0.05 cell low | 0.30 to 0.90, 1 cell |
| a = 0.5, 10 cells: slope (body) | 45.0 (34.8) (45) | 62.3 (45.6) (27.9) | 37 |
| a = 2, 10 cells: run-out | 1.6 | 1.8 | 1.2 to 3.6 |
| a = 2, 10 cells: slope (body) | 59.4 (53.3) | 48.0 (33.4) (33.9) | 37 |
| a = 0.5, 20 cells: run-out; slope | 0.43; 45.0 (27.6) | 0.45; 46.6 (34.4) | |
| a = 2, 20 cells: run-out; slope | 2.05; 51.7 (48.5) | 2.45; 39.8 (28.6) | |
| Mass | 5e-14 | 5e-14 | 1e-10 |

The surface of the deposit is now smooth below 34 deg everywhere except the toe; at 20
cells the a = 2 run-out moves towards the fit (2.45 against 2.4). The collapse test
still fails, on the toe only, and stays in `tests/manual/`. `tests/test_proto_void.sh`
passes unchanged. The bands and the measurement were not changed.

## Consequences

- WP 3.2s: give partly filled cells a stress (frictional strength) without the
  stiffness problem, then re-measure the toe. The per-cell slope at the toe moves by
  15 deg between variants that look alike (33 to 62 deg at a = 0.5), so it measures the
  last cell more than the deposit. If 3.2s cannot make the toe pass, the question for
  Sia is whether the 0019 slope check should exclude the front cell or use a two-cell
  baseline; that choice is not made here.
- 3.3p (discharge) has the same stress-free partly filled cells at the free surface and
  the jet; expect the free surface over the outlet to creep.
- C++ WP 3.2 follows rules 1 to 9; rule 9 needs the volume fraction gradient across the
  sweep, which an unsplit C++ scheme has anyway (isoAdvector's interface normal).
