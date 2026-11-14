# WP 3.2p: Prototype free surface (material plus void)

State: done (2026-11-14); column collapse split to WP 3.2q
Branch: `wp-3.2p-void`

## Goal

The 2D prototype carries one material plus void: cells can be full, partly filled or
empty, the free surface moves with the material and stays sharp, and empty space carries
no stress. A decision note fixes the void treatment for the C++ code.

## Read first

- `docs/decisions/0013`, `0014`, `0016`; `docs/wp/WP-3.1p.md` log; `prototype/silo2d.py`
- `docs/spec/silo-reference.md` (column collapse values), if it exists
- Benson (1992), sections on multi-material cells and interface tracking

## Build

- A volume fraction `alpha` per cell (material volume over cell volume), advected in the
  remap with a sharp scheme. Try in this order and keep the first that passes: the
  existing swept-region remap with a compressive limiter (superbee); a 1D-per-sweep
  geometric (PLIC or SLIC) reconstruction. Record numbers for each tried.
- Partly filled cells (the void closure). Candidates to test, with the numbers recorded
  in the decision note: (a) all volumetric strain of a partly filled cell goes to the
  void until the cell is full, deviatoric strain goes to the material; (b) cell stress =
  `alpha` times material stress, strain shared equally; (c) a threshold on `alpha` below
  which a cell is treated as void. State variables are per unit material mass; empty
  cells hold no state and are skipped.
- Free surface in the Lagrange step: zero traction where material meets void; no
  velocity or stress "leaks" into empty cells; a rule for the velocity of newly filled
  cells (momentum-conserving).
- Tension: with tension cut-off 0 the material separates instead of pulling; separated
  material must be able to open void inside the domain.
- `tests/test_proto_void.sh` and the decision note
  `docs/decisions/NNNN-void-treatment.md`, including what the C++ version (WP 3.2)
  should do and whether isoAdvector fits (it advects with a face flux, which `phiL`
  provides).

## Done when

`tests/test_proto_void.sh` passes:

- Lake at rest with a free surface inside the domain (flat, and tilted mesh-diagonally
  is not required): velocities below 1e-8 m/s for 1000 steps; no stress above the surface.
- A block of material translating through void at constant velocity for 50 cells keeps
  its mass to round-off, its edges within 3 cells, and its stress at zero.
- A block falling freely under gravity accelerates at g within 0.5 percent and lands on
  a rigid floor without losing mass or passing through.
- 2D granular column collapse, aspect ratios 0.5 and 2, cohesionless Mohr-Coulomb,
  friction 35 degrees, tolerances from `docs/spec/silo-reference.md` section 6
  (decision 0019): a = 0.5 keeps its top (H_f = H_i within 1 cell); run-out within 0.5
  to 1.5 times the quasi-2D fit 1.2 a (0.30 to 0.90 at a = 0.5, 1.2 to 3.6 at a = 2);
  no deposit slope steeper than 37 degrees.
- Mass conserved to 1e-10 in all cases. Existing tests still pass unchanged.

## Out of scope

- Several materials in a cell, the outlet (3.3p), C++.

## Log

### 2026-11-14 (thread "Add free surface to tank prototype")

Done except the column collapse, which is split to WP 3.2q (`docs/wp/WP-3.2q.md`,
recorded in `docs/work-packages.md`). `prototype/void2d.py` (class `VoidBox`, sweep
`void_sweep`, built on `silo2d.SiloBox`; `silo2d.py` unchanged so the 3.1p test is
unchanged), `tests/test_proto_void.sh` passes in about 70 s. Decision 0021. Branch was
the thread branch `claude/project-thread-8wfjkq`, started from the 3.1p branch (PR #17,
open at the time) merged with main.

- Lake at rest, 0.1 x 0.6 m box filled to 0.455 m (top cell half full), K0, 1000
  steps: 2.5e-16 m/s, zero stress in empty cells, mass exact.
- Block 10 x 10 cells at 5 m/s for 50 cells (1830 steps): mass 2e-14, edges exact, no
  partly filled cells left, stress 3.5e-9 Pa, centroid exact.
- Elastic block falling 0.1 m: acceleration g to 1e-13, stress 2.6e-9 Pa while falling;
  lands and bounces (vy +0.22 m/s after 0.17 s), centroid never below 0.049 m, mass 5e-14.
- Cohesionless block pulled apart at +-0.5 m/s: void 2.000 cells per row = 2 v t, one
  fully empty column, no tension, momentum 1e-15.
- Tried and dropped, with numbers in decision 0021: walls as full for the slab side
  (lake at rest 6e-4 m/s); material stress and impedance in partly filled cells
  (blow-up at al 1e-3); s_hat uniform in all first-order cells (falling block 40 kPa,
  0.1 % of its mass shed); no separation at pulling faces (the pulled block stretched
  over 6 cells instead of cracking); second order next to partly filled cells with a
  "holds material" mask (falling block blew up).
- Column collapse (T0, 10 cells across, floor 35 deg, E 1 MPa): a = 0.5 run-out 0.25
  (band 0.30 to 0.90; 20 cells also 0.25), top kept; a = 2 run-out 1.2 (band 1.2 to
  3.6); toes 50 to 60 deg per cell (limit 37), main slopes 25 to 33 deg. Run-out grows as
  E falls (a = 0.5: 0.35 at E 0.1 MPa), so the Riemann solver's dissipation is part of it;
  glass-bead friction (22 deg) gives 0.67 to 0.71 of the experimental fit. About 1 to 2
  minutes per case at 10 cells. Details and the reading in decision 0021.

Next: WP 3.2q (column collapse) and WP 3.3p (discharge) can both start; 3.3p should
measure the discharge rate's dependence on E first (same dissipation).
