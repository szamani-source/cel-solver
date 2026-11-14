# WP 3.2p: Prototype free surface (material plus void)

State: blocked on WP 3.1p
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
  friction 35 degrees: final run-out and height within the tolerance stated in
  `docs/spec/silo-reference.md` (if absent: within 25 percent of the published
  experimental fits, with the source cited in the log).
- Mass conserved to 1e-10 in all cases. Existing tests still pass unchanged.

## Out of scope

- Several materials in a cell, the outlet (3.3p), C++.

## Log

(none yet)
