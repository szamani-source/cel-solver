# WP 3.2s: Prototype collapse toe: stress in partly filled cells (split from WP 3.2r)

State: not started
Branch: thread branch

## Goal

Partly filled surface cells in `prototype/void2d.py` carry stress (and so frictional
strength), the deposit stops creeping, and the column collapse toe is no steeper than
37 deg per cell, with the WP 3.2q run-outs and the WP 3.2r smooth surface kept.

## Read first

- `docs/decisions/0024-top-layers-spread-sideways.md` (what was found, consequences),
  `0021` rules 2 to 4, `0023`, `0019`; `docs/wp/WP-3.2r.md` log

## Build

- Rule 4 (void first) leaves a partly filled cell stress-free until it is full, so a
  thin surface layer has no friction, slides below the friction angle and piles up
  behind the front (a = 0.5 still creeps at t = 2 s with 10 cells). Find a closure
  that lets a slab with a free side carry its own weight without the stiffness of a
  slab of length al h (wave transit al h / c below the time step). Ideas not yet tried:
  give the material the strain of the cell's homogenised length (dh / h) only in the
  gravity sweep and only for slabs lying on material; or carry the overburden of the
  cell's own weight as a stress floor in the yield check of the face shear cap.
- Rejected in 3.2r (0024): all compression to the slab (unstable), a share al of it
  where the free side faces void (films sprayed along the floor, a = 2 run-out 3.3).
- Check the creep first: a = 0.5 should stop (kinetic energy below 1e-3 of peak)
  before t = 1 s at 10 cells.
- If the toe still fails at 10 and 20 cells, write the measurement question (exclude
  the front cell, or two-cell slope) as a decision for Sia. Do not change the band or
  the measurement silently.

## Done when

`tests/manual/test_proto_collapse.sh` passes and is moved to `tests/` (CI runs it):
a = 0.5 keeps its top within 1 cell; run-out 0.30 to 0.90 (a = 0.5) and 1.2 to 3.6
(a = 2); no deposit slope steeper than 37 deg; mass to 1e-10. Existing tests unchanged.

## Out of scope

- The normal jump term's E dependence (0023); C++; the outlet (3.3p).

## Log

(none yet)
