# WP 3.2r: Prototype column collapse deposit slope (split from WP 3.2q)

State: done (staircase); toe split to WP 3.2s
Branch: thread branch

## Goal

The column collapse deposit in `prototype/void2d.py` has no surface slope steeper than
37 deg (decision 0019), with the WP 3.2q run-outs kept in their bands.

## Read first

- `docs/decisions/0023-plastic-face-shear-cap.md` (consequences), `0021` (rules 2 and 6),
  `0019`; `docs/wp/WP-3.2q.md` log

## Build

- The 1D slab rule (0021 rule 2) leans a partly filled top cell against its fuller
  neighbour in the x sweep, so thin top layers cannot slide sideways and the surface is a
  staircase of whole cells (45 deg per step). Try a slab position from a 2D normal
  (Youngs-type gradient of al, or "spread when al changes more across the sweep than
  along it"). The simple across/along test smoothed the surface without the shear cap
  but blew up with it (2026-11-14): find why (likely capped slip velocities in nearly
  empty cells that now exchange material sideways) before anything else.
- Then the blunt front (2 to 4 cells high).
- If 37 deg per cell cannot be met at 10 cells across H, show it with a 20-cell run and
  write the measurement question (per cell, two cells, fitted) as a decision for Sia.
  Do not change the band or the measurement silently.

## Done when

`tests/manual/test_proto_collapse.sh` passes and is moved to `tests/` (CI runs it):
a = 0.5 keeps its top within 1 cell; run-out 0.30 to 0.90 (a = 0.5) and 1.2 to 3.6
(a = 2); no deposit slope steeper than 37 deg; mass to 1e-10. Existing tests unchanged.

## Out of scope

- The normal jump term's residual E dependence (0023), unless it turns out to set the
  slope; C++; the outlet (3.3p).

## Log

### 2026-11-14

- Rule 9 (decision 0024): in the x sweep, layers (al >= 0.1, al gradient larger across
  x than along it) are spread instead of leaning on the fuller neighbour. The body of
  the deposit is smooth: steepest per-cell slope without the toe 27.9 deg (a = 0.5) and
  33.9 (a = 2), was 45 and 59. Run-outs 0.40 and 1.8 (20 cells: 0.45 and 2.45); top
  kept. `tests/test_proto_void.sh` passes unchanged.
- Why the 2D rule blew up in 3.2q: spreading in the y sweep (vertical films and pores),
  with or without the shear cap. Thin films (al < 0.1) spread in x sprayed along the
  floor at 20 cells (run-out 1.08), hence the 0.1 limit. Fixed a NaN path in the remap
  slopes of near-empty cells.
- Not passed: the toe. Per-cell slope at the front 62.3 deg (a = 0.5) and 48.0 (a = 2)
  at 10 cells, 46.6 and 39.8 at 20. Cause: rule 4 (void first) leaves partly filled
  cells without stress, so surface layers have no friction, creep downhill and pile up
  behind the front. Two ways of straining them were unstable or sprayed (0024). Split
  to WP 3.2s; the test stays in `tests/manual/`.
