# 0025 Column collapse toe: why it stays steep, and how to measure it

Date: 2026-11-15. Status: accepted. Sia chose option 1, the fitted flank angle
(2026-11-15); the collapse test uses it and runs in CI. The prototype's behaviour is unchanged (two study knobs added,
both off).

## Context

After WP 3.2r (decision 0024) the column collapse deposit is smooth below 34 deg except
at the toe: per-cell slope at the front 62.3 deg (a = 0.5) and 48.0 deg (a = 2) at 10
cells across the short side, against 37 (decision 0019). The WP 3.2s brief blamed
partly filled cells carrying no stress (rule 4 of 0021), so that surface layers have no
friction and creep into a stuck front.

## What was found

- **Partly filled cells are frictionless**, as suspected: rule 4 keeps their normal
  stress at zero, and Mohr-Coulomb with c = 0 then returns any shear to zero.
- **Giving them strength makes the toe worse.** `FLOOR` (rule 10, study knob): a
  supported partly filled cell is returned to the yield surface with the overburden
  pressure of its column added (not stored, no force on neighbours). Stable, but the
  surface layers stick: a = 2 run-out 1.8 -> 1.5, toe 48.0 -> 59.9 deg, body
  33.9 -> 37.4. At a = 0.5 the per-cell slope read 29.1 only because the last column fell
  just under one cell (0.86) and so left the measurement; the deposit still crept.
- **The a = 0.5 creep is at faces held at yield, not in partly filled cells.** At
  t = 0.6 to 2 s full cells in the second row slide over the first at about 1 cm/s near
  the front (the first row is held by floor friction) and the front column grows a cell
  at a time. Turning the shear cap (rule 8, decision 0023) off at t = 0.6 s stops the
  deposit at once (KE 6e-4 of peak). Near the front the material is at yield in its
  normal stresses alone (active state), so every face there is capped and slips with no
  resistance; nothing rate-dependent is left to stop it.
- `VSTICK` (rule 11, study knob): a capped face sticks while the Riemann slip speed is
  below VSTICK. It stops the creep (a = 0.5 stops at 0.9 s with 0.003 m/s) but the toe
  slope it gives is noise, not a trend:

  | VSTICK (m/s) | 0 | 0.001 | 0.002 | 0.003 | 0.005 | 0.01 | 0.03 |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | a = 0.5 toe (deg) | 62.3 | 30.3 | 32.6 | 38.3 | 47.7 | 49.5 | 53.5 |
  | a = 0.5 run-out | 0.40 | 0.40 | 0.35 | 0.35 | 0.30 | 0.30 | 0.25 |
  | a = 2 toe (deg) | 48.0 | 33.2 | 37.5 | 36.3 | 38.5 | 40.8 | 42.8 |
  | a = 2 run-out | 1.8 | 1.8 | 1.9 | 1.8 | 1.7 | 1.9 | 1.7 |

  The two "passes" at 0.001 are runs whose last column ended just under one cell
  (0.92, 0.82) and so dropped out of the check. A velocity threshold also brings a
  scale that depends on E and the mesh. Not adopted.
- **The toe is a resolution limit.** With the creep stopped (cap off at 0.6 s) the toe
  still reads 39.7 deg: the flowing layer is one to two cells thick, and when it stops
  its nose is a step of one to two cells. Per-cell toe slope falls with refinement
  (62.3 -> 46.6 deg at a = 0.5, 48.0 -> 39.8 at a = 2, 10 -> 20 cells, decision 0024) but
  slowly, and the 20-cell runs take 10 to 20 minutes each.

## Measurements of the same 10-cell deposits (no physics change)

| | a = 0.5 | a = 2 |
| --- | --- | --- |
| Per cell, all columns (current check) | 62.3 | 48.0 |
| Per cell, without the last pair (toe) | 27.9 | 33.9 |
| Over two cells | 45.6 | 33.4 |
| Line fitted from 10 to 90 % of H_f | 22.2 | 22.3 |

At 20 cells (WP 3.2r deposits, same code): per cell 46.6 / 39.8, without the last pair
37.8 / 39.8 (the steepest single-cell step is no longer always the last one), over two
cells 34.4 / 28.6, fitted 24.9 / 17.8. Single-cell values move by up to 30 deg between
runs that look alike; the fitted angle moves by a few degrees between meshes.

The fitted angle is printed by `--collapse` from this package on, for information.

## Options put to Sia

1. Fitted flank angle from 10 to 90 % of H_f, 37 deg or less (recommended): passes at
   10 and 20 cells, the test moves to CI. It does not see a single steep step; the
   staircase of 3.2r would still show as a fitted angle near 45 deg only if it covered
   the flank.
2. Exclude the last pair of columns from the per-cell check: passes at 10 cells, fails
   at 20 (37.8, 39.8), so refining the mesh would break the test.
3. Keep the check as it is and add physics (rate-dependent friction such as mu(I)); the
   test stays manual.
(A two-cell baseline was also considered: a = 0.5 fails at 10 cells, 45.6.) A first
card recommending option 2 was withdrawn when the 20-cell numbers came in.

## Consequences

- The a = 0.5 deposit still creeps at about 1 cm/s at t = 2 s (KE 2e-3 of peak). 3.3p
  (discharge) has the same capped faces at yield next to its free surface; a creeping
  free surface over the outlet is expected. A rate-dependent friction (mu(I)) is the
  physical cure and would also replace `VSTICK`.
- `FLOOR` and `VSTICK` stay in `prototype/void2d.py` as study knobs, off by default,
  so the runs above can be repeated; default results are unchanged bit for bit.

## Decision (Sia, 2026-11-15)

The collapse test judges the deposit by a line fitted to the surface from 10 to 90 % of
the final height, at most 37 deg. Single-cell slopes (all columns, without the toe, over
two cells) are printed for information only. `tests/test_proto_collapse.sh` is in
`tests/` and runs in CI: 22.2 and 22.3 deg at 10 cells, run-outs 0.40 and 1.8, top kept,
mass 5e-14. Revisit if a later package needs a local check on single steep steps.
