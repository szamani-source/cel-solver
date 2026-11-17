# WP 3.6c2: stiff inclined layer (floor cell and top cell)

State: done (2026-11-17)
Branch: `claude/project-thread-qho7rn`

## Goal

Split from WP 3.6c (decision 0040): make the mu(I) inclined layer follow the Bagnold
profile at E = 10 MPa (stuck floor cell, slow body) and stop the top cell next to void
sliding over the body at 1 and 10 MPa.

## Read first

- Decisions 0039, 0040, 0044; `docs/wp/WP-3.6c.md`
- `tests/layer/inclined_layer.py`; `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C`
  (wall faces, `faceStates`), `eulerianSolidRemap.C` (`remapFreeSurface`, alpha)

## Done when

The inclined layer within 5 % of the exact solution (surface speed and profile, steady
within 1 % over the last second) at E = 0.1, 1 and 10 MPa, 5 mm cells, free surface;
`tests/test_collapse.sh`, `tests/test_outlet.sh`, `tests/test_discharge.sh` still pass.
Test: `tests/manual/test_inclined_layer.sh` (about 2.5 minutes). CI part:
`tests/test_inclined_layer.sh` (0.1 and 10 MPa, about 50 s).

## Out of scope

- The silo with mu(I) (WP 3.6d), shear-band mesh dependence (decision 0038).

## Log

- 2026-11-17: One session. Traced face velocities and tractions. Floor: the wall
  reconstruction mirrored the cell velocity, so the wall jump doubled the drag; fixed
  with a wall ghost at the wall's last velocity. Top cell: void of 1e-9 to 1e-5 from
  noise kept switching it to the partly filled rules; fixed by counting void below 1e-4
  as full. Both behind `flowImpedance yes` (decision 0044). Profile errors 3.6 / 0.9 /
  0.8 % at 0.1 / 1 / 10 MPa (was 7.8 / 22 / 136). Dropped: tilted-gravity support factor
  (no effect), one-sided traction under the top cell (unstable). Left: at 0.1 MPa the
  layer creeps faster over time (+0.7 %/s at 5 mm, +1.05 %/s at 2.5 mm). Next: WP 3.6d.
