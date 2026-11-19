# WP 2.6: Taylor bar on the Eulerian wedge

State: done
Branch: `claude/project-thread-53vkqx`

## Goal

The Phase 2 benchmark: the copper Taylor bar (32.4 mm, 227 m/s) on the Eulerian solid,
run axisymmetric on the wedge (3.8a/b).

## Read first

- `docs/decisions/0063-eulerian-taylor-bar.md`, `docs/decisions/0061-axisymmetric-wedge.md`
- `docs/decisions/0071-viscous-hourglass.md` (the Lagrangian Taylor bar, WP 5.2)

## Build

- No solver changes: `vonMises` with `freeSurface yes` on a wedge, rod filled with
  -227 m/s, smooth floor, no gravity (celVoidCheck).

## Done when

- `tests/test_taylor_wedge.sh` passes (about 7.5 min): at 0.2 mm, shortening within
  5 % and widening within 20 % of the axisymmetric Abaqus benchmark (13.11 / 6.006
  mm); errors fall from 0.4 to 0.2 mm; mass kept to 1e-10; bar at rest.

## Out of scope

- Johnson-Cook hardening; finer runs (0.1 mm, about 50 min).
