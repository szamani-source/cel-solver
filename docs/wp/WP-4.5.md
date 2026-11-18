# WP 4.5: Silo and hopper discharge with STL walls

State: not started. Needs WP 4.3a.
Branch: the thread's `claude/project-thread-*` branch

## Goal

The T0 silo discharge with its walls from STL, then a hopper (inclined STL walls) that
the box mesh cannot follow.

## Read first

- `docs/decisions/0033`, `0041`, WP 4.3a's log; `docs/spec/silo-reference.md`

## Build

- Tutorial `tutorials/stlHopper` (STL walls, outlet still a mesh patch).

## Done when

`tests/test_stl_silo.sh` passes: flat-bottom tank with STL side walls, discharge rate
within 5 % of the mesh-boundary tank at the same settings, no leakage through walls;
a 60 deg hopper runs to empty, rate recorded against a hopper reference (Beverloo with
the hopper correction, source opened) without a pass band until Sia agrees one.

## Out of scope

- Moving inserts (WP 4.3b).

## Log

(none yet)
