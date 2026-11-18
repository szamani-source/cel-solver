# WP 4.5: Silo and hopper discharge with STL walls

State: done (2026-11-18). Needs WP 4.3a.
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

- 2026-11-18: Done in one session (decision 0052), after WP 4.3c (STL walls with void).
  `tests/test_stl_silo.sh` (`tests/stlwalls/stl_silo_checks.py`, about 6 min): the
  stlSilo tank (STL side walls and floor, cut cells, gap over the mesh outlet) against
  the WP 3.3 mesh-wall tank, both T0 with mu(I) at E = 0.1 MPa, three runs each (maxCo
  0.40/0.38/0.36, decision 0046 bands): mesh 65.5/66.3/65.7, STL 68.2/68.0/70.1 kg/s/m,
  mean +4.4 % (limit 5 %). The offset is systematic (every STL run above every mesh run);
  not traced (candidates: cut floor cells beside the gap, the gap below the floor line,
  corner cells). No mass lost while closed (4e-14). 60 deg hopper (half angle 30 deg,
  walls from the gap edges to the side walls): mass flow, empties 5.0 s after opening,
  rate 66.2 kg/s/m = 1.01 times the flat tank; Rose-Tanaka (tan 30)^-0.35 = 1.21. No pass
  band (brief). Tutorial `tutorials/stlHopper` (0.1 MPa default, 2 min). Outlet is still
  a mesh patch. Next: WP 4.3b (moving bodies).
