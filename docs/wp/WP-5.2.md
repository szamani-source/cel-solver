# WP 5.2: Lagrangian-only Taylor bar

State: done (2026-11-18). Needs WP 5.1.
Branch: the thread's `claude/project-thread-*` branch

## Goal

The hex8 body (WP 5.1) in large plastic flow against a published impact result: the
Taylor bar of the Abaqus Benchmarks ("Impact of a copper rod", after Ferencz 1989), as a
3D quarter model.

## Read first

- `docs/decisions/0070`, `0071`; `tests/hex8/hex8Check/hex8Check.C` (taylorTest)

## Build

- `celHex8Check -taylor` (options `-nq`, `-nz`, `-hourglass`, `-hourglassType`) and
  `tests/test_taylor_bar.sh`.
- `hourglassType viscous` in hex8Body (decision 0071).

## Done when

`tests/test_taylor_bar.sh`: copper rod 32.4 mm long, radius 3.2 mm, 227 m/s onto a
rigid frictionless wall, E 110 GPa, nu 0.3, rho 8970, von Mises perfectly plastic
314 MPa, 80 microseconds, 6 x 6 x 45 quarter model: shortening within 5 % of 13.10 mm
and widening (radial displacement of the impact face edge) within 10 % of 5.528 mm (the
Abaqus C3D8R quarter model), volume within 1 %, no element inverted.

## Out of scope

- Hardening and rate laws (Johnson-Cook), friction at the wall, the Eulerian Taylor bar
  (WP 2.6, needs the axisymmetric wedge of WP 3.8).

## Log

- 2026-11-18: Done in one session (decision 0071). 13.14 mm (+0.3 %) and 5.76 mm
  (+4.2 %; the Abaqus axisymmetric run gives 6.006), largest plastic strain 4.2, about
  1.5 min. Learned: the stiffness hourglass control of WP 5.1 locks the impact face in
  large plastic flow (hourglass 0.05: widening 3.7 mm, 0.2: 2.4 mm; 0.005: 5.1 mm), the
  viscous form does not (0.05: 5.76, 0.1: 5.67). Next: WP 5.3.
