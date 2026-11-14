# WP 1.1: 1D elastic Lagrange-plus-remap prototype

State: not started
Branch: `wp-1.1-1d-elastic-remap`

## Goal

A small Python code (`prototype/lr1d.py`) that propagates an elastic wave in a 1D bar on a
fixed mesh using an explicit Lagrange step followed by a second-order monotone remap.

## Read first

- `CLAUDE.md`, `docs/decisions/0001-architecture.md`
- Benson (1992), sections on the Lagrange step and advection (van Leer MUSCL)

## Build

- `prototype/lr1d.py`: staggered grid (velocity at nodes, stress and density at cells),
  central-difference Lagrange step, MUSCL remap of mass, momentum and stress.
- `tests/test_proto_1d_elastic.sh`

## Done when

- Smooth Gaussian pulse: L2 error converges at order 1.8 or better over 3 mesh levels.
- Step wave: no new extrema (monotone), wave speed within 1 percent of sqrt(M/rho).
- Mass and momentum conserved to round-off.

## Out of scope

- Plasticity (WP 1.2), more than one material, 2D.

## Log

(none yet)
