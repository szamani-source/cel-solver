# WP 1.1: 1D elastic Lagrange-plus-remap prototype

State: done (waiting for review)
Branch: `claude/project-thread-as2tdz` (the session's assigned branch; the planned name was `wp-1.1-1d-elastic-remap`)

## Goal

A small Python code (`prototype/lr1d.py`) that propagates an elastic wave in a 1D bar on a
fixed mesh using an explicit Lagrange step followed by a second-order monotone remap.

## Read first

- `CLAUDE.md`, `docs/decisions/0001-architecture.md`
- Benson (1992), sections on the Lagrange step and advection (van Leer MUSCL)

## Build

- `prototype/lr1d.py`: all fields at cell centres, MUSCL-Hancock Godunov Lagrange step
  with an acoustic Riemann solver, MUSCL remap of mass, momentum and stress.
  (The brief first asked for a staggered grid with a central-difference step; decision
  0005 explains why that was replaced.)
- `tests/test_proto_1d_elastic.sh`

## Done when

- Smooth Gaussian pulse: L2 error converges at order 1.8 or better over 3 mesh levels.
- Step wave: no new extrema (monotone), wave speed within 1 percent of sqrt(M/rho).
- Mass and momentum conserved to round-off.

## Out of scope

- Plasticity (WP 1.2), more than one material, 2D.

## Log

### 2026-11-14

Done:
- Wrote the checks first (`lr1d.py --verify`, run by `tests/test_proto_1d_elastic.sh`).
  All three cases run on a periodic 1 m bar, rho 2000 kg/m^3, M 50 MPa, with the material
  streaming at 0.25 c so the remap moves real amounts. Test runs in about 1 s.
- Built the staggered version from the brief. Smooth order and conservation passed, but
  no artificial viscosity made the step wave monotone without losing second order
  (table in decision 0005). Kept in git history.
- Replaced it with a cell-centred MUSCL-Hancock step plus Riemann solver. All checks pass:
  orders 1.88 and 1.83 (200/400/800 cells), largest new extremum 1.8e-6 Pa on a 10 Pa
  jump, front speed error 0.04 percent, mass and momentum change below 1e-15.

Learned:
- Characteristic slopes must be built from differences with the cell's own impedance.
  Building them from s - z v mixed density noise times the streaming velocity into the
  stress and dropped the order to 1.2.
- A square pulse set up as sig = f, v = -f/(rho c) is not a discrete right-going wave on
  a staggered grid. It launches a spurious left-going wave of about 1 percent.
- The "new extrema" floor is physics plus round-off: the hypoelastic bar is nonlinear at
  O(strain), and the strain update rounds at about M * 1e-16 per step. The test uses a
  10 Pa pulse and a tolerance of 1e-6 of the jump. Its comments give the reasons.
- Smooth order is limited by peak clipping. It is 1.78 at 0.5 c streaming
  (`--v0 0.5` fails the 1.8 bar). An extremum-preserving limiter would fix it if needed.

Update: Sia approved the switch to cell-centred fields with a Riemann-solver step
(decision 0005 accepted). Decision renumbered from 0004 because WP 1.4 took 0004.

Next: Sia merges PR #3. WP 1.2 (plasticity and history remap) builds
on this file.
