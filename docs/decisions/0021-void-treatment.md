# 0021 Void and free surface in the Lagrange-plus-remap step

Date: 2026-11-14. Status: accepted for the prototype (WP 3.2p); the C++ WP 3.2 follows it.

## Context

Target T0 needs one material plus void: a free surface that moves, stays sharp and
carries no stress, material that can separate (cohesionless, tension cut-off 0), and
empty space the material can fall into (the outlet jet in 3.3p). The prototype is the
WP 3.1p scheme (`prototype/silo2d.py`, decision 0020): symmetric dimensional splitting,
acoustic Riemann solvers, gravity by subtracting the discrete weight G. The brief named
three closures for partly filled cells: (a) volume change goes to the void first,
(b) homogenised cell (stress al s, shared strain), (c) an al threshold.

## Decision

Implemented in `prototype/void2d.py` (class `VoidBox`, sweep `void_sweep`).

1. **State.** `rho` stays the bulk density (cell mass / cell volume), so mass and every
   mass-weighted remap are unchanged; `al` is the material volume fraction. Velocity,
   stress and history are per unit material mass. al <= 1e-12: empty (no mass, v = 0,
   s = 0).
2. **Slab position, 1D SLIC per sweep.** A partly filled cell holds its material as one
   slab against its fuller neighbour; equal neighbours mean a layer parallel to the
   sweep (spread evenly). Walls mirror the cell next to them (treating walls as full
   broke the lake at rest: the end cells of a half-full row leaned onto the walls and
   left gaps between cells); free ends count as empty.
3. **Lagrange step: homogenised partly filled cells (closure b for forces).** Traction
   al s, impedance al z, mass al rho_m h: the wave speed is the material's for any al.
   Using the material's own stress and impedance instead blew up at al about 1e-3
   (its stress acted on a thousandth of the mass; amplification CFL^2 (G/M)/al per
   step).
4. **Faces.** A contact face has material touching it on both sides (not empty, slab
   not leaning away); it gets the Riemann solution. Every other face is free: zero
   traction, velocity from the one-sided relation of the side that touches it (or
   that holds material; zero between empty cells). For material without tensile
   strength (c = 0 or sigma_t = 0) a contact face whose solution would pull is made free
   too: this is what lets a block crack and open void inside the domain.
5. **Strain partition: closure (a).** Compression of a partly filled cell closes its
   void first; expansion first unloads the material elastically to zero normal stress,
   then opens void (for a material with tensile strength, a full cell stretches
   instead). Shear goes to the material. Closure (c) is not needed as a rule; only
   round-off dust (al < 1e-9 with mass) is made inert (s = 0, v = 0, al from its mass).
6. **Remap.** Material volume flux = the material part of the donor's swept region from
   the slab position (exact for a slab moving into void); mass flux = volume flux times
   the donor's Lagrangian material density (MUSCL as before where the donor and its
   neighbours are full). A donor left with less than 1e-9 of a cell gives all of it
   (no dust behind a moving edge); fluxes below 1e-9 of a cell are dropped (round-off
   velocities must not seed dust). Reconstruction is first order in and next to partly
   filled cells.
7. **Gravity in first-order cells.** The s_hat difference across a first-order cell is
   -w phi, w the cell's own weight, phi = 1 - |s_nn| / (w / 2) clipped to [0, 1]:
   phi = 0 keeps s_hat uniform (exact for a supported cell at rest, whose stress
   carries at least half its own weight); phi = 1 keeps the actual stress uniform
   (exact in free fall). With phi = 0 everywhere, a falling elastic block picked up
   stresses of 1 to 40 kPa and shed 0.1 percent of its mass sideways; with the rule it
   falls with stress 3e-9 Pa and acceleration g to 1e-13.

Superbee on al (the first option in the brief) was not tried: the geometric slab flux
was needed anyway for exact mass at moving edges and keeps a translating block at zero
partly filled cells per row after 50 cells, which no limiter does.

Measured (`tests/test_proto_void.sh`, 67 s): lake at rest with a half-full top cell
2.5e-16 m/s over 1000 steps; block translating 50 cells: mass 2e-14, edges exact (0.5
cell = the block itself), stress 3.5e-9 Pa; free fall acceleration g to 1e-13, landing
without loss (centroid stays above half the block less one cell); cohesionless block
pulled apart: void 2.000 cells per row for 2 v t = 2.000, no tension.

## Column collapse: not passed, split to WP 3.2q

T0 material, 10 cells across the short side, floor friction 35 deg, E 1 MPa:

| Case | Run-out (L_f - L_i)/L_i | Band (0019) | Top kept | Steepest slope per cell / two cells |
| --- | --- | --- | --- | --- |
| a = 0.5 | 0.25 (20 cells: 0.25, front unchanged from 0.33 s to 0.51 s) | 0.30 to 0.90 | yes | 50 / 47 deg |
| a = 2 | 1.2 | 1.2 to 3.6 | n/a | 60 / 53 deg |

Sensitivity (a = 0.5 unless noted): floor friction 20 deg 0.25 to 0.30; E 0.1 MPa 0.35
(a = 2: 1.4); Riemann jump terms halved 0.25 (a = 2, scaled by 0.3: 1.5); glass-bead
friction phi = delta = 22 deg 0.40 (a = 2: 1.7), that is 0.67 to 0.71 of the
experimental fit 1.2 a. The main deposit slopes are 25 to 33 deg; the steep part is the
front, a blunt nose 2 to 4 cells high.

Reading: the run-out is mesh-converged at a = 0.5 (10 and 20 cells agree) but depends on
E (lower sound speed, longer run-out), which a rate-independent model should not: the
acoustic Riemann solver's dissipation (about rho c h / 2, a numerical viscosity) slows
slow granular flow. Part of the shortfall may be real (phi = 35 deg against 22 deg glass
beads, decision 0019 already expects "at or below the fit"), but this cannot be
separated from the numerics yet. The blunt toe is the first-order free-surface cells.

So WP 3.2p closes with the void treatment and its four checks; the column collapse
checks (unchanged bands) move to WP 3.2q with the test `tests/manual/test_proto_collapse.sh`,
which joins `tests/` when it passes. Nothing about the bands was changed.

## Consequences

- **3.3p (discharge) is affected by the same dissipation.** Expect the rate to depend on
  E; 3.3p already plans an E study (its cost note). Measure it before trusting C.
- C++ WP 3.2 should follow rules 1 to 7. isoAdvector fits rule 6 in spirit (geometric
  flux of al with the face flux `phiL` the Lagrange step provides) and is unsplit, but
  it advects a passive al with a velocity; here the volume flux is tied to the
  Lagrangian cell lengths and slab positions, and the mass flux to the material density.
  Use isoAdvector's face-cut volume as the material volume flux, then rule 6 for mass;
  check against the translating-block and lake-at-rest checks before anything else.
- Rule 7's phi is a heuristic between two exact states; a block resting on another
  with little load sits between them. Its error is O(g h / c) in velocity in those
  cells only.
- The slab rule is 1D per sweep; a tilted surface is a staircase. Not tested at an angle.
