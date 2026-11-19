# 0100: Conical hopper validation case: what we compare and how

Date: 2026-11-19 (WP H.1)

## Context

Sia supplied two papers for the validation case they said they would set with us
(accuracy pass, 2026-11-19): Zheng and Yu (2015), Abaqus Eulerian runs of a conical
hopper, and Couto et al. (2013), the wheat silo experiment that paper compares with.
The runs use our axisymmetric wedge with Mohr-Coulomb and Drucker-Prager. Several
values are not in the papers or differ between them.

## Decision

- Three cases (`docs/spec/conical-hopper.md`): A, the Couto silo (wall pressure at rest
  and at 0.5 s of free discharge, against the experiment and P1's FEM); B, P1's basic
  hopper (outlet velocity plateau, mass flow); C, P1's flat-bottomed silo (Beverloo
  C = 0.59).
- Couto's discharge was throttled by a screw conveyor (0.6 to 4.9 kg/s). We compare
  the discharge pressures with free discharge at 0.5 s, as P1 did, because P2 found
  the peak almost independent of the rate.
- Dilation 1 deg (Sia: 0.5 to 1) is the main setting; P1's 5 deg is a sensitivity run.
- Rate-independent models (no mu(I)), as in P1. Drucker-Prager fitted to triaxial
  compression.
- Case B's cylinder height (1.701 m) follows from H = 2.48 m, alpha 30 deg and D0 0.1 m;
  the fill is to the top. P1 does not state either.
- Wheat grain about 4 mm (typical published kernel sizes), so the shear band is 40 to
  80 mm and the cells 10 mm (5 mm as the check).
- Pass bands are set before any run (spec page). The upper cylinder cells (GA4 to GA7)
  read 11 to 31 % above Janssen at rest while the lower ones match it; there either the
  experiment or Janssen counts.
- The reference numbers live in `docs/spec/hopper-data/*.csv`: Table 1 of P2 copied by
  hand, P1's figures digitised by eye (the precision is written in each file).
  `tests/test_hopper_spec.sh` checks them against the geometry and P1's closed forms.

## Consequences

- The conical hopper needs walls on the wedge that are not mesh lines: WP 3.8d (STL
  walls on a wedge) becomes H.2 and comes first.
- Decisions for this series are numbered from 0100.
- The papers are not in the repository (copyright); they are on Sia's laptop.
