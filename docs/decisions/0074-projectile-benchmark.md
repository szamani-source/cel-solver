# 0074 Projectile benchmark: free slug into clay, steel against rigid (WP 5.4)

Date: 2026-11-18. Status: accepted.

## Context

WP 5.4 (target problem T3, projectile into soil, deformable FE body) asks for a match
with a reference run. Published projectile results (Forrestal and Luk; Young's
equation) are for ogive noses at hundreds of m/s, axisymmetric; the axisymmetric
wedge is not ready (WP 3.8a/b) and the solver has no rate effects. The plane-strain
machinery that has a published check is the smooth strip footing in Tresca clay
(WP 4.4, Prandtl 2 + pi, -8 % at 10 mm).

## Decision

1. **Case.** A flat-nosed slug (0.095 x 0.1 m, plane strain, tungsten mass 1.83 kg in
   the 10 mm slab) falls 7.5 mm then strikes a 0.15 m layer of Tresca clay (cu 10 kPa,
   E 3 MPa, nu 0.45 for undrained, gravity) at 1.07 m/s, smooth, free (no drive).
2. **Reference run.** The same slug as a free rigid body (decision 0058). The steel
   hex8 slug (decisions 0070, 0072, 0073) must stop within 5 % of its penetration.
3. **Sanity band.** Both against an energy estimate: impact energy plus the work of
   gravity equals the work of Prandtl's resistance with Skempton's depth factor
   (1 + 0.2 d / B), 32.3 mm. Band 30 %.
4. A free deformable body starts with its `velocity` at t = 0 (was at rest).

## Consequences

- `tests/test_projectile.sh` (about 2.5 min): steel 23.75 mm, rigid 23.84 mm (-0.3 %);
  both 26 % short of the estimate; at rest; density 1.021 rho0 (band 1.03: elastic
  compression under the impact, bulk modulus 10 MPa; with nu 0.3 it reads 1.11).
- Known issue: the slug penetrates less than the quasi-static estimate, more so
  slowly (-33 % at 0.5 m/s, -26 % at 1 m/s, -21 % at 2 m/s), so not inertia of the
  clay. The static footing (WP 4.4) reads 8 % low. Not isolated; left for the accuracy
  work after the features (Sia, 2026-11-18).
