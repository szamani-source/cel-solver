# 0050 Column collapse with mu(I): the T0 rule set (WP 3.6e)

Date: 2026-11-18. Status: accepted (WP 3.6e done at 10 mm).

## Context

WP 3.6e: the granular column collapse with mu(I) against the experiments, run-out and
final height within the decision 0019 bands (run-out 0.5 to 1.5 times Lube's 1.2 a:
0.30 to 0.90 at a = 0.5, 1.2 to 3.6 at a = 2; top kept at a = 0.5; slopes not steeper
than phi) at a = 0.5 and 2, at T0's real E (10 MPa). Same set-up as
`tests/test_collapse.sh` (10 mm, smooth left wall, floor 35 deg, from K0 rest, shear cap
and rule 9), with the decision 0037 placeholders.

## Findings

1. **Bug: material stuck against a ceiling.** At 10 MPa the release throws a little
   material up to the domain's top wall. A wall mirrors the cell for the slab side
   (rule 2 of decision 0021), so a partly filled cell under a wall with void below
   leaned up against it, never touched its lower face and fell in place forever (16 m/s
   at 2 s; the kinetic energy at the end was the run's peak, so "stops" was never
   reached). Fix: for the slab side only, a wall or prescribed-velocity face whose
   outward normal points against gravity counts as empty. The smooth flag keeps the
   mirror, so a full cell under a lid stays second order (inclined layer with a lid
   unchanged). Applies to every free-surface run; the CI tests pass unchanged.
2. **mu(I) + flowImpedance alone** (the tank's rule set, decision 0045): 1 MPa run-outs
   0.35 / 1.50 in band, but 10 MPa 0.25 (below 0.30) / 1.20 (edge).
3. **Adding `normalCap yes`** (rule 12, decision 0035): the normal Riemann jump is the
   remaining E-dependent resistance in the collapse (decision 0032 said so for plain MC).

| Rule set, 10 mm | a = 0.5, 1 MPa | a = 2, 1 MPa | a = 0.5, 10 MPa | a = 2, 10 MPa |
| --- | --- | --- | --- | --- |
| Mohr-Coulomb (CI test) | 0.35 | 1.40 | 0.20 | (splash, 4.3 before the fix) |
| mu(I) + flowImpedance | 0.35 | 1.50 | 0.25 | 1.20 |
| mu(I) + flowImpedance + normalCap | 0.40 | 1.70 | 0.30 | 1.30 |

   Top kept to 0.1 mm at a = 0.5, slopes 24 to 26 deg, mass 3e-15. The 10 MPa a = 0.5
   column still creeps at 2 s (kinetic energy 1e-3 of peak), so its run-out is a lower
   bound.
4. **The same rule set on the other two problems:** inclined layer passes (profile
   error 3.5 / 0.8 / 0.8 % at 0.1 / 1 / 10 MPa); tank bands 73.6 / 61.0 / 58.8 kg/s/m
   (C 0.89 at 0.1 MPa, was 0.79; 1 and 10 MPa 17 and 20 % slow, inside the 20 %
   uncertainty accepted for now; 10 MPa band 51.6 to 67.6).

## Decision

- WP 3.6e is done at 10 mm: `tests/manual/test_collapse_mui.sh` checks all four cases
  against the bands and passes (about 2.5 min on 4 cores).
- **The T0 rule set is mu(I) + `flowImpedance yes` + `normalCap yes`** (with shear cap
  and rule 9, as before): the one setting found that keeps the layer exact, puts the
  collapse in band at T0's E and gives the tank its highest C. Defaults in the code are
  unchanged (all switches off); the tutorial is not changed in this package.
- The run-outs sit at the bottom of the bands. Expected in part: T0's mu_s = 0.70 is
  rougher than the glass beads and sand of the experiments (mu_s about 0.4 to 0.55),
  and the 10 mm band is mesh dependent (decision 0038). Not a sharp validation.

## Consequences

- `normalCap` cell-by-cell parallel agreement is still not exact (decision 0035).
- `tests/collapse/collapse_checks.py` reads `CEL_COLLAPSE_CHECK_ALL=1` to check every E.
- Next: switch `tutorials/siloDischarge` to the T0 rule set (small package), then WP 3.5
  (cohesion).
