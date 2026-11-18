# 0057 Sand footing (N_gamma) does not yet match the exact value (WP 4.4c)

Date: 2026-11-18. Status: accepted (finding; package closed without its done-when).

## Context

WP 4.4 and 4.4b checked moving bodies in Tresca clay (Prandtl 5.14 within 2 % at 5 mm,
Randolph and Houlsby 9.14 within 2 %). WP 4.4c adds weight and friction: a smooth strip
footing B = 0.1 m pushed at 0.01 m/s into cohesionless Mohr-Coulomb sand (phi 30, c 0,
rho 1800, E 10 MPa, nu 0.3, K0 1 - sin phi, `fullBlend yes`, frictionless mesh walls,
0.6 x 0.3 m layer). Reference: Martin (2005) exact smooth N_gamma 7.653 (rough 14.754),
plus the surcharge of the settlement, q_ref = 0.5 gamma B N_gamma + gamma s N_q
(N_q 18.40).

## Findings

| Run | Cells | Implied N_gamma | Against 7.653 |
|---|---|---|---|
| psi 30 (associated) | 10 mm | about 21 at s 0.03 B, then blows up | n/a |
| psi 0 | 10 mm | 10.8 to 11.3 | +44 % |
| psi 0 | 5 mm (678 s) | 9.6 to 9.85 | +27 % |
| psi 0, shearCap and normalCap | 10 mm | 9.1 to 9.7 | +19 to +27 % |

A non-associated answer should lie below the associated exact value, so the reading is
high by more than the table shows. It falls with the mesh (11 to 9.7 from 10 to 5 mm),
so part is resolution (B is 10 cells; the failure wedge under a smooth footing is
shallow, about 0.5 B), and part may be the surface cells beside the footing, where the
free-surface rules (decision 0051) carry the heave.

## Decision

1. Close WP 4.4c without its done-when, as WP 3.2t was; the check stays at 15 % and
   associated flow, in `tests/manual/test_sand_footing.sh` (options `--fine`, `--psi`).
2. Not in CI (a failing test).
3. Associated flow with `fullBlend` is unstable here; use psi well below phi for now.

## Consequences

- Moving bodies are checked in cohesive (Tresca) material only. Bearing results in sand
  are about 30 to 45 % high at 5 to 10 mm cells; treat them as upper estimates.
- Next steps if a sand footing is needed: run 2.5 mm (B 40 cells, about 1.5 h), check
  the surface cells beside the footing, and the associated blow-up.
