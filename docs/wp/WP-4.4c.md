# WP 4.4c: Strip footing on sand (N_gamma)

State: closed without its done-when (2026-11-18), decision 0057. Needs WP 4.3b, 4.4.
Branch: the thread's `claude/project-thread-*` branch

## Goal

The strip footing of WP 4.4 on frictional, heavy sand against a published exact value:
the self-weight bearing factor N_gamma of a smooth strip footing, 7.653 for phi 30
(C. M. Martin, "Exact bearing capacity calculations using the method of characteristics",
Proc. IACMAG, Turin, 2005; tabulated in the Optum GX verification "Bearing capacity").
Sand values from published typical ranges (loose sand: phi 30, gamma 17.7 kN/m3,
E 10 MPa, nu 0.3), as Sia asked ("use the sources online").

## Done when

`tests/manual/test_sand_footing.sh` (phi 30, psi 30, 10 mm cells): the peak of
q / (0.5 gamma B N_gamma + gamma s N_q) within 15 % of 1, mass kept to 1e-10.

## Out of scope

- Rough footings, the N_gamma of other angles, embedded footings.

## Log

- 2026-11-18: Closed without its done-when (decision 0057). Associated flow (psi 30) blows
  up at a settlement of about 0.03 B; psi 0 reads N_gamma about 11 at 10 mm and 9.7 at
  5 mm (exact 7.65; a non-associated answer should be lower, not higher). Kept as a
  manual test, not weakened. Next: find why the reading is high (mesh, the surface
  rules near the footing edge) before a rough footing.
