# WP 3.5: cohesion in the silo

State: done at 10 mm (2026-11-18)
Branch: `claude/project-thread-50m2qy`

## Goal

The T0 tank with a cohesive material (c > 0): the material arches over a narrow
outlet and the flow stops; above a critical width it flows.

## Read first

Decisions 0015 (tension cut-off), 0021 (void rules 3 and 5), 0050 (T0 rule set), 0054.

## Done when

`tests/test_cohesion.sh` passes: c = 500 Pa stops at 80 and 100 mm and flows at 120 and
140 mm (Jenike 109 mm); c = 0 flows at 80 mm; mass to 1e-10, no stress in empty cells.

## Out of scope

Wall adhesion; cohesion that grows with consolidation (Jenike flow function); 5 mm
mesh check (next, if wanted).

## Log

- 2026-11-18: One session. The constitutive side already had c; the free-surface rules
  treated a cohesive material as unbreakable. Added one tensile strength (sigmaT, else
  c cot(phi)) at contact faces and in void opening (decision 0054). Bracket 100 / 120 mm
  at c = 500 Pa, same at E = 10 MPa; 40 / 60 mm at c = 250 Pa (Jenike 109 / 55 mm).
  Tutorial `tutorials/siloCohesion`. Left: stopped tanks creep at under 1 % of Beverloo
  (4.5 % at c = 250 Pa, 40 mm).
