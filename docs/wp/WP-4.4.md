# WP 4.4: Strip footing bearing capacity

State: done (2026-11-18). Needs WP 4.3b. Split (2026-11-18): the pile goes to WP 4.4b.
Branch: the thread's `claude/project-thread-*` branch

## Goal

The first geotechnical benchmark with a moving STL body: a rigid strip footing pushed
into undrained clay reaches Prandtl's bearing capacity, q = (2 + pi) cu.

## Read first

- `docs/decisions/0053` (moving bodies); `tests/stlwalls/stl_push_checks.py`;
  `docs/spec/validation-references.md` (Qiu, Henke, Grabe 2011 for later CEL cross-checks)

## Build

- `tests/test_footing.sh` and `tests/footing/footing_checks.py` on celVoidCheck
  (no solver change expected).

## Done when

`tests/test_footing.sh` passes: smooth strip footing B = 0.1 m (10 cells) on weightless
Tresca clay (phi 0, cu 10 kPa, E/cu 300), pushed at 0.025 m/s: q / cu over a settlement of
0.10 B to 0.18 B within 10 % of 5.14 and flat to 5 %; mass kept, nothing lost.

## Out of scope

- The pile (WP 4.4b): shaft friction and base capacity need the footing first.
- Frictional soil (N_gamma, N_q), rough footing, rate effects beyond a check of speed.

## Log

- 2026-11-18: Done in one session (decision 0055), no solver change. 10 mm: 4.73 (-8.1 %);
  5 mm: 5.04 (-1.9 %); 0.05 m/s adds 4 % (inertia), so the test pushes at 0.025 m/s.
  Force dips each time the footing bottom crosses a mesh line. Plot
  `WP-4.4-footing.png`. Next: WP 4.4b (pile), reference to agree with Sia.
