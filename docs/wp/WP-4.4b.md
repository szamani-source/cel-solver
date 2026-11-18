# WP 4.4b: Pile section moving through clay

State: done (2026-11-18). Needs WP 4.3b, 4.4. Split from WP 4.4 (2026-11-18).
Branch: the thread's `claude/project-thread-*` branch

## Goal

A curved moving STL body against a published limit: the lateral limiting pressure on a
deep pile in undrained clay, P = 9.14 cu D for a smooth pile and 11.94 cu D fully rough
(M. F. Randolph and G. T. Houlsby, "The limiting pressure on a circular pile loaded
laterally in cohesive soil", Geotechnique 34(4), 613-623, 1984). The reference was chosen
from published work as Sia asked (2026-11-18: "use the sources online"); its plane-strain
section fits the 2D solver, unlike the axial pile capacity.

## Read first

- `docs/decisions/0053`, `0055`; `tests/footing/footing_checks.py`

## Build

- `tests/test_pile.sh`, `tests/pile/pile_checks.py`; `sigma0` (uniform initial stress) in
  celVoidCheck.

## Done when

`tests/test_pile.sh` passes: smooth pile D = 0.1 m (10 cells, 64-sided prism) moving
sideways through weightless Tresca clay (cu 10 kPa, E/cu 300) confined at 10 cu:
P / (cu D) over a travel of 0.15 D to 0.3 D within 10 % of 9.14 and flat to 5 %; mass kept.

## Out of scope

- Rough pile (11.94: needs adhesion at the wall, decision 0054 item 4), axial pile
  capacity and installation (axisymmetry), sand.

## Log

- 2026-11-18: Done in one session (decision 0056). P / (cu D) = 9.29 (+1.7 %), flat from
  0.1 D. Learned: the plain void rules let cells beside a moving curved wall carry a
  trace of void and shake (force 1 to 5, blow-up near 0.3 D); `fullBlend yes` fixes it.
  Next: rough pile needs wall adhesion; axial pile needs axisymmetry.
