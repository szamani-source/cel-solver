# WP 4.3e: Rough walls (wall adhesion)

State: done (2026-11-18). Needs WP 4.3a, 4.4b.
Branch: the thread's `claude/project-thread-*` branch

## Goal

STL walls whose shear strength has an adhesion part, so a wall can be fully rough in
clay (decision 0054 left walls without adhesion). Checked on the rough pile of Randolph
and Houlsby (1984): P / (cu D) = 11.94.

## Build

- `adhesion <Pa>` per body in `constant/rigidBodies` (default 0).
- The wall Riemann solution caps the wall shear at adhesion + friction x pressure; the
  material still separates at zero normal stress.

## Done when

`tests/test_rough_pile.sh`: the pile section of WP 4.4b with adhesion cu, P / (cu D) over
a travel of 0.15 D to 0.3 D within 15 % of 11.94 at 10 mm cells, steady within 5 %, mass
kept. Within 10 % at 5 mm (`tests/manual/test_rough_pile_fine.sh`).

## Out of scope

- Adhesion on mesh-boundary walls and outlets; tension at walls (sticky walls).

## Log

- 2026-11-18: Done in one session (decision 0059). 10.47 at 10 mm (-12 %), 11.26 at 5 mm
  (-5.7 %): the rough limit converges more slowly than the smooth one (9.29 at 10 mm),
  as the shear layer at a rough wall is thin. Next: Taylor bar.
