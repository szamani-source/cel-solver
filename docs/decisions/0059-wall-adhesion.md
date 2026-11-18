# 0059 Rough walls: adhesion in the wall shear limit (WP 4.3e)

Date: 2026-11-18. Status: accepted.

## Context

STL walls carried Coulomb friction only (decision 0047): the wall shear was capped at
mu times the wall pressure. A wall in clay is often rough, its shear limited by the
clay's strength (adhesion alpha cu) rather than by friction. Decision 0054 left wall
adhesion out. The pile section of WP 4.4b checked the smooth limit (9.14); the rough
limit of the same reference, 11.94 (Randolph and Houlsby 1984), needs adhesion.

## Decision

1. `adhesion <Pa>` per body (default 0). The shear limit of the wall Riemann solution is
   `adhesion - tn mu` (tn the normal traction, negative in compression), in both the
   predictor and the Lagrange step.
2. Separation is unchanged: the wall face carries no traction when the material pulls
   away (tn >= 0). A sticky wall (tension) is out of scope.
3. Mesh-boundary walls and outlets stay without adhesion.

## Consequences

- `tests/test_rough_pile.sh` (CI, about 3 min): rough pile, adhesion cu, 10 mm cells:
  P / (cu D) = 10.47 against 11.94 (-12 %, band 15 %). At 5 mm 11.26 (-5.7 %, band 10 %,
  manual, 25 min). The smooth pile reads 9.29 at 10 mm (+1.7 %): the rough limit needs
  the thin shear layer at the wall, so it converges more slowly with the mesh.
- Bit-for-bit unchanged without adhesion.
