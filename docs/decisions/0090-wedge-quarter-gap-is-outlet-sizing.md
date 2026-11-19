# 0090: The wedge against quarter gap comes from how the square outlet was sized

Date: 2026-11-19 (accuracy work, item 7 of the ranking in the thread "Fix accuracy
issues before packaging")

## Context

Decision 0062 compared a cylindrical silo on the axisymmetric wedge with the same
tank as a 3D quarter model and found C = 0.678 (wedge) against 0.635 (quarter), 6.8 %
apart. The quarter model's outlet is a square (full side a = 0.06 m); its C was taken
from the circle of equal area, D_eq = 2 a / sqrt(pi) = 0.0677 m. Sia asked where the
6.8 % comes from.

## Finding

The gap is in the comparison, not in the solver. For a non-circular outlet the
Beverloo form uses the area and the hydraulic diameter D_h = 4 A / P:

    W = C rho (4/pi) A sqrt(g D_h)          (circle: C rho sqrt(g) D^2.5)

(Zhao, Sharma, Leonelli, Meiburg and Sauret, "Steady granular discharge across
aperture shapes", arXiv 2610.02543, eq. 15: Q = A rho phi sqrt(g D_h)). For the square, D_h = a, which is smaller
than D_eq, so the equal-area circle reads C low by about sqrt(a / D_eq) = 0.94 even
when the two flows are the same. With the empty annulus k d taken off each side:

    C_square = pi W / (4 rho sqrt(g) (a - k d)^2.5)

On the same runs (2026-11-19, 10 mm): wedge 0.678, quarter by D_h 0.680, gap +0.2 %.
The 2-cell outlet of decision 0062 still differs (0.69 against about 0.53 by D_h, recomputed from the rates in 0062), so the
3-cell rule of 0062 stands.

## Decision

`tests/wedgesilo/wedge_silo_checks.py` keeps its 10 % equal-area check and adds a
check by the hydraulic diameter, within 3 %.

## Consequences

- The wedge and the 3D quarter agree to 0.2 % at 10 mm: the wedge can be trusted for
  axisymmetric silos.
- Rates from square or slot outlets should be compared with Beverloo by D_h, not by
  an equal-area circle.
