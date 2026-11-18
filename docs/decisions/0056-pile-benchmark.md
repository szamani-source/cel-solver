# 0056 Pile benchmark: Randolph and Houlsby, confined clay, blended rules (WP 4.4b)

Date: 2026-11-18. Status: accepted.

## Context

WP 4.4 named a pile. Sia asked (2026-11-18) for references from published sources. The
solver is 2D, so an axial pile (axisymmetric) is out of reach; the plane section of a
deep pile moving sideways is plane strain and has a classical exact limit: Randolph and
Houlsby (1984, Geotechnique 34(4), 613-623), P / (cu D) = 9.14 smooth, 11.94 fully rough,
for full flow of Tresca clay around the pile.

## Decision

1. **Case.** Weightless Tresca clay (phi 0, cu 10 kPa, E 3 MPa, nu 0.3) fills a 0.6 x
   0.5 m box (6 D by 5 D), frictionless mesh walls. Smooth pile D = 0.1 m, a 64-sided prism
   of the circle's area, moved along x at 0.025 m/s (ramped over 0.05 s) for 0.3 D.
2. **Full flow.** A uniform confining pressure of 10 cu (`sigma0` in celVoidCheck) keeps
   the clay against the back of the pile, as the reference assumes; it gives no net
   force on the pile. Without it the clay leaves a gap behind the pile (walls carry no
   tension, decision 0054 item 4).
3. **Rules.** `fullBlend yes` (decision 0049). With the plain void rules the cells next
   to the moving curved wall carry a trace of void (alpha about 1 - 1e-6), take the
   partly filled rules, compress without stress and shake: the force stays at 1 to 5 cu D
   and the run blows up after a travel of 0.3 D. The blend treats them as full; the force
   is then steady.
4. **Measure.** P = force of the pile on the clay per depth; P / (cu D) averaged over a
   travel of 0.15 D to 0.3 D, within 10 % of 9.14, flat to 5 %.

## Results

10 mm cells (D = 10 cells): P / (cu D) = 9.29 (+1.7 %), flat from a travel of 0.1 D;
largest speed in the clay 0.05 m/s (twice the pile's); mass to 1e-14, nothing lost.

## Consequences

- `tests/test_pile.sh` is the done-when (CI, about 3.5 min).
- Moving bodies in materials without mu(I) (no nearly-full snap, decision 0044) should
  use `fullBlend yes`. The footing (decision 0055) passed without it because its clay
  never carried void next to the footing.
- Not done: rough pile (needs wall adhesion), finer mesh, the axial pile.
