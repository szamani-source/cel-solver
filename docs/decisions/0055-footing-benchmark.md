# 0055 Strip footing benchmark: weightless Tresca clay, quasi-static push (WP 4.4)

Date: 2026-11-18. Status: accepted.

## Context

WP 4.4's done-when is "bearing capacity factor near 5.14": Prandtl's limit pressure of a
strip footing on undrained clay, q = (2 + pi) cu, for a weightless Tresca material and
a footing at the surface. It is the first check of a moving STL body (decision 0053)
against theory. The work package also named a pile; that is split off as WP 4.4b.

## Decision

1. **Case.** Weightless Mohr-Coulomb with phi 0 (Tresca), cu = 10 kPa, E = 3 MPa
   (E/cu 300), nu 0.3; 0.6 m wide and 0.3 m deep (6 B by 3 B), frictionless mesh walls,
   void above (`freeSurface yes`). Smooth rigid footing B = 0.1 m, an STL box whose sides
   cut cells at a quarter and whose bottom starts 2.5 mm below the surface (D/B 0.025).
2. **Quasi-static push.** Speed ramped from 0 to 0.025 m/s over 0.05 s (wave speed about
   50 m/s). At 0.05 m/s q is 4 % higher (inertia), at 0.025 m/s it is flat after 0.1 B.
3. **Measure.** q = footing force over B times depth; q / cu averaged over s = 0.10 B to
   0.18 B; pass within 10 % of 5.14 and flat to 5 % over that range.
4. **Mesh.** CI runs 10 mm (B = 10 cells, about 1 min); 5 mm is a manual check
   (`footing_checks.py WORK --fine`, 9 min).

## Results

| Cells | Speed | q / cu (0.10 to 0.18 B) | vs 5.14 |
| --- | --- | --- | --- |
| 10 mm | 0.025 m/s | 4.73 | -8.1 % |
| 10 mm | 0.05 m/s | 4.90 (0.08 to 0.15 B) | -4.8 % |
| 5 mm | 0.025 m/s | 5.04 (5.10 at 0.15 B) | -1.9 % |

The limit converges on Prandtl's value from below as the mesh is refined (B = 10 cells
-8 %, 20 cells -2 %). Plot: `docs/wp/WP-4.4-footing.png`.

## Consequences

- `tests/test_footing.sh` is the done-when (CI). Its 10 % band holds at 10 mm with a
  margin of under 2 %; a change that lowers the 10 mm limit by more will fail it, and
  the 5 mm run is then the arbiter.
- The force dips for a step or two (up to 25 %) each time the footing bottom crosses a
  mesh line (a row of small cut cells becomes covered and hands over). It does not move
  the average; a smoother hand-over would remove it.
- Not checked: rough footing, frictional soil (N_gamma, N_q), gravity with depth, the
  pile (WP 4.4b), comparison with Abaqus CEL (Qiu et al. 2011).
