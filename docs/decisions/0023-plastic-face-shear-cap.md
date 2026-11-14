# 0023 Face shear traction capped at yield (column collapse run-out)

Date: 2026-11-14. Status: accepted for the prototype (WP 3.2q); the C++ Lagrange step
(decision 0018) should follow it before the silo runs in celFoam.

## Context

WP 3.2p (decision 0021) left the column collapse short: a = 0.5 ran out 0.25 (band 0.30
to 0.90), and the run-out grew as E fell. A rate-independent model should not care
about E. Measurements in WP 3.2q (`prototype/void2d.py`, a = 0.5, 10 cells across H,
E in Pa, run-out (L_f - L_i)/L_i, resolution 0.05):

| Change from WP 3.2p | E 1e4 | 3e4 | 1e5 | 1e6 | 3e6 | 1e7 |
| --- | --- | --- | --- | --- | --- | --- |
| none | 0.50 | 0.40 | 0.35 | 0.25 | 0.25 | 0.20 |
| this decision | | | 0.45 | 0.40 | | 0.25 |

What it is not:

- Reconstruction order next to the surface: second order in every full cell (also next
  to partly filled cells) gave the same 0.25 and 0.35.
- The time step: CFL 0.5, 0.25, 0.125 at E = 0.1 MPa gave 0.35, 0.40, 0.40 (a per-step
  error would go the other way).
- Mesh: 10 and 20 cells agree (0.25 before, 0.40 after). A numerical viscosity of the
  form rho c h dv/dx would fall with h; this one does not, because the shear band is one
  cell wide at any mesh, so the face velocity jump dv stays finite and the viscous stress
  z dv / 2 is mesh independent and grows with z = sqrt(rho G), that is with sqrt(E).

What it is: the jump terms of the acoustic Riemann solver. Scaling them by 0.3 at E = 1
MPa raised the run-out from 0.25 to 0.40 (0.1: 1.2, but the column fell apart). The
jump term z_s (vR - vL) / 2 adds a face shear stress that a plastic material cannot
carry: in a band at yield the face traction can exceed the Mohr-Coulomb strength by the
order of rho g H.

## Decision

Rule 8 of the void prototype (`shear_cap` in `prototype/void2d.py`, on by default,
`SHEAR_CAP`): at each interior face the shear jump term is scaled back so that the face
stress stays inside the yield surface.

- Face stress for the check: the Riemann normal traction (with its jump term) and shear
  traction, with the mean of the two cells' tt and zz stresses. Mohr-Coulomb yield
  function on its principal stresses (c and phi from the material parameters).
- The largest fraction beta in [0, 1] of the shear jump with f no larger than without
  the jump (or zero) is found by bisection (20 steps, per face, vectorised).
- A capped face slips: the cell on each side takes the tangential velocity of its own
  characteristic at the capped traction, v_L + (t - s_L) / z_L and
  v_R - (t - s_R) / z_R, exactly as a frictional wall does (decision 0020). The shear
  strain of each cell uses the velocity it sees. Momentum is exact (one traction per
  face); the slip is plastic work at the face.
- The normal jump term is not capped.

Tried and rejected:

- A Coulomb limit on the face plane only (|t| <= tan(phi) |s_n|): no change (0.25). The
  shear bands are inclined, so the face planes are not at their limit when the band is.
- Capping the normal jump term too, with one normal face velocity (all faces, or only
  extending ones), and capping the whole shear traction rather than its jump term: the
  column fell apart or the top sank 1.5 cells at E = 1 MPa. The normal jump term is the
  only damping of odd-even velocity modes in a perfectly plastic cell.

## Result (`tests/manual/test_proto_collapse.sh`, E = 1 MPa, 10 cells)

| Check | Before | After | Band |
| --- | --- | --- | --- |
| a = 0.5 run-out | 0.25 | 0.40 (20 cells: 0.40) | 0.30 to 0.90 |
| a = 0.5 top kept | yes | 0.06 cell lower | within 1 cell |
| a = 2 run-out | 1.2 | 1.6 | 1.2 to 3.6 |
| Steepest deposit slope per cell (two cells) | 50 (47), 60 (53) deg | 45 (35), 59 (53) deg | at most 37 deg |
| Mass | 1e-14 | 5e-14 | 1e-10 |

`tests/test_proto_void.sh` (WP 3.2p) passes unchanged.

## Consequences

- Run-out passes both bands and is mesh converged. It still depends on E (0.45, 0.40,
  0.25 at 0.1, 1, 10 MPa) through the uncapped normal jump term, so stiffer material
  under-predicts flow. Measure the 3.3p discharge rate against E before trusting C; use
  E = 1 MPa as the collapse does if T0's 10 MPa is slower.
- The slope check still fails, so the collapse test stays in `tests/manual/`. The steep
  part is the front and the 1-cell steps of the surface: the 1D slab rule of decision
  0021 stands a partly filled top cell against its fuller sideways neighbour in the x
  sweep, so a thin top layer cannot move sideways until the cell is full, and column
  heights snap to whole cells (45 deg per step). A 2D-normal slab rule (spread when the
  volume fraction changes more across the sweep than along it) gave a smoother surface
  without the cap (0.35, slopes 45 deg) but blew up with it. Split to WP 3.2r.
- The C++ Lagrange step (decision 0018) has the same acoustic face solver and will show
  the same over-damping of slow plastic flow; it needs this cap (the yield function is
  available from the constitutive kernel).
