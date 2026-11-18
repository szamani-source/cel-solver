# 0044 Stiff inclined layer: wall ghost and nearly full cells (WP 3.6c2)

Date: 2026-11-17. Status: accepted (both behind `flowImpedance yes`, default off).

## Context

After WP 3.6c (decision 0040) the mu(I) inclined layer (31 deg, 5 mm cells, free
surface) was within 7.8 / 22 / 136 % of the Bagnold profile at E = 0.1 / 1 / 10 MPa.
Two faults were left: at 10 MPa the floor cell barely moved and the body was 25 to 35 %
slow; at 1 and 10 MPa the top cell slid over the body.

Diagnosis (traces of the face velocities and tractions, and of the top cell step by step):

1. **Floor cell**: walls are zeroGradient for U, so the reconstruction mirrored a
   flowing cell's own velocity and its face state at the wall was the whole cell
   velocity u0. The wall shear jump (eta / distance to the wall, decision 0040 item 2)
   then acted on u0 on top of the cell's mu(I) stress, which already carries the rate:
   the wall drag was counted twice. At 10 MPa the floor cell stick-slipped near zero
   speed at the static friction angle.
2. **Top cell**: its stress is tens of Pa, an elastic strain of about 4e-6 at 10 MPa.
   Wave and round-off noise open void of 1e-9 to 1e-5 in it. Any void above 1e-12
   switched the cell to the partly filled rules (first order, no flow viscosity, a
   leaning slab, shear cap at the faces, decision 0021) and back. Each switch kicked its
   stress (the shear changed sign within one step), and then it sat stress-free and slid.

## Decision

With `flowImpedance yes` only:

1. **Wall ghost of a flowing cell**: the velocity reconstruction at a wall face uses the
   wall's last face velocity (vStar, 0 where it sticks, the slip velocity where it
   slips) as the patch value, so the mirror ghost is 2 v_wall - u0 and the face state is
   the cell's extrapolated wall velocity. The jump then acts only on the discrepancy.
2. **Nearly full is full**: after the remap a cell with void fraction below 1e-4 is set
   to alpha = 1 (was 1e-12). Mass is untouched; the implied material density changes by
   at most 1e-4.

Tried and dropped:

- A support factor (rule 7) using |sigma & e| for tilted gravity: right in principle, no
  effect on the result once 2 was in.
- A one-sided tangential traction at a face between a flowing cell and a partly filled
  top cell: unstable (the top cell blew up), because the face lost the stress-difference
  coupling that the upwind velocity needs.

## Results (inclined layer, 4 s from the exact state, 5 mm, free surface)

| E (MPa) | WP 3.6c profile error | WP 3.6c2 profile error | u / u_exact by cell |
| --- | --- | --- | --- |
| 0.1 | 7.8 % | 3.6 % | 1.029 to 1.074 |
| 1 | 22 % | 0.9 % | 0.985 to 0.996 |
| 10 | 136 % | 0.8 % | 0.982 to 0.995 |

Smooth lid: 2.0 / 0.6 / 0.7 %. 2.5 mm at 0.1 MPa: 5.1 % (free surface), 3.4 % (lid). The
at-rest layer below the stopping angle: 1.2e-6 m/s (limit 1e-4). The stiff cases are now
the most accurate: the error that is left at 0.1 MPa grows slowly in time.

## Consequences

- WP 3.6c2 done-when met: `tests/manual/test_inclined_layer.sh` checks all three E at
  5 mm. CI (`tests/test_inclined_layer.sh`) now runs 0.1 and 10 MPa with the 5 % checks.
- At 0.1 MPa the layer creeps faster over time: +0.7 % of the surface speed per second
  at 5 mm, +1.05 % at 2.5 mm (just over the 1 % steadiness limit, recorded only); all
  cells run 3 to 11 % fast. Not found in this package; the soft case has the largest
  elastic strains, so the cause is likely in the elastic part of the predictor.
- The 1e-4 void tolerance also applies at a real free surface in mu(I) runs: a cell
  with less than 1e-4 void counts as full. Watch the silo surface in WP 3.6d.
- The wall ghost uses last step's wall velocity; the first step uses 0.
- Silo, collapse and outlet runs keep `flowImpedance` off: unchanged bit for bit.
