# 0091: Sand footing, hopper and dilation angle: what the accuracy pass found

Date: 2026-11-19 (accuracy work before packaging, thread "Fix accuracy issues before
packaging"; items 3, 4 and 1 of the ranking Sia replied to)

## Context

Sia's reply to the ranked accuracy issues (2026-11-19): mass flow and mesh dependence
are secondary until a validation case is agreed; discharge should depend on a small
dilation angle (0.5 to 1 deg); cells must resolve the shear band (10 to 20 grain
diameters) and mesh quality matters; find where errors come from. Decision 0090 covers
the wedge item.

## Findings

### Sand footing (decision 0057), 10 mm cells, psi 0

| Run | Implied N_gamma (s 0.01 to 0.06 B) |
| --- | --- |
| smooth base (as 0057) | 11.3 to 10.2 |
| rough base (friction 30 deg) | 14.1 to 13.2 |
| smooth, edges and base on cell faces (x0 0.25, gap 0.01) | 12.0 to 10.6 |

Exact associated values (Martin 2005): smooth 7.653, rough 14.754. A non-associated
(psi 0) answer lies below these; with the Davis reduction (tan phi* = sin phi cos psi /
(1 - sin phi sin psi), phi* 26.6 deg) Martin's table gives smooth about 4.4. So the
smooth footing reads high by more than the 27 to 44 % recorded, and the smooth and rough
readings are closer (1.34 times) than the exact ratio (1.93): part of the excess is not
roughness.

**A body just above the free surface disturbs the material at rest.** A footing held
2.5 mm above a sand layer, not moving (`tests/manual/test_body_above_surface.sh`):

- psi 0: the top cells at the footing's edges move at 0.1 to 0.2 m/s by 0.2 s, with
  stresses near 20 kPa where the geostatic stress is about 90 Pa; the top row under the
  footing becomes partly filled (alpha 0.8) as its small cut cells are mixed with the
  material below.
- psi 30: the top cell at the footing's right corner (below a cut cell 81 % open)
  gains stress exponentially (5 kPa at 0.03 s, 500 kPa at 0.09 s) and the run blows up
  near 0.12 s, before any contact. This is the "associated flow blows up" of 0057: it
  is not the footing push.
- The same layer without the footing stays at rest to 1e-17 m/s, for psi 0 to 30. With
  the footing's base on a cell face (no small cut cells at the start) the run still
  blows up once the footing enters, and one cell packs to 2.3 rho0.

These cells confine the surface beside the footing exactly where the failure wedges come
out, which raises a frictional bearing reading and not a Tresca one (clay footing within
2 to 8 %). That fits the sand footing being high and the clay footing being right.

### Hopper (decision 0052), 5 mm cells

| Cells | Flat STL tank | 60 deg hopper | Ratio |
| --- | --- | --- | --- |
| 10 mm | 68.2 (maxCo 0.40) | 66.2 | 0.97 (1.01 against the mesh-wall mean) |
| 5 mm | 81.3 | 87.0 | 1.07 |

Rose-Tanaka predicts 1.21 for 3D orifices. The ratio moves toward it with refinement;
at 5 mm the outlet is 12 cells and a 10 to 20 mm band 2 to 4 cells. 2.5 mm (about 1.5 h
per run) not run.

### Dilation angle in the T0 tank, 10 mm, three runs each (kg/s/m)

| psi | E 0.1 MPa | E 1 MPa | 1 MPa against 0.1 MPa |
| --- | --- | --- | --- |
| 0 (decision 0046) | 65.8 | 52.0 | -21 % |
| 0.5 deg | 63.0 | 55.7 | -12 % |
| 1 deg | 61.8 | 52.5 | -15 % |

A small dilation angle lowers the rate at 0.1 MPa by 4 to 6 % and narrows the stiffness
gap; band scatter is 1 to 2 %. psi was set through `CEL_EXTRA_MATERIAL` (the later
entry overrides the driver's psi 0).

## Decision

1. Record the hovering-body test as a failing manual test; fixing it is the next
   package (body near a free surface), done when `tests/manual/test_body_above_surface.sh`
   passes. The sand footing is then re-run, with psi 30 against Martin's exact values.
2. The T0 tutorials keep psi 0 for now; psi 0.5 to 1 deg is recommended for bulk solids
   and will be the default once Sia's validation case fixes the material.
3. Hopper and STL rates at 10 mm are not resolved; use 5 mm or finer for silo results
   (cells no larger than half the shear band, 5 mm for 1 mm grains).

## Consequences

- The sand footing and the associated-flow blow-up share one cause that can be tested
  in 45 s, instead of a 1.5 h run at 2.5 mm.
- The projectile (decision 0074) also starts above the surface; it may share the cause.
