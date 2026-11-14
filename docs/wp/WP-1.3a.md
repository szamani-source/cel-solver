# WP 1.3a: 2D plane-strain prototype, elastic waves

State: not started
Branch: `wp-1.3a-2d-elastic`

WP 1.3 is split in two (recorded in `docs/work-packages.md`): 1.3a builds the 2D scheme
for elastic material; 1.3b adds rotation, plasticity and writes the formulation decision.

## Goal

`prototype/lr2d.py`: the 1D cell-centred Godunov Lagrange step plus remap of WP 1.1,
extended to 2D plane strain on a fixed uniform Cartesian mesh, for linear elastic
material. It shows whether the 1D choices carry over to more dimensions, where shear
waves and transverse stress appear.

## Read first

- `docs/decisions/0005-cell-centred-lagrange-step.md`, `0008-plastic-wave-numerics.md`
- `prototype/lr1d.py` (class `Bar1D`), `prototype/lr1d_plastic.py`
- `docs/wp/WP-1.1.md` and `docs/wp/WP-1.2.md` logs

## Build

Suggested scheme (change it with a decision note if it fails; record what was tried):

- Fields per cell: density, velocity (vx, vy), stress (xx, yy, zz, xy), all cell-centred.
- Direction-split cycle (Strang or alternating x-y / y-x): each sweep is the 1D Lagrange
  step followed by the 1D remap, so cell volume change is exact per sweep.
- Face Riemann solver in the sweep direction: normal velocity and normal stress use the
  P-wave impedance `rho c_p`; tangential velocity and shear stress use the shear
  impedance `rho c_s`. Slopes limited in characteristic variables, MC limiter.
- Strain increment in a sweep from the face velocities (normal strain from the normal
  velocity difference, shear strain from the tangential velocity difference); stress
  update with the elastic matrix. No rotation correction yet (1.3b).
- Periodic and fixed-wall boundaries. Vectorised numpy; a 200 x 200 case should run in
  well under a minute.
- `tests/test_proto_2d_elastic.sh`

## Done when

`tests/test_proto_2d_elastic.sh` passes:

- Plane P-wave along x, along y, and at 45 degrees (periodic box): speed within 1 percent
  of `sqrt((K + 4G/3)/rho)`; smooth pulse converges at order 1.7 or better over 3 meshes.
- Plane S-wave along x and at 45 degrees: speed within 1 percent of `sqrt(G/rho)`; same
  order criterion.
- A wave along x on the 2D mesh matches `lr1d.py` on the same resolution to 1e-10
  (the 2D code reduces to the 1D code).
- x-y symmetry: a problem mirrored about the diagonal gives mirrored fields to 1e-12.
- Step P-wave with the material streaming at 0.25 c diagonally: no new extrema above
  1 percent of the jump; mass and both momentum components conserved to round-off.
- Existing tests still pass unchanged.

Record in the log (not pass criteria): the 45-degree wave errors compared with the
axis-aligned ones, and the run time per cell per step.

## Out of scope

- Rotation and objective stress rate, plasticity, the formulation decision (all 1.3b).
- More than one material, void, contact.

## Log

(none yet)
