# WP 1.3a: 2D plane-strain prototype, elastic waves

State: done (2026-11-14)
Branch: `claude/project-thread-npbl39` (thread branches replace the `wp-` names)

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

### 2026-11-14

Done:
- Wrote the checks first (`lr2d.py --verify`, run by `tests/test_proto_2d_elastic.sh`),
  then `prototype/lr2d.py` as the brief suggested: per sweep a P system, an S system and
  stationary s_tt, s_zz; MC characteristic slopes; acoustic Riemann solvers; 1D remap.
  Material rho 2000, M 50 MPa (as WP 1.1), nu 0.3. Smooth runs stream the material at
  0.25 c_p along the diagonal. All criteria pass; test takes about 5 minutes.
- Smooth Gaussian (width 0.1 m along the normal), 100/200/400 cells per metre:
  P 0 and 90 deg orders 1.85, 1.86; P 45 deg 1.84, 1.85; S 0 deg 1.78, 1.77;
  S 45 deg 1.87, 1.87. Speeds within 0.002 percent.
- Reduces to 1D bit-for-bit (difference 0.0). Mirror about the diagonal bit-for-bit.
- Step P-wave along x with diagonal streaming: no new extrema (0.000 percent); mass and
  momentum change below 4e-16. Fixed-wall box keeps its mass.

Recorded (not pass criteria):
- 45 deg vs axis-aligned L2 error at 400 cells: P 1.46e-4 vs 2.62e-4; S 1.43e-4 vs
  3.93e-4. The diagonal is more accurate (more cells across the pulse along the normal).
- Cost about 2.5 microseconds per cell per step on 400 x 400 (symmetric split, four
  sweeps per step), 3 to 5 on the thin 400 x 4 strips (numpy overhead). A 200 x 200 step
  case of 106 steps runs in about 10 s.
- Step P-wave at 45 degrees: largest new extremum 13 percent of the jump (s_xy, start-up),
  4.9 percent at the end of the run. Independent of CFL, limiter and streaming: the
  direction-by-direction Riemann solver. The brief's step criterion is read as the wave
  along x with diagonal streaming; the 45 degree number is printed only. Decision 0012.

Learned:
- Symmetric splitting (average of x-y and y-x) is needed for exact mirror symmetry.
- Matching 1D to 1e-10 needs identical arithmetic: the strain `(h_new - h)/h` cancels
  when the material streams, so one ulp grows to 1e-9. `lr1d.py` density and remap are
  now written as increments (exact identity when nothing moves); 1D results unchanged
  beyond round-off.

Next: WP 1.3b (rotation, plasticity, formulation decision, now also item 6: oblique
fronts).
