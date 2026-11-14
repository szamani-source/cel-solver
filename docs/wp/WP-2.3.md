# WP 2.3: Lagrange step in C++ (face Riemann solver, objective stress update)

State: not started
Branch: `wp-2.3-lagrange-step`

## Goal

`eulerianSolid::lagrangeStep(dt)` implements the scheme of decision 0013: unsplit,
cell-centred, acoustic Riemann solver per face, Hughes-Winget rotation, one constitutive
call per cell per step. It runs on any `fvMesh`, serial and parallel.

## Read first

- `docs/decisions/0013-eulerian-formulation.md` (the scheme), `0014-lagrange-remap-contract.md`
  (what this step must hand to the remap), `0010`, `0006`
- `prototype/lr2d_plastic.py` and `prototype/lr2d.py`: the reference for the sweep
  physics (impedances, characteristic slopes, strain and spin from face velocities)
- `src/celEulerian/eulerianSolid/`, `src/celConstitutive/constitutiveModel/constitutiveModel.H`
- `CLAUDE.md`: run `tools/setup_openfoam.sh` first

## Build

In `src/celEulerian/eulerianSolid/eulerianSolidLagrange.C` (new file; move the stub there):

1. Limited cell gradients of `U` and `sigma` (second order; first order as a switch).
2. Per face, left and right states at the face centre (half-step predictor as in the
   prototype if needed for second order in time), then the acoustic Riemann solution in
   the face-normal frame: normal velocity and normal traction with `rho c_p`, tangential
   velocity and shear traction with `rho c_s`. Result: face velocity `v_f` and face
   traction `t_f`.
3. `phiL_ = v_f & Sf` (decision 0014).
4. Momentum: `m dU = dt * sum(t_f |Sf|)` (cell mass constant in the step).
5. Velocity gradient by Gauss from `v_f`: strain increment `dt sym(grad v)`, spin
   increment `dt skew(grad v)`; rotate `sigma` and `epsP` with Hughes-Winget; call the
   model once per cell; a cell with zero strain increment keeps its state.
6. Lagrangian density: `rho = m / V_L`.
7. Boundary faces, chosen from the patch type of `U` (record the mapping in a decision
   note): fixed wall and prescribed velocity (`fixedValue`: face velocity given, traction
   from the interior characteristic), free surface (zero traction, face velocity from the
   interior characteristic), `symmetry`/`empty`/`processor` handled correctly. Processor
   faces must give bit-identical face states on both sides.
8. A dictionary switch to hold `U` fixed (`solveMomentum no;`) for prescribed-motion tests.

Tutorials and `tests/test_lagrange_step.sh`. Use one-cell-thick meshes (`empty` patches)
for 2D checks.

## Done when

`tests/test_lagrange_step.sh` passes in CI. Remap is still empty in this package, so
wave tests use strain amplitude 1e-6 (decision 0014) and the rest use prescribed motion:

- P-wave and S-wave, along x and at 45 degrees, periodic (`cyclic`) box: speeds within
  1 percent; smooth pulse converges at order 1.7 or better over 3 meshes.
- Measured stable Courant number in 2D and in 3D (largest `maxCo` in steps of 0.05 with
  no growth of a 45 degree, respectively body-diagonal, step wave over 2 box crossings);
  the test asserts stability at the value then written as the default `maxCo`.
- Uniform pre-stress under prescribed rigid rotation, one full turn: principal stresses
  within 0.1 percent, direction within 0.1 degree.
- Uniform simple shear, von Mises, to shear strain 1: shear stress at Y/sqrt(3) within
  0.1 percent after yield, no oscillation. Same with Drucker-Prager on its cone.
- Step P-wave along x reflecting from a fixed wall and from a free surface: reflected
  amplitude and sign correct within 2 percent.
- Total momentum conserved to 1e-12 relative in the periodic cases.
- The 45 degree wave case on 2 and on 4 processors matches the serial fields to 1e-12.
- Existing tests still pass unchanged.

Record in the log: run time per cell per step, and the overshoot of a 45 degree step
(compare decision 0013: 2.7 percent at start-up, below 0.3 percent after).

## Out of scope

- `remap` (WP 2.4, running in parallel in `eulerianSolidRemap.C`; both touch
  `Make/files` and `docs/STATUS.md`). Large-motion checks that need both (rotating disc,
  streaming waves, plastic impact): WP 2.4b.
- More than one material, void, contact, gravity.

## If it does not fit in three sessions

Split and record it: 2.3a interior scheme and periodic tests; 2.3b boundary conditions
and parallel.

## Log

(none yet)
