# 0016 Lagrange step in C++: face states, boundaries, time step (WP 2.3)

Date: 2026-11-14. Status: accepted (Sia's standing rule: proceed and merge on green CI).

Number 0015 is taken by WP 2.5c (Mohr-Coulomb tension cut-off).

## Context

WP 2.3 implements `eulerianSolid::lagrangeStep(dt)` to decision 0013 (unsplit,
cell-centred, acoustic face Riemann solver, Hughes-Winget rotation, one model call per
cell) and hands the remap what decision 0014 asks for. Decision 0013 left open how to
build second-order face states on a general `fvMesh`, which boundary condition each
patch type means, how coupled faces stay identical on both sides, and the Courant
number. All numbers below are from `tests/test_lagrange_step.sh` (rho 2000 kg/m^3,
E 20 MPa, nu 0.3, strain 1e-6) or from the measurements recorded in `docs/wp/WP-2.3.md`.

## Decision

### 1. Face states: per-face MC along the cell-centre line, Hancock predictor

- Gauss cell gradients (linear interpolation) of the nine primitive components
  (U and sigma). On each side of a face, the state is the cell value plus half the
  MC-limited slope `minmod(2a, d.g, 2(2 d.g - a))` along `d`, the vector between the
  two cell centres, with `a` the jump across the face. On a uniform mesh this is exactly
  the 1D MC MUSCL slope of the prototypes, and needs no division.
- Half-step (Hancock) predictor per cell, from the cell's own face states: U gets
  `dt/2 div(sigma)/rho`, sigma gets `dt/2 C:sym(grad U)` with the elastic `C` (as the
  prototype, decision 0008), added to every face state of the cell.
- Limiting is on primitive components, not on characteristics as in the 1D prototype.
  It gives orders 1.9 to 2.3 (L1, 25/50/100 cells) for P and S waves along x and at
  45 degrees, speed errors below 0.13 percent; the prototype got about 1.8.
- Non-coupled boundary faces reconstruct through a mirror ghost cell at `2 (Cf - C)`
  holding `2 q_b - q_P`, `q_b` the field's patch value. A linear velocity with exact
  `fixedValue` patch values is then reconstructed exactly (rigid rotation keeps its
  principal stresses to 3e-13); on `zeroGradient` patches the face state is the cell value.
- `secondOrder no;` in `system/controlDict` switches to first order (order 0.92).

### 2. Boundary conditions from the patch type of U

| Patch / U patch type | Face velocity | Face traction |
| --- | --- | --- |
| `cyclic`, `processor`, `processorCyclic` | Riemann solver | Riemann solver |
| `empty` | none (no faces) | none |
| `symmetry`, `symmetryPlane` | normal 0, tangential from the interior characteristic | shear 0, normal from the interior characteristic |
| U fixes its value (`fixedValue` and derived) | the patch value | from the interior characteristic, `t* = tL + z (v* - vL)` (P impedance normal, S tangential) |
| anything else (`zeroGradient`, `calculated`, ...) | from the interior characteristic | zero (free surface) |

`wedge` and other coupled types (AMI, ...) stop with a fatal error. A step P-pulse
reflects from a fixed wall, a free surface and a symmetry plane, and an S-pulse from a
symmetry plane, with the right sign and amplitude to 1e-4.

### 3. Coupled faces are solved on the owner side only

For `cyclic` and `processor` patches the owner side (`cyclicPolyPatch::owner()`,
`processorPolyPatch::owner()`, which `processorCyclic` takes from its cyclic) runs
the Riemann solver; the other side receives the face velocity as is and the face force
and `phiL` negated (cyclic: copied locally; processor: one message per neighbour
processor, keyed by patch tag). Both sides hold bit-identical values whatever the
rounding of their own face geometry, so momentum crosses processor boundaries exactly.
Serial and 2 or 4 processors agree to 3e-15.

### 4. Default maxCo 0.4

`maxCo` keeps its meaning from WP 2.1: `dt = maxCo * min(cbrt(V)/c_p)`, the
Courant number along one mesh direction on a uniform mesh. Measured with a step P-wave
(45 degrees in 2D, 50 x 50; body diagonal in 3D, 24^3), periodic, steps of 0.05:

| | stable over 2 crossings (brief's test) | over 8 crossings | first unstable |
| --- | --- | --- | --- |
| 2D, 45 degrees | 0.65 | 0.60 and 0.65 bounded (slow creep of 1 to 2 percent) | 0.70 (explodes in the first crossing) |
| 3D, body diagonal | 0.50 | 0.45 bounded; 0.50 grows after 6 crossings | 0.55 |

The default is 0.4, one step below the 3D value that stays bounded over 8 crossings, so
it also covers mildly distorted cells. The Hancock predictor's cross terms make the
unsplit scheme more stable than the "about 0.3 in 2D" decision 0013 expected from the
split proxy. `tests/test_lagrange_step.sh` reads the default from `eulerianSolid.C`
and asserts no growth at it in 2D and 3D. Overshoot of the 45 degree step at 0.4:
0.46 percent (2D; the split prototype had 2.7 percent at start-up), 2.4 percent for the
3D body-diagonal step.

### 5. Other choices

- `solveMomentum no;` (controlDict) holds U fixed for prescribed-motion tests; the
  Riemann solver still gives the face velocities.
- The kernel `Model` gets `shearModulus()` (the S impedance needs G; the wrapper
  exposes it). Each model returns its `Elastic::G`.
- A cell whose strain increment is exactly zero is rotated but not passed to the model.
- The plastic strain field `epsP` stores tensor shear components; the kernel history
  uses engineering shear, converted at the call.
- The cell volume of the Lagrange step uses the fixed mesh volume `V` for the mass
  and `V_mid = (V + V_L)/2` for the velocity gradient, as the prototype's `h_mid`.

### 6. The WP 2.1 skeleton test changes

`tests/test_eulerian_skeleton.sh` case B (half the smokeBox moving at 1 m/s inside
fixed walls) asserted that the time step, mass and momentum stay constant: true only
while the Lagrange step was empty. With the step in, the walls stop the moving half
(momentum changes, physically) and, with no remap yet, the printed mass uses the
Lagrangian density. The test now checks case B's start values against the hand
calculation, compares every printed value of the 2 processor run with serial (it
compared only the last), and keeps the full constancy checks for case A (at rest).
Nothing it checked about field I/O, defaults and parallel agreement is dropped.

## Consequences

- The remap (WP 2.4) receives `phiL` and the Lagrangian state as decision 0014.
  WP 2.4b re-runs rotation and streaming checks with both steps.
- The test driver `tests/lagrange/lagrangeStepCheck` keeps the cell mass fixed after
  each step (a "null remap", linear acoustics), so it measures the Lagrange step alone.
- `stableDeltaT` uses `cbrt(V)`: a 2D one-cell-thick mesh must be as thick as its
  in-plane cells, or the time step is wrong. The tests do this. A face-based Courant
  number (sum of `|Sf|` over non-empty faces) would remove the restriction; left for a
  later package because it changes `maxCo`'s meaning and the WP 0.1 and 2.1 tests.
- The Courant limit on distorted or graded meshes is not measured.
