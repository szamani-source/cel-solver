# 0013 Eulerian formulation for the C++ solver (WP 1.3b)

Date: 2026-11-14. Status: accepted (Sia's standing rule: proceed and merge on green CI).
Revisit item 2 when WP 2.3/2.4 have run the same checks on the C++ unsplit scheme.

Numbers 0009 to 0012 are taken by WP 2.5a, 2.1, 2.5b and 1.3a.

## Context

WP 1.3b adds rotation and plasticity to the 2D prototype (`prototype/lr2d_plastic.py`)
and has to fix the formulation that the C++ Lagrange step (WP 2.3) and remap (WP 2.4)
implement. The brief asks six questions. All numbers below come from
`tests/test_proto_2d_rotation.sh` (`lr2d_plastic.py --verify`), from
`tests/test_proto_2d_elastic.sh`, or from `lr2d_plastic.py --study` (measurements, not
pass criteria). Material: rho 2000 kg/m^3, E 20 MPa, nu 0.3; von Mises Y 10 kPa;
Drucker-Prager phi 30, c 5 kPa, psi 0 (compression fit).

### A correction to decision 0012 first

Decision 0012 reported that a step P-wave at 45 degrees overshoots by 13 to 14 percent of
the jump at start-up and 5 to 7 percent at the end, and put that down to the
direction-by-direction Riemann solver. Most of it was the test, not the scheme. The step
was set with `s = x + y` computed from cell centres in floating point, so on about 1900
of the 10000 cells of a 100 x 100 box the threshold `s >= 0.2` fell on the wrong side for
some cells of a diagonal and not for others. The initial front was jagged at the cell
scale, and a jagged front radiates real 2D waves (tangential velocity up to 7 percent).
With the step set from integer indices (every cell of a diagonal on the same side):

| 45 degree step P-wave (lr2d.py) | 100 x 100 | 200 x 200 |
| --- | --- | --- |
| Jagged step (decision 0012): over the run / at the end | 14.1 / 6.6 % | 13.1 / 4.9 % |
| Clean step: over the run / at the end | 2.7 / 0.20 % | 2.7 / 0.12 % |
| Clean step, after 3 cells of travel / after 10 | | 0.28 / 0.22 % |

So an oblique step behaves like the 1D impulsive start-up of decision 0008 (1.7 to
2.7 percent in the first cells, then below 1 percent). `lr2d.check_step` now sets the
step from integer indices and checks the 1 percent criterion at 45 degrees after a
start-up of 10 cells (it was reported only); the along-x check is unchanged.

## Decision

### 1. Cell-centred Godunov, not node-based staggered

Keep every field at the cell centre, with an acoustic Riemann solver per face
(decision 0005).

- Oblique fronts are fine (item 6), so the one argument left for a nodal or staggered
  layout is gone.
- The staggered 1D version rang by 4 to 7 percent at a step with every limited
  artificial viscosity tried (decision 0005); the cell-centred one has no new extrema
  along an axis and 0.22 percent at 45 degrees after start-up.
- OpenFOAM's `fvMesh` is cell-centred: face loops, face fluxes, processor-boundary
  faces and isoAdvector all work per face. A staggered scheme needs point fields,
  point-to-cell loops and point synchronisation across processors, which is more of our
  own code to write and keep correct in parallel.
- Abaqus/Explicit CEL and LS-DYNA ALE use the staggered layout with bulk viscosity and
  accept its ringing; we are not bound to their layout.

### 2. Unsplit in C++; the prototype stays split

The C++ solver runs on general polyhedral `fvMesh` meshes, where there are no mesh
directions to split along. Its Lagrange step is unsplit: one Riemann solve per face from
the same state, all faces of a cell at once, one strain increment per cell per step,
one constitutive call per cell per step, then one remap.

Evidence that an unsplit face-based update keeps oblique fronts clean: in the prototype,
an "additive" proxy (x and y sweeps both from the same state, conservative increments
added, as a face-based unsplit update would) on the 45 degree elastic impact
(`--study`, item C):

| Elastic impact at 45 deg, 100 x 100, extra over along x | start-up | after start-up | at the end |
| --- | --- | --- | --- |
| Symmetric split (prototype), CFL 0.5 | 0.44 % | 0.24 % | 0.24 % |
| Unsplit proxy, CFL 0.5 | 3.85 % | 14.3 % (growing) | 12.9 % |
| Unsplit proxy, CFL 0.4 | 1.01 % | 1.34 % | 0.69 % |
| Unsplit proxy, CFL 0.3 | 0.80 % | 0.48 % | 0.45 % |

(Run to 0.5 of the period, longer than the test's 0.22, so the compression and tension
structures meet; hence 0.24 % for the split scheme here and 0.000 % in the test.)

The price of unsplit without corner coupling is the time step: the proxy grows
oscillations at CFL 0.5 and is clean at 0.3 in 2D. WP 2.3 sets the Courant number from
its own stability test (expect about 0.3 in 2D and lower in 3D), or adds corner
transport later if the cost matters. The proxy is not meaningful for plasticity (it
returns each sweep's increment to the yield surface separately and adds the results,
which a real unsplit update does not do), so only its elastic numbers are used.

### 3. Jaumann rate, integrated with the Hughes-Winget rotation

Stress (and the plastic strain tensor) is rotated with
`Q = (I - dW/2)^-1 (I + dW/2)`, `dW = dt skew(grad v)` at mid-step, then the model adds
the response to the strain increment `dt sym(grad v)`. In the prototype each sweep adds
its half of the spin (`dW_tn = d gamma_nt / 2`).

| Check (test_proto_2d_rotation.sh) | Result | Criterion |
| --- | --- | --- |
| Pre-stressed disc, rigid rotation, 1 turn, 100 cells across: principal stresses | 0.027 % | 2 % |
| Same: principal direction | 0.0006 deg | 1 deg |
| Same: stress more than 6 cells outside the disc | 0.015 % | 1 % |
| Without the rotation correction: direction error at a quarter turn | 90 deg | |
| Simple shear, elastic, against the Jaumann closed form G sin(gamma) | 2.1e-4 G | 1e-3 G |
| Simple shear, von Mises, to gamma 1: shear stress vs Y/sqrt(3) after yield | 3e-7 | 1 % |
| Same: largest drop (oscillation) | 2.5e-7 of k | 0.1 % |
| Simple shear, Drucker-Prager at p = -50 kPa: sqrt(J2) vs the cone | 2e-14 | 1 % |

Reasons: an Eulerian code keeps no deformation gradient, so the Green-Naghdi rate
(which needs the rotation R of F) is not available without carrying F as nine more
remapped fields. Jaumann's known flaw, the oscillating shear stress of a hypoelastic
material in simple shear (0.84 G at gamma 1 instead of G gamma), needs elastic shear
strains of order one. Soils and other bulk solids stay below about 1 percent elastic
strain, and with plasticity the von Mises and Drucker-Prager shear stresses stay on the
yield surface with no oscillation at all. Hughes-Winget is an exact rotation for any
step and costs a 3 x 3 product per cell.

### 4. Impulsive start-up: accept, no elastic-plastic Riemann solver now

Start-up overshoot is 1.7 percent (elastic) to 2.7 percent (von Mises) at a 1D impact
(decision 0008), 2.7 percent for a clean step at 45 degrees, and 0.4 to 0.8 percent extra
for an impact at 45 degrees over the same impact along x. It decays below 1 percent
within 3 to 10 cells of travel at every angle and plasticity tested. The C++ code does
not add an elastic-plastic Riemann solver or contact viscosity now. Contact (WP 4.3)
must measure the overshoot at a penalty contact in its own test and add a ramp or local
viscosity only if that test needs it.

### 5. Remap limiters per field

| Field | Remap | Limiter |
| --- | --- | --- |
| Density (mass) | volume flux | MC |
| Velocity | mass-weighted (mass coordinates) | MC |
| Stress (all components) | mass-weighted | MC |
| History: plastic strain tensor, equivalent plastic strain | mass-weighted | superbee |

As decisions 0005 and 0008, now confirmed in 2D. The C++ remap (WP 2.4) takes the limiter
as a per-field option with these defaults. A cell whose strain increment is exactly zero
keeps its state (the model is not called); this matters only for the split prototype,
where it makes a 1D problem on the 2D mesh match `lr1d_plastic.py` to 1e-11.

### 6. Oblique fronts: accept the face Riemann solver

No nodal (multidimensional) Riemann solver and no staggered grid. With the corrected
initial data the 45 degree step overshoots only during start-up (2.7 percent, then below
0.3 percent), and an elastic-plastic impact at 45 degrees adds at most 0.08 percent over
the same impact along x after start-up (von Mises 0.000, Drucker-Prager 0.082, elastic
0.000 percent; checked at 1 percent by test_proto_2d_rotation.sh). Using the P
impedance for the shear system (more dissipation, `--study` item D) does not help:
it makes the elastic start-up worse (1.31 against 0.44 percent) and gains nothing after
it (0.08 against 0.00 percent; von Mises 0.37 and 0.00, same as the elastic impedance).

## Consequences

- WP 2.3 (stress update) implements: face Riemann solver (acoustic, elastic impedances
  `rho c_p` normal, `rho c_s` tangential), strain and spin increments from the face
  velocities by the Gauss theorem, Hughes-Winget rotation, one model call per cell per
  step. Its "rigid rotation produces no spurious stress" test is the disc check above.
- WP 2.4 (remap) implements the per-field limiters above, mass-weighted for velocity,
  stress and history.
- The C++ time step is lower than the prototype's (unsplit, item 2). WP 2.3 measures it.
- Risk: the unsplit C++ scheme is not the scheme prototyped here. WP 2.3/2.4 re-run the
  disc, simple-shear and 45 degree checks on the C++ code; if oblique fronts ring there,
  revisit items 1, 2 and 6.
- Decision 0012 item 4 is corrected by this note.
