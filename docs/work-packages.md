# Work packages

Each package fits in one to three working sessions and is done when its test passes.
Phases and gates are explained in the plan doc (link in `STATUS.md`).
Each package gets a brief in `docs/wp/` before work starts (see `docs/MANAGEMENT.md`).
Constitutive models for bulk solids are tracked in `docs/constitutive-models.md`.

| ID | Work package | Done when |
| --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Empty solver builds and runs in CI |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Files exist and describe the plan |
| 0.3 | Target problems, units, material data | Spec file agreed with Sia |
| 0.4 | Fast OpenFOAM setup for fresh cloud sessions | New session builds and tests in under 10 minutes |
| 1.1 | 1D elastic Lagrange-plus-remap prototype | Converges at second order on a smooth wave |
| 1.2 | 1D plasticity and history remap | Elastic-plastic wave matches analytic |
| 1.3a | 2D plane-strain prototype, elastic waves | P and S waves at 0 and 45 degrees converge; reduces to 1D |
| 1.3b | 2D rotation, plasticity and formulation decision | Rotating disc keeps its stress; decision note written |
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Triaxial strengths match closed form |
| 2.1 | Eulerian library skeleton and solver fields | Runs empty cycles serial and on 2 processors, fields written |
| 2.2 | Constitutive library interface, elastic and J2 | Single-point unit tests pass |
| 2.3 | Lagrange step in C++: face Riemann solver, objective stress update | Waves, rotation, shear, walls, parallel checks pass |
| 2.4 | Remap in C++, single material | Reduces to 1D; conservative, monotone, second order |
| 2.4b | Lagrange step plus remap, integration checks | 2D prototype checks reproduced by celFoam |
| 2.5a | Drucker-Prager in C++ | Matches the 1.4 reference on all paths |
| 2.5b | Mohr-Coulomb in C++ (corners, apex) | Matches the 1.4 reference on all paths |
| 2.5c | Mohr-Coulomb tension cut-off (Python reference, then C++); split from 2.5b, decision 0011 | Cut-off returns match the reference; uniaxial tension stops at the cut-off |
| 2.6 | Phase 2 benchmarks. Done 2026-11-18 on the wedge, decision 0063 (`tests/test_taylor_wedge.sh`) | Taylor bar within tolerance (shortening 5 %, widening 20 % against the axisymmetric Abaqus benchmark) |
| 3.0 | Silo benchmark definition and reference values | Janssen, Beverloo and column-collapse values sourced |
| 3.1p | Prototype gravity, walls and settling | Settled tank follows Janssen |
| 3.2p | Prototype free surface (material plus void) | Lake at rest, translating and falling blocks, separation |
| 3.2q | Prototype column collapse (split from 3.2p, decision 0021) | Column collapse run-out, top and deposit slope within decision 0019 |
| 3.2r | Prototype column collapse deposit slope (split from 3.2q, decision 0023) | Deposit body without whole-cell steps (done, decision 0024); toe split to 3.2s |
| 3.2s | Prototype collapse toe: stress in partly filled cells (split from 3.2r, decision 0024) | Same test as 3.2q: run-out, top and deposit slope within decision 0019 (slope as a fitted line, decision 0025; done) |
| 3.3p | Prototype silo discharge (target T0) | Outlet opens, no leak before, mass balance exact, rate independent of head (rate law split to 3.3q, decision 0026) |
| 3.3q | Prototype silo discharge rate (split from 3.3p, decision 0026) | Steady rate, D^1.5 law with C 1.0 to 1.9, independent of E and mesh |
| 3.3r | Prototype discharge rate at T0 stiffness (split from 3.3q, decision 0027) | Rate at E = 0.1, 1, 10 MPa within 10 percent, no noise in the dead zones |
| 3.2t | One rule set for collapse and discharge (split from 3.3r, decision 0028) | Collapse test passes with rules 12 and 13 on, run-out at 1 and 10 MPa within 0.1 |
| 3.1 | C++ gravity and wall boundaries (smooth, frictional) | celFoam settles a tank, matches the prototype |
| 3.2a | C++ void and free surface (split from 3.2) | celFoam passes the 3.2p checks (`tests/test_void.sh`) serial and on 2 processors |
| 3.2b | C++ column collapse: shear cap (rule 8), rule 9, neither rule 12 nor 13 (decision 0029; split from 3.2) | Column collapse run-out at E = 0.1, 1, 10 MPa recorded against decision 0019 |
| 3.3 | C++ outlet boundary that opens | Material leaves only after the opening time, mass balance exact |
| 3.4 | Silo discharge in celFoam (target T0), serial and parallel | Matches the prototype and the reference values; rate at E = 0.1, 1, 10 MPa within 10 % with the same code as the collapse (settles rules 12 and 13, decision 0029) |
| 3.4b | C++ discharge rate independent of E (split from 3.4, decision 0034) | `tests/manual/test_discharge_stiffness.sh` passes at 10 mm and at 5 mm (0.1, 1 MPa) while `tests/test_collapse.sh` still passes |
| 3.6a | mu(I) friction (Barker-Gray regularised where a well-posed range exists, decision 0037) in the constitutive library, Python then C++ (decision 0036) | Single-point simple shear at fixed pressure gives mu(I) to 1e-6 over I = 1e-4 to 1; mu_2 = mu_s reproduces the present MC/DP paths bit for bit; C++ matches Python |
| 3.6b | Steady inclined layer (Bagnold profile) in celFoam with mu(I), E = 0.1, 1, 10 MPa | Surface speed and profile within 5 % of the exact solution at 0.1 MPa; the error against E measures the numerical viscosity (recorded either way) |
| 3.6s | Mesh dependence of the collapse shear band (split from 3.6b, 2026-11-17) | `tests/manual/test_shear_band.sh` runs a = 0.5 at 10, 5, 2.5 mm with MC and mu(I), mass kept, widths and run-outs recorded (decision 0038) |
| 3.6c | Plastic-wave impedance at flowing faces in the Riemann solver (replaces rules 8 and 12 if it matches them). Closed 2026-11-17 without its done-when; rest split to 3.6c2 (decision 0040) | 3.6b profile within 5 % at all three E, `tests/test_collapse.sh` still passes |
| 3.6c2 | Inclined layer at 10 MPa and the top cell (split from 3.6c, decision 0040): floor cell stuck and slow body at 10 MPa; top cell next to void slides (stress-free at 10 MPa). Done 2026-11-17 (decision 0044) | 3.6b profile within 5 % at all three E, 5 mm, `tests/test_collapse.sh` still passes |
| 3.6d | Silo discharge with mu(I) and the new impedance. Closed 2026-11-17 without its done-when (1 MPa about 20 % slow; decision 0045) | `tests/manual/test_discharge_stiffness.sh` passes at 10 mm and 5 mm (the 3.4b done-when), C within decision 0019 |
| 3.6f | Why the silo drains about 20 % slower at 1 MPa (split from 3.6d). Closed 2026-11-18, finding recorded (decision 0046) | Each E as a band of three runs; the cause named with runs that switch it on and off |
| 3.6g | Dilatancy in the mu(I) law (compressible mu(I)/phi(I), decision 0046 proposal). Closed 2026-11-18 without its done-when (decision 0048) | Silo bands at 1 and 10 MPa within 10 % of 0.1 MPa at 10 mm, inclined layer still passes; kernel `tests/test_dilatancy.sh` |
| 3.6h | Smooth blend of the full and partly filled rules (proposed in decision 0048). Closed 2026-11-18 without its done-when (decision 0049) | Silo bands within 10 % at 10 mm, inclined layer still passes |
| 3.6e | Column collapse with mu(I) against experiments (Lube 2005, Lajeunesse 2005). Done 2026-11-18 at 10 mm (decision 0050: T0 rule set mu(I) + flowImpedance + normalCap) | Run-out and final height within the decision 0019 bands at a = 0.5 and 2, at T0's real E |
| 3.5 | Cohesion in the silo (arching, flow or no flow) | Flow stops below a critical outlet width |
| 3.7a | Solid inlet: a patch that feeds fresh material (speed over time, density, stress), decision 0060. Done 2026-11-18 | `tests/test_inlet.sh`: nothing in while closed, fed mass exact, mass balance 1e-12, the fed material lands and comes to rest, 2 processors match serial |
| 3.7b | Silo fed from the top while it discharges (tutorial `siloFed`, first with `inlets`). Done 2026-11-18 | `tests/test_fed_silo.sh`: outflow with a feed of R/2 and R within 10 % of the unfed R; matched feed holds the tank mass within 3 %; tank = initial + fed - discharged to 1e-10 |
| 3.8a | Axisymmetric wedge: wedge patches in the Lagrange step and remap, time step from the in-plane cell size (no free surface), decision 0061. Done 2026-11-18 | `tests/test_wedge.sh`: Lamé thick cylinder within 1 % at 2.5 mm beyond two wall cells (wall cells first order); lake at rest with axis cells; cylindrical Janssen tank within 5 % of Janssen for a circular bin (was 2 %; 5 % as the Cartesian check, measured 2.0 %); 2 processors match serial |
| 3.8b | Free surface on the wedge: hoop part of the radial flux in the void closure and slab lengths, decision 0061 addendum. Done 2026-11-18 | `tests/test_wedge_void.sh`: lake at rest, falling plug and T0 void-top tank (within 5 % of circular Janssen) on a wedge, mass 1e-12, 2 processors match serial (the done-when said "as `tests/test_void.sh`"; its translating and separating blocks are Cartesian-only motions) |
| 3.8c | Cylindrical silo discharge on the wedge, central circular outlet. Done 2026-11-18, decision 0062 (outlet D 0.06 m: the quarter needs 3 cells across) | Rate within the 3D Beverloo band (published C, k); same tank as a coarse 3D quarter model within 10 % |
| 3.8d | STL walls on a wedge (conical hopper), later | Conical hopper discharges; rate recorded against Rose-Tanaka |
| 3B.1 | N-material volume fractions with isoAdvector | Pure advection keeps interfaces sharp |
| 3B.2 | Mixed-cell closure | Two-material impact conserves mass and momentum |
| 3B.3 | Per-material history remap | Plastic strain stays with its material |
| 3B.4 | General meshes (polyhedral, snappyHexMesh): PLIC void on general cells, void closure along face normals, face extrapolation to the face centre, time step from cell inradius (not planned, decision 0041) | Silo settle and outlet tests pass on a polyhedral mesh |
| 4.0 | Investigate Lagrangian-Eulerian contact schemes (LS-DYNA constraint and penalty coupling, Abaqus CEL contact, cut-cell, immersed boundary) | Sourced comparison, prototype experiments, recommendation agreed with Sia |
| 4.1 | Rigid body: STL input and kinematics (decision 0041) | Body moves on prescribed path |
| 4.2 | Contact search near the body | Cut cells found, unit tested |
| 4.3a | Fixed STL walls with Coulomb friction (cut-cell moving walls, decision 0042) | Janssen tank with STL walls, aligned and tilted 10 deg, within 2 / 5 % of mesh-boundary walls, no leakage |
| 4.3c | STL walls with the free surface and outlets (split from 4.3a) | Void and slab rules at cut cells; a block falling onto a tilted STL floor stops without leaking; silo settle with STL walls and void top |
| 4.3b | Moving bodies: coupling forces both ways and reaction output (decision 0042) | No leakage in a push test |
| 4.3d | Free bodies: motion from the material force and gravity (split from 4.3b) | Dropped block comes to rest carrying its weight, momentum kept |
| 4.3e | Rough walls: wall adhesion (split from 4.4b) | Rough pile within 15 % of 11.94 at 10 mm |
| 4.5 | Silo and hopper discharge with STL walls | Flat tank rate within 5 % of mesh-boundary walls; hopper runs and its rate is recorded |
| 4.4 | Strip footing benchmark (split 2026-11-18: pile to 4.4b) | Bearing capacity factor near 5.14 |
| 4.4b | Pile section moving sideways through clay (Randolph and Houlsby 1984) | Smooth limit 9.14 cu D within 10 % |
| 4.4c | Strip footing on sand (Martin 2005 N_gamma, phi 30) | Smooth N_gamma 7.653 within 15 % (closed without it, decision 0057) |
| 5.1 | Hex8 element with hourglass control | Patch test passes |
| 5.2 | Lagrangian-only Taylor bar | Matches published result |
| 5.3a | Two-way deformable coupling (split 2026-11-18) | Steel footing push follows the rigid one within 5 %, density within 2 % |
| 5.3b | Soft deformable bodies | Soft footing shortens as an elastic column within 20 %, density within 2 % |
| 5.4 | Projectile benchmark | Matches reference run |
| 6.1 | MPI with Lagrangian bodies | 1 and N cores agree (decision 0080: exact phase, then outcomes within the round-off band) |
| 6.2a | Restart (split 2026-11-18) | Restarted run matches continuous run, bit for bit |
| 6.2b | Mass scaling of deformable bodies (split from 6.2; later, on request) | Quasi-static push with scaled elements gives the unscaled force within 5 % |
| 6.3a | Parallel scaling (split 2026-11-18) | Speed-up on 1, 2, 4 processors recorded for the T0 tank at 10 and 5 mm |
| 6.3b | Validation report (split 2026-11-18) | Report published: every benchmark with reference, result, tolerance, test, known issues |

Note: the skeleton written in 0.1 already covers part of 2.1 (explicit loop and CFL step).

Note (2026-11-14): 1.3 split into 1.3a and 1.3b to keep each within three sessions.

Note (2026-11-14, WP 1.3b): decision 0013 fixes the C++ formulation for 2.3 and 2.4:
cell-centred, unsplit face Riemann solver, Jaumann rate by Hughes-Winget, one
constitutive call per cell per step, per-field remap limiters (MC; superbee for
history). 2.3's rigid-rotation test is the WP 1.3b disc check; 2.3 also re-runs the
simple-shear and 45 degree impact checks and measures its stable Courant number.

Note (2026-11-14, decision 0016): Sia's primary target is 2D silo discharge (T0 in
`docs/spec/target-problems.md`). Phase 3 is now the silo: packages ending in "p" are
Python prototypes that start at once; 3.1 to 3.4 are their C++ counterparts and follow
WP 2.3, 2.4 and 2.4b. The old 3.1, 3.3, 3.4 (several materials per cell) are now 3B;
the old 3.2 (void) is 3.2p and 3.2; the old 3.5 (column collapse) was inside 3.2p and
is 3.2q since 2026-11-14 (split: the run-out depends on the Riemann solver's dissipation,
decision 0021). 3.2q fixed the run-out (face shear cap, decision 0023); the deposit
slope is split to 3.2r (2026-11-14). 3.2r removed the whole-cell staircase (decision
0024); the toe, which needs partly filled cells to carry stress, is split to 3.2s
(2026-11-14).
3.3p ran T0 end to end with an outlet; its rate is set by the Riemann solver's
dissipation, so the rate law is split to 3.3q (2026-11-15, decision 0026).
3.3q found the low rate was an outlet bug (a partly filled outlet cell shut the face,
decision 0027); with it fixed the rate law holds at 5 mm and E = 0.1 MPa, but the rate
still falls with E (81 / 73 / 54 kg/s/m at 0.1 / 1 / 10 MPa, 10 mm), split to 3.3r
(2026-11-15).
3.3r made the rate independent of E (normal jump cap and dispersed void, rules 12 and
13, decision 0028), but the rules make the column collapse run too far, so they are on
only in the discharge; reconciling them is split to 3.2t (2026-11-15).
3.2t found no single rule set in the split prototype (decision 0029): every variant of
rule 12 that keeps the 10 MPa collapse stopping also slows the discharge. Closed without
its done-when; the C++ silo packages start with the full jump terms and measure both
cases with one code (2026-11-15).
3.4b (decision 0035) added rule 12 in C++ (`normalCap`): 10 mm met, 5 mm still -16 %, closed without its done-when; Sia to choose E = 0.1 MPa or mu(I). 3.4 found the C++ rate falls with E with the collapse code (10 mm: 74 / 66 / 71 kg/s/m
at 0.1 / 1 / 10 MPa; 5 mm: 78 / 55 at 0.1 / 1 MPa) and parallel matches serial; closed
without its done-when, the E independence split to 3.4b (2026-11-16, decision 0034).
Phases 4 to 6 (Lagrangian bodies, contact, FE, parallel scaling) are unchanged but come
after T0. 2.6 (Taylor bar) done on the wedge (decision 0063).

Note (2026-11-15): WP 4.0 added at Sia's request: a thorough comparison of contact
schemes before any coupling code is written. It is the first package of Phase 4 and
can be started earlier if wanted, since it needs only the 2D prototype. Tank walls are
mesh boundaries, not contact, so the silo track does not depend on it.

Note (2026-11-15): 3.2 split into 3.2a (void, decision 0031) and 3.2b (shear cap and
column collapse) to keep each within three sessions.

Note (2026-11-17): Sia asked for wall interactions from STL and about mesh types.
Decision 0041: the Eulerian mesh stays Cartesian hex and STL walls are immersed rigid
surfaces (as in Abaqus CEL); polyhedral or snappyHexMesh meshes would need a PLIC void
on general cells first (3 to 4 packages, not planned). 4.3 split into 4.3a (fixed STL
walls with friction, what Sia asked for) and 4.3b (moving bodies); 4.5 added (silo and
hopper with STL walls). Order: 4.1 (done), 4.0 and 4.2 (may run side by side), 4.3a,
4.5, 4.3b, 4.4.

Note (2026-11-18): WP 4.3a done without the free surface (decision 0047). 4.3c added:
STL walls with void and outlets, which 4.5 needs. Order: 4.3c, 4.5, 4.3b, 4.4.

Note (2026-11-18): WP 4.3c done (decision 0051, tutorial `stlSilo`). Its tank check
compares with the mesh-wall tank with void (not the lid version): the free-surface tank
itself is up to 17 % off the no-void one (STATUS known issues). Next: 4.5, 4.3b, 4.4.

Note (2026-11-18): WP 4.5 done (decision 0052, tutorial `stlHopper`): STL tank +4.4 % against
mesh walls (bands of three runs at 0.1 MPa); hopper rate recorded. Next: 4.3b, 4.4.

Note (2026-11-18): WP 4.3b done (decision 0053): prescribed-motion bodies push the material;
free bodies (motion from the force) are left for when a case needs them. Next: 4.4.

Note (2026-11-18): WP 4.4 split: the strip footing is 4.4 (done, decision 0055, Prandtl
within 8 % at 10 mm and 1.9 % at 5 mm); the pile is 4.4b: the lateral limit of Randolph and Houlsby (1984), done (decision 0056).

Note (2026-11-18): Sia asked for two things Abaqus CEL cannot do and asked for both
("Please make these two possible in the solver"): a solid inlet (3.7a, 3.7b) and
axisymmetric runs on a wedge so a cylindrical silo needs no quarter model (3.8a to 3.8d).
A wedge mesh keeps every internal face along r or z, so the Cartesian rule of decision
0031 holds inside the mesh; only the two wedge sides are tilted and no material crosses
them. The work is in the wedge faces (rotated cell state, hoop terms), the time step
(cbrt(V) collapses at the axis), the slab rules (cells widen with r, triangles at the
axis) and the outlet on the axis. Order: 3.7a, 3.7b, 3.8a, 3.8b, 3.8c; 3.8d later.


Note (2026-11-18, phase 6): 6.2 split into 6.2a (restart) and 6.2b (mass scaling, kept
for later: the material's step, not the elements', sets the cost of the coupled runs,
since deformable bodies are subcycled, decision 0072), and 6.3 into 6.3a (scaling) and
6.3b (validation report). Phase 6 decisions are numbered from 0080.
