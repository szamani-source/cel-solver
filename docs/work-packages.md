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
| 2.6 | Phase 2 benchmarks | Taylor bar within tolerance |
| 3.0 | Silo benchmark definition and reference values | Janssen, Beverloo and column-collapse values sourced |
| 3.1p | Prototype gravity, walls and settling | Settled tank follows Janssen |
| 3.2p | Prototype free surface (material plus void) | Lake at rest, translating and falling blocks, separation |
| 3.2q | Prototype column collapse (split from 3.2p, decision 0021) | Column collapse run-out, top and deposit slope within decision 0019 |
| 3.2r | Prototype column collapse deposit slope (split from 3.2q, decision 0023) | Deposit body without whole-cell steps (done, decision 0024); toe split to 3.2s |
| 3.2s | Prototype collapse toe: stress in partly filled cells (split from 3.2r, decision 0024) | Same test as 3.2q: run-out, top and deposit slope within decision 0019 |
| 3.3p | Prototype silo discharge (target T0) | Steady rate, D^1.5 law, rate independent of head |
| 3.1 | C++ gravity and wall boundaries (smooth, frictional) | celFoam settles a tank, matches the prototype |
| 3.2 | C++ void and free surface | celFoam passes the 3.2p checks |
| 3.3 | C++ outlet boundary that opens | Material leaves only after the opening time, mass balance exact |
| 3.4 | Silo discharge in celFoam (target T0), serial and parallel | Matches the prototype and the reference values |
| 3.5 | Cohesion in the silo (arching, flow or no flow) | Flow stops below a critical outlet width |
| 3B.1 | N-material volume fractions with isoAdvector | Pure advection keeps interfaces sharp |
| 3B.2 | Mixed-cell closure | Two-material impact conserves mass and momentum |
| 3B.3 | Per-material history remap | Plastic strain stays with its material |
| 4.1 | Rigid body: STL input and kinematics | Body moves on prescribed path |
| 4.2 | Contact search near the body | Cut cells found, unit tested |
| 4.3 | Penalty forces and reaction output | No leakage in a push test |
| 4.4 | Footing and pile benchmarks | Bearing capacity factor near 5.14 |
| 5.1 | Hex8 element with hourglass control | Patch test passes |
| 5.2 | Lagrangian-only Taylor bar | Matches published result |
| 5.3 | Two-way deformable coupling | Deformable push test stable |
| 5.4 | Projectile benchmark | Matches reference run |
| 6.1 | MPI with Lagrangian bodies | 1 and N cores agree |
| 6.2 | Restart and mass scaling | Restarted run matches continuous run |
| 6.3 | Scaling tests and validation report | Report published |

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
Phases 4 to 6 (Lagrangian bodies, contact, FE, parallel scaling) are unchanged but come
after T0. 2.6 (Taylor bar) waits for Phase 5.
