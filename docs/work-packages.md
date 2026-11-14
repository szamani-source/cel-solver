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
| 2.3 | Stress update with objective rate | Rigid rotation produces no spurious stress |
| 2.4 | Single-material remap of all fields | 1D wave reproduced in 3D mesh |
| 2.5a | Drucker-Prager in C++ | Matches the 1.4 reference on all paths |
| 2.5b | Mohr-Coulomb in C++ (corners, apex) | Matches the 1.4 reference on all paths |
| 2.5c | Mohr-Coulomb tension cut-off (Python reference, then C++); split from 2.5b, decision 0011 | Cut-off returns match the reference; uniaxial tension stops at the cut-off |
| 2.6 | Phase 2 benchmarks | Taylor bar within tolerance |
| 3.1 | N-material volume fractions with isoAdvector | Pure advection keeps interfaces sharp |
| 3.2 | Void cells | Free surface holds with zero stress |
| 3.3 | Mixed-cell closure | Two-material impact conserves mass and momentum |
| 3.4 | Per-material history remap | Plastic strain stays with its material |
| 3.5 | Phase 3 benchmarks | Column collapse run-out within tolerance |
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
