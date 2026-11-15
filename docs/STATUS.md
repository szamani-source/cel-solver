# Status

Last updated: 2026-11-15 (WP 3.3p)

Plan doc: https://claude.ai/artifact/DTT6kSikz9ZUTBkUxY6Grr

## Current phase

Phase 2 (the C++ Eulerian solver: 2.1 to 2.5c and 2.4b done; 2.6 waits for Phase 5) and Phase 3 (silo discharge, Sia's primary target, decision 0016) in parallel.

## Work packages

| WP | Title | State | Notes |
| --- | --- | --- | --- |
| 0.1 | Repository, build against pinned OpenFOAM, CI | Done | Repository szamani-source/cel-solver; CI builds celFoam in opencfd/openfoam-dev:2512 and test_smokeBox passes (2026-11-14) |
| 0.2 | CLAUDE.md, STATUS.md, decisions folder | Done | |
| 0.3 | Target problems, units, material data | Done | Sia named the primary target (2026-11-14): 2D silo discharge, cohesionless, phi 35 deg. `docs/spec/target-problems.md` (T0). Tank dimensions are placeholders |
| 0.4 | Fast OpenFOAM setup for fresh cloud sessions | Done | `tools/setup_openfoam.sh` installs the CI image's OpenFOAM from registry layers; setup + build + tests in 66 s (decision 0007). PR #5 merged (2026-11-14) |
| 1.1 | 1D elastic Lagrange-plus-remap prototype | Done | `prototype/lr1d.py`, `tests/test_proto_1d_elastic.sh` passes (order 1.88/1.83, monotone, conservative). Cell-centred instead of staggered: decision 0005. PR #3 merged (2026-11-14) |
| 1.2 | 1D plasticity and history remap | Done | `prototype/lr1d_plastic.py`, `tests/test_proto_1d_plastic.sh` passes (wave speeds and HEL within 0.2 %, history edges 2 cells). Decision 0008 accepted. PR #6 merged (2026-11-14) |
| 1.3a | 2D plane-strain prototype, elastic waves | Done | `prototype/lr2d.py` (symmetric dimensional splitting), `tests/test_proto_2d_elastic.sh` passes: P and S waves at 0 and 45 degrees, orders 1.77 to 1.87, speeds within 0.002 %; matches `lr1d.py` bit-for-bit; mirror-exact. PR #10 merged (2026-11-14). The 45 degree overshoot was a jagged test step (decision 0013) |
| 1.3b | 2D rotation, plasticity and formulation decision | Done | `prototype/lr2d_plastic.py` (Jaumann via Hughes-Winget, WP 1.4 models vectorised in `constitutive/batch.py`); `tests/test_proto_2d_rotation.sh` passes (disc 0.03 %, simple shear 3e-7, 1D match 1e-11). **Formulation: `docs/decisions/0013-eulerian-formulation.md`**: cell-centred, unsplit in C++, Jaumann/Hughes-Winget, accept start-up overshoot, MC + superbee for history, no nodal solver |
| 1.4 | Single-point constitutive driver: elastic, J2, Drucker-Prager, Mohr-Coulomb | Done | `prototype/constitutive/`; `tests/test_proto_constitutive.sh` passes in CI; PR #2 merged (2026-11-14) |
| 2.1 | Eulerian library skeleton and solver fields | Done | `src/celEulerian` (class `eulerianSolid`), celFoam is a thin driver; `tests/test_eulerian_skeleton.sh` passes serial and on 2 processors (decision 0010). PR #8 |
| 2.2 | Constitutive library interface in C++, elastic and J2 | Done | `src/celConstitutive/` (plain C++ kernel + OpenFOAM wrapper, decision 0006); `tests/test_constitutive_cpp.sh` passes in CI (kernel matches Python to 3e-15; wrapper and celFoam checks pass). PR #4 merged (2026-11-14) |
| 2.3 | Lagrange step in C++ | Done | `eulerianSolidLagrange.C`: MUSCL-Hancock face states, acoustic face Riemann solver, Hughes-Winget, boundaries from the U patch type, default maxCo 0.4 (decision 0018). `tests/test_lagrange_step.sh` passes: wave orders 1.9 to 2.3, rotation 3e-13, shear 3e-7, reflections 1e-4, 2/4 processors 3e-15 (2026-11-14) |
| 2.4 | Remap in C++, single material | Done | `eulerianSolidRemap.C`: unsplit face-based remap, geometry from `phiL` only, Barth-Jespersen form of MC/superbee per field (`remapLimiters` in materialProperties), decision 0017. `tests/test_remap.sh` passes: 1D prototype match 7e-14, order 2.7 (2D) / 2.3 (3D), no new extrema, history edge 2 cells, slotted disc, 2 and 4 processors match serial to 3e-14, conservation 1e-13. Stable to summed Courant about 0.9. Merged on green CI (2026-11-14) |
| 2.4b | Lagrange step plus remap, integration checks | Done | `tests/test_eulerian_cycle.sh` passes (about 3 min): streaming P/S waves along x and at 45 deg, speeds within 0.2 %, orders 1.80 to 2.11; disc one turn 0.11 % (prototype 0.03 %); plastic step on a 3D mesh, speeds and HEL within 0.12 %; 45 deg impact 0.022 % over along x; 2 processors 1e-13. Fixed: prescribed-motion face velocity, remap divide by -0 (decision 0022). Stable to maxCo 0.8 (2D), at least 0.7 (3D); 2.3e-6 s per cell per step |
| 2.5a | Drucker-Prager in C++ | Done | `kernel/celDruckerPrager`, wrapper `druckerPrager` (`fit` word option, decision 0009); 14 new reference paths match Python to 2e-15 with apex returns; celFoam selects it. Merged on green CI (2026-11-14) |
| 2.5b | Mohr-Coulomb in C++ | Done | `kernel/celMohrCoulomb`, Jacobi eigen-solver `kernel/celEigen` (decision 0011); `tests/test_constitutive_cpp.sh` matches Python on 10 MC paths to 1e-14, edge and apex returns required. Tension cut-off split to 2.5c |
| 2.5c | Mohr-Coulomb tension cut-off (Python, then C++) | Done | Optional `sigma_t` / `sigmaT` (Rankine planes, active-set return, decision 0015). Python checks pass (strengths, caps, unloading, 4 x 10,000 tension-biased increments, associated return = independent closest point to 1e-15); C++ matches Python on 10 new paths to 2e-14 with every tension return exercised; celFoam reads `sigmaT`. Merged on green CI (2026-11-14) |
| 3.0 | Silo benchmark definition and reference values | Done | `docs/spec/silo-reference.md` (10 sources opened 2026-11-14; unread values marked not verified), `tools/silo_reference.py`, `tests/test_silo_reference.sh` passes. T0: Janssen K 0.29 to 0.43, bottom stress 80 to 86 percent of hydrostatic; Beverloo C = 1.45 (continuum mu(I) 1.40 to 1.55) gives 65, 120, 185 kg/s/m for D = 0.04, 0.06, 0.08 m. Tolerances: decision 0019 |
| 3.1p | Prototype gravity, walls and settling | Done | `prototype/silo2d.py`; `tests/test_proto_silo_settle.sh` passes (about 3 min): lake at rest 3e-18 m/s, smooth walls 0.13 % from rho g H with K0 = 0.4286 exact, frictional tank (H = 3 W, 10 mm cells) within 1.7 % of Janssen with K = 0.425 measured at the walls, wall ratio within 1.8 % of tan 20 deg, force balance 0.08 %. Gravity by subtracting the discrete weight, walls as one-sided Riemann solutions, ramp 0.1 s plus damping 250 1/s (decision 0020). Mohr-Coulomb vectorised in `constitutive/batch.py` |
| 3.2p | Prototype free surface (material plus void) | Done | `prototype/void2d.py` (`VoidBox`); `tests/test_proto_void.sh` passes (about 70 s): lake at rest with a half-full top cell 2.5e-16 m/s, block translating 50 cells with exact edges and 3.5e-9 Pa stress, free fall at g to 1e-13 and landing, cohesionless block cracks with void = 2 v t. Slab-per-sweep remap, homogenised partly filled cells, void-first strain, separation at pulling faces (decision 0021). Column collapse split to 3.2q |
| 3.2q | Prototype column collapse (split from 3.2p) | Done (run-out); slope split to 3.2r | Cause of the short run-out: the Riemann solver's shear jump term acts as a mesh-independent viscosity in one-cell shear bands (grows with sqrt(E)). Fix: face shear capped at yield, capped faces slip (`shear_cap`, decision 0023). a = 0.5 run-out 0.25 -> 0.40 (10 and 20 cells; band 0.30 to 0.90), top kept; a = 2 1.2 -> 1.6 (band 1.2 to 3.6). `tests/test_proto_void.sh` unchanged and passing |
| 3.2r | Prototype column collapse deposit slope | Done (staircase); toe split to 3.2s | Rule 9 (decision 0024): top layers (al >= 0.1) spread sideways in the x sweep. Deposit body smooth: steepest slope without the toe 27.9 deg (a = 0.5) and 33.9 (a = 2), was 45 and 59; run-outs 0.40 and 1.8 (20 cells 0.45 and 2.45), top kept. The 3.2q blow-up came from spreading in the y sweep, not the cap. Toe still 62.3 / 48.0 deg per cell (10 cells), 46.6 / 39.8 (20 cells); test stays in `tests/manual/` |
| 3.2s | Prototype collapse toe | Done | Decision 0025 (Sia chose the fitted flank angle, 2026-11-15). `tests/test_proto_collapse.sh` now in CI and passing: fitted slope 22.2 / 22.3 deg, run-outs 0.40 / 1.8. Partly filled cells are frictionless, but giving them strength (`FLOOR`) made the toe worse (a = 2: 59.9 deg, run-out 1.5). The a = 0.5 creep is full cells sliding along faces held at yield by the shear cap; making slow faces stick (`VSTICK`) stops it but the toe then reads 30 to 54 deg depending on whether the last column is just over one cell. With creep stopped the toe is still 39.7 deg: the nose of a 1 to 2 cell flowing layer (resolution). Both knobs off, defaults unchanged. Single-cell slopes printed for information only |
| 3.3p | Prototype silo discharge (target T0) | Done (outlet, end-to-end run); rate law split to 3.3q | `prototype/discharge2d.py`, outlet = rule 10 of `void2d.py` (decision 0026). `tests/test_proto_silo_discharge.sh` passes (about 3 to 5 min): nothing leaks, mass 1e-14, head independence 0.975. Results, figures and animation: `docs/results/silo-prototype.md`. Steady rate, stagnant corners and crater are right, but the rate is a third to a half of Beverloo (5 mm, E = 0.1 MPa: C 0.32 to 0.54, exponent 2.26; 5 against 2.5 mm 23 percent apart) and depends on E (6 / 17 / 38 kg/s/m at 10 / 1 / 0.1 MPa): Riemann dissipation |
| 3.3q | Prototype silo discharge rate (split from 3.3p) | Next | Brief: `docs/wp/WP-3.3q.md`. Damping that does not resist plastic flow, and/or mu(I) friction. Test `tests/manual/test_proto_silo_rate.sh` (fails now) |

## Standing instruction from Sia (2026-11-14)

Proceed without asking for permission, and merge each pull request automatically once
CI is green. Sia is away; report results, do not wait for review. Still ask before
anything that cannot be undone outside this repository.

## Waiting on Sia

- Optional: real tank dimensions and material data for the silo problem (placeholders
  are in `docs/spec/target-problems.md`).
- Optional: paste the setup line from `docs/wp/WP-0.4.md` into the cloud environment's
  setup script so every session starts with OpenFOAM installed.
- Confirm GPL-3.0 (decision 0003). Needed before the code is shared with anyone.
- Abaqus or LS-DYNA access for reference runs? Not needed for the silo problem.

## Next

The C++ Eulerian cycle is done and checked against the prototypes: WP 2.3 (Lagrange
step, decision 0018), WP 2.4 (remap, decision 0017), WP 2.4b (both together, decision
0022). The C++ silo packages 3.1 to 3.4 can start once their prototype counterparts
are done (briefs not written yet).

Silo track (decision 0016), one thread each:

1. WP 3.0: reference values. Done (`docs/spec/silo-reference.md`).
2. WP 3.1p: prototype gravity, walls, settling. Done (decision 0020).
3. WP 3.2p: prototype free surface. Done (decision 0021). WP 3.2q fixed the collapse
   run-out (decision 0023), WP 3.2r the deposit staircase (decision 0024); WP 3.2s
   closed the collapse with a fitted flank angle (decision 0025); the collapse test is in CI.
4. WP 3.3p: prototype silo discharge. Done (decision 0026): the tank fills, settles,
   opens and empties (`docs/results/silo-prototype.md`), but the rate is set by the
   Riemann solver's dissipation.
5. WP 3.3q: discharge rate (`docs/wp/WP-3.3q.md`). Start now.

Then the C++ counterparts 3.1 to 3.4 (2.4b is done; briefs not written yet) and 3.5
(cohesion).

After the silo (Phase 4, Lagrangian bodies): WP 4.0 comes first, a thorough
investigation of body-to-material contact schemes (LS-DYNA constraint and penalty
coupling, Abaqus CEL contact, cut-cell, immersed boundary), requested by Sia on
2026-11-15 (`docs/wp/WP-4.0.md`). Its recommendation is shown to Sia before the coupling
is built. Silo walls are mesh boundaries with friction, not contact.

## Known issues

- C++ cycle (WP 2.4b): a 45 degree S step overshoots 0.9 percent after start-up in
  streaming material (test limit 1 percent, thin margin) and 1.3 to 1.4 percent
  without streaming. Cause: the Lagrange step limits primitive components (decision
  0018 item 1); characteristic limiting would fix it if contact or the silo needs it.
- C++ cycle (WP 2.4b): over many crossings a 45 degree P step's overshoot creeps up and
  levels off (0.5 percent first crossing, 2.7 percent after 16, 5.7 percent after 40;
  2D, maxCo 0.4), with or without streaming (decision 0022).
- C++ remap costs about 1.4e-6 s per cell per step, the Lagrange step 8e-7 (2D, serial).
- Void prototype (decisions 0021, 0023): the face shear cap (WP 3.2q) removed most of
  the over-damping of slow granular flow, but the normal jump term of the Riemann
  solver is not capped (capping it was unstable), so run-out still falls with E: 0.45,
  0.40, 0.25 at 0.1, 1, 10 MPa (a = 0.5). Stiff material under-predicts flow. The silo
  discharge rate shows it strongly (WP 3.3p, decision 0026): 6 / 17 / 38 kg/s/m at
  E = 10 / 1 / 0.1 MPa against Beverloo's 120, 23 percent mesh dependence (5 against
  2.5 mm), and a 20 to 40 percent drop from WP 3.2r's free-surface rule. Reducing the jump terms makes the flow noisy; WP 3.3q. Surfaces were staircases of whole cells (fixed in
  WP 3.2r); fronts are still blunt (WP 3.2s). The C++ Lagrange step has the same face solver and
  needs the shear cap before the silo runs in celFoam.
- Void prototype (decisions 0024, 0025): partly filled cells carry no stress (void-first
  closure). The a = 0.5 collapse still creeps at about 1 cm/s at t = 2 s with 10 cells
  (stops at 0.5 s with 20): full cells near the front are at yield in their normal
  stresses (active state), so the shear cap lets their faces slip with no resistance.
  A rate-dependent friction (mu(I)) is the physical cure. The toe is a 1 to 2 cell step
  (the stopped nose of the flowing layer). Expect a creeping free surface in 3.3p.

- WP 3.1p prototype: the free top of a full tank is a zero-traction face with no mass
  flux, so top cells lose density as the column settles (7 % in the top cell for 1.2 m).
  WP 3.2p replaces it. After settling the frictional tank creeps slowly (3e-4 m/s, 1.5 %
  of rho g H over a further 10 % of the run undamped) because the walls sit at the limit.
  A fast ramp (0.05 s) with light damping lets the column rebound and unload the upper
  wall (ratio 10 to 25 % below tan delta_w): use the 0.1 s ramp (decision 0020).

- Silo reference (WP 3.0): every published continuum discharge constant used mu(I)
  with friction 0.32 to 0.38; T0 is Mohr-Coulomb with tan 35 = 0.70. The test band for C
  is therefore wide (1.0 to 1.9); the exponent is the real check. Lube et al. 2005 and
  Lajeunesse et al. 2005 prefactors were not read in the original papers (decision 0019).

- WP 1.1 prototype: smooth-wave order is limiter-bound at about 1.8. It passes at the test
  setting (0.25 c streaming) but gives 1.78 at 0.5 c. Only matters if later packages need
  a sharper order; an extremum-preserving limiter would fix it.
- WP 1.2 prototype: an impulsive start (velocity jump at a face) overshoots by 2.7
  percent of the jump at the impact plane for the first 5 to 10 cells of wave travel
  (1.7 percent elastic). Matters for contact (WP 4.x); revisit in WP 1.3 (decision 0008).
- WP 1.2 prototype: history is remapped with superbee, stress with MC (decision 0008).
- Mohr-Coulomb tension cut-off (WP 2.5c) is optional: without `sigmaT` cohesive
  material still carries tension up to c cot(phi). The vectorised 2D prototype
  (`prototype/constitutive/batch.py`) and Drucker-Prager have no cut-off.
- 2D prototype: a step at 45 degrees overshoots 2.7 percent during its first cell of
  travel, then below 0.3 percent (decision 0013; the 13 percent in decision 0012 was a
  jagged test step). Like the 1D impact start-up, accepted; contact (WP 4.3) re-measures.
- C++ remap (WP 2.4) has no corner transport: stable up to a summed outflow material
  Courant number of about 0.9 (2D 0.45 and 3D 0.30 per direction on a diagonal), far
  above what the acoustic step gives except for near-sonic material speeds. Density
  is not returned to uniform exactly by a compress-then-reverse flux (1e-3 at Courant
  0.2); velocity, stress and history are (decision 0017).
- C++ Lagrange step (WP 2.3): stable Courant number 0.65 in 2D, 0.45 to 0.50 in 3D
  on uniform meshes (default maxCo 0.4); not measured on graded or distorted meshes.
  WP 2.4b found that measure too strict: the full cycle is stable to 0.8 (2D) and at
  least 0.7 (3D); above about 0.6 in 3D the start-up overshoot grows past 5 percent
  (decision 0022). The disc check passed in WP 2.4b.
- `stableDeltaT` uses cbrt(V): 2D one-cell-thick meshes must be as thick as their
  in-plane cells (decision 0018, consequences).
- WP 1.3a test takes about 5 minutes (400 x 400 meshes for the 45 degree convergence runs);
  WP 1.3b test about 1.5 minutes.
