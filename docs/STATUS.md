# Status

Last updated: 2026-11-18 (WP 4.3a; 4.1, 4.0, 4.2; 3.6c2, 3.6c, 3.6b, 3.6s, validation references)

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
| 3.3q | Prototype silo discharge rate (split from 3.3p) | Done (rate law at E = 0.1 MPa); stiffness split to 3.3r | The low rate was an outlet bug, not dissipation: a partly filled outlet cell leaned away from the face, which then got the wall velocity 0 and shut (rule 10b, decision 0027). `tests/test_proto_silo_rate.sh` now in CI and passing (5 mm, four processes, about 20 to 25 min): C 0.99 / 0.99 / 1.06 for D = 0.04 / 0.06 / 0.08 m, exponent 1.59, fitted C 1.035, steady within 6 %, averaged corner speed 0.8 %, head 0.992. Outlet test gained C >= 0.8 (now 0.98). Not met: E (81 / 73 / 54 kg/s/m at 0.1 / 1 / 10 MPa, 10 mm), and 2.5 against 5 mm (93 against 82 kg/s/m, 14 %) |
| 3.3r | Prototype discharge rate at T0 stiffness (split from 3.3q) | Done; collapse with the new rules split to 3.2t | Decision 0028. Causes: the uncapped normal Riemann jump term (numerical viscosity ~ rho c h on plastic flow) and voids opened by one sweep that made interior faces free surfaces. Rule 12 caps the normal jump at yield (per-side strain, volume kept), rule 13 spreads void inside the bulk; both on in the discharge only. `tests/test_proto_silo_stiffness.sh` in CI and passing (about 16 min): 91.2 / 87.1 / 84.0 kg/s/m at E = 0.1 / 1 / 10 MPa (8.7 %, band 10, thin), windows 4 to 5 %, corners 0.1 to 0.6 %. Rate test still passes (exponent 1.48, C 1.133); 5 against 2.5 mm 3 % (was 14). Collapse with the rules runs too far (0.8 / 1.05 / 1.2): WP 3.2t |
| 3.2t | One rule set for collapse and discharge (split from 3.3r) | Closed, done-when not met | Decision 0029: no rule set found in the split prototype. Rule 13 makes stress-free porosity diffuse through the collapsing bulk (cured by keeping isochoric stretch in the material, which slows the discharge, 76 against 87 kg/s/m at 1 MPa). Rule 12 turns the column hydrostatic; capped only at yield the 1 MPa collapse is fine (0.35) but the 10 MPa column never stops (time-step dependent) and the discharge falls to 48 to 52 kg/s/m. Nine more variants tabled. The 0028 E independence leans on capping faces of still-elastic cells. Defaults unchanged, knobs kept off. Stiffness margin not improved |
| 3.1 | C++ gravity and wall boundaries (smooth, frictional) | Done | Decision 0030. `tests/test_silo_settle.sh` passes (57 s): lake at rest 3e-17 m/s; smooth column 0.12 % from rho g H, K0 exact; frictional tank (0.4 x 1.2 m, 10 mm) Janssen within 1.7 % with K 0.424 at the walls (prototype 0.425), wall ratio 0.4 % where it slips, force balance 0.05 %; 2 processors match serial to 7e-18. Gravity by hydrostatic reconstruction per face; walls named in `materialProperties` `walls`; the Hancock predictor uses the wall Riemann solution (without it a one-cell shear layer formed at frictional walls). First step of the C++ rule 12/13 question (decision 0029): needs 3.2 and 3.3 before collapse and discharge run on one code |
| 3.2a | C++ void and free surface (split from 3.2) | Done | Decision 0031: `freeSurface yes` and field `alpha` in eulerianSolid; prototype rules 1 to 7 face by face on Cartesian hex meshes (slab side per axis, no isoAdvector). `tests/test_void.sh` passes in 17 s: lake at rest 2.5e-17 m/s, block translating 50 cells with exact edges, free fall at g, cohesionless block cracks with void 2.000 = 2 v t, 2 processors match serial to 2e-15. Old code paths unchanged bit for bit. Shear cap and collapse are 3.2b |
| 3.2b | C++ column collapse: shear cap and rule 9 | Done (recorded; creep and E dependence remain) | Decision 0032: `shearCap yes` (rule 8, face by face, capped faces slip) and `spreadLayers yes` (rule 9) in materialProperties, both default off. `tests/test_collapse.sh` passes (about 2.5 min, 10 mm): a = 0.5 run-out 0.55 / 0.35 / 0.20 at E = 0.1 / 1 / 10 MPa (prototype 0.55 / 0.40 / 0.25), top kept, fitted slope 17 to 24 deg; a = 2 at 1 MPa 1.40 (prototype 1.8), stops at 1.5 s. The a = 0.5 columns still creep at 2 s (as the prototype at 10 cells). 10 MPa run-out below the 0019 band (normal jump term uncapped, the rule 12 question for WP 3.4). Without the cap 0.25 at 1 MPa. 2 processors match serial to 8e-13 at 0.1 s; later round-off flips cap/rule 9 switches and the runs part |
| 3.3 | C++ outlet boundary that opens | Done | Decision 0033: `outlets { outlet { friction 20; openTime 0.3; } }` in materialProperties (needs `freeSurface yes`): a closed wall until it opens, then a traction-free outflow face with rule 10b; damping off on opening; `discharged()` counts the mass out. `tests/test_outlet.sh` passes (about 30 s, 10 mm, E = 0.1 MPa, shear cap and rule 9 on): nothing through the closed outlet (0 exactly), tank plus discharged constant to 4e-14, rate 74.1 kg/s/m, C 0.89 (prototype without rules 12/13 81.35, -9 %); without rule 10b 27.4. Not tested in parallel (WP 3.4) |
| 3.4 | C++ silo discharge across E, serial and parallel | Closed, done-when not met; E independence split to 3.4b | Decision 0034. `tests/test_discharge.sh` passes (about 50 s): mass 4e-14, no stress in void, 2 processors with the outlet split give the serial rate (74.09 both). Rate with the collapse code falls with E: 10 mm 74.1 / 65.7 / 70.9 kg/s/m at 0.1 / 1 / 10 MPa (1 MPa -11 %), 5 mm 78.0 / 55.0 at 0.1 / 1 MPa (-29 %). The E sweep is `tests/manual/test_discharge_stiffness.sh` (fails now) |
| 3.4b | C++ discharge rate independent of E (rule 12) | Closed, done-when not met; Sia chose mu(I) (decision 0036) | Decision 0035: rule 12 in C++ behind `normalCap yes` (default off). With it: 10 mm 78.6 / 73.8 / 84.8 kg/s/m at 0.1 / 1 / 10 MPa (met), 5 mm 82.6 / 69.7 at 0.1 / 1 MPa (-16 %, was -29 %: not met). Collapse unharmed, unlike the prototype (a = 0.5: 0.50 / 0.35 / 0.25; a = 2: 1.50). 2-processor discharge rate +0.6 %; collapse cell-by-cell parallel check fails with the switch on (bisection round-off). Drivers read `CEL_EXTRA_MATERIAL` |
| 3.6a | mu(I) friction in Drucker-Prager and Mohr-Coulomb (Python, C++) | Done | Decision 0037. Optional `mu2`, `I0`, `d`, `rhoS` (and `muIAlpha`): the return uses phi_m = atan(mu(I)), I from the plastic shear rate and the pressure after the return, solved per point (Illinois). `tests/test_mu_i.sh` passes (about 55 s): steady simple shear at fixed pressure gives mu(I) to 2e-14 for I = 1e-4 to 1; mu2 = tan(phi) bit for bit (Python and C++); C++ matches Python on 30 paths to 1e-12; celFoam passes dt and prints the law. Barker-Gray regularisation impossible for mu_s = tan 35 deg (no well-posed range), so T0 uses the plain Jop law. Placeholders: mu_2 - mu_s, I0 from glass beads |
| 3.6s | Mesh dependence of the collapse shear band (split from 3.6b) | Done (finding recorded) | Decision 0038, `docs/wp/WP-3.6s.md`. a = 0.5, E = 1 MPa on 10 / 5 / 2.5 mm: band 1 to 4 cells wide on every mesh, peak shear rate 10 / 23 / 30 1/s, run-out (5 mm front) 0.30 / 0.33 / 0.43 at 1 s: **mesh dependent, not converged**; deposit body and top converge. mu(I) with the placeholders changes nothing at 10 and 5 mm (I about 0.05). Cure deferred (compressible mu(I)/phi(I) preferred). Fixed a floating-point trap in `inertialNumber` at P about 1e-320 Pa. `tests/manual/test_shear_band.sh` (about 9 min) |
| 3.6b | Steady inclined layer (Bagnold) with mu(I), E = 0.1, 1, 10 MPa | Closed, done-when not met; measurement recorded | Decision 0039. Layer 0.05 m on 31 deg, x periodic, tilted gravity, from the exact steady state, 4 s. The body flows at 0.14 to 0.20 of the exact speed at 0.1 MPa, 0.05 to 0.09 at 1 MPa, 0.01 to 0.02 at 10 MPa (5 mm; 2.5 mm doubles it): numerical viscosity about 0.06 rho c_s h (2.8 / 6.5 / 29 Pa s) against a mu(I) tangent viscosity of 0.8 Pa s, so mu(I) alone cannot remove the E dependence. Fixed: the hydrostatic face reconstruction for gravity off the mesh axes (the layer slid uphill; axis-aligned runs unchanged bit for bit). The top cell slides over the body at 1 and 10 MPa. `tests/test_inclined_layer.sh` (CI, 10 s: at rest at 25 deg, 2.4e-5 m/s); done-when `tests/manual/test_inclined_layer.sh` (3 min, fails now) |
| 3.6c | Flow impedance at flowing faces | Closed, done-when not met at 5 mm; rest split to 3.6c2 | Decision 0040: `flowImpedance yes` (MC with mu(I), default off). Three sources of resistance, not one: the elastic Hancock predictor (G dt gamma_dot / 2 on the faces), the elastic shear jump, and the static-phi shear cap (flowing faces slip). Now: mu(I) tangent viscosity in the shear jump between flowing cells and at walls, plastic predictor along the flow direction, cap skipped between flowing cells and at phi_m elsewhere. Layer body 0.65 to 0.97 of exact at 0.1 MPa (was 0.14 to 0.20), 0.41 to 0.92 at 1 MPa, 0.03 to 0.74 at 10 MPa (5 mm); 2.5 mm within 5 % at 0.1 and 1 MPa. Not met: 5 mm (floor cells), 10 MPa (floor cell stuck), top cell slides at 1 and 10 MPa. CI guard in `tests/test_inclined_layer.sh` |
| 3.6c2 | Stiff inclined layer: floor cell and top cell (split from 3.6c) | Done | Decision 0044, both behind `flowImpedance yes`: the wall reconstruction of a flowing cell mirrors about the wall's last velocity (the wall jump counted the drag twice; the 10 MPa floor cell stuck), and void below 1e-4 counts as full (noise void switched the top cell to the partly filled rules, which kicked its stress). `tests/manual/test_inclined_layer.sh` passes: profile error 3.6 / 0.9 / 0.8 % at 0.1 / 1 / 10 MPa, 5 mm (was 7.8 / 22 / 136). CI runs 0.1 and 10 MPa. Left: the 0.1 MPa layer creeps faster over time (+0.7 %/s) |
| 4.0 | Contact scheme investigation (LS-DYNA, Abaqus, cut cell) | Done; Sia agreed 2026-11-17 | Decision 0042: cut-cell moving walls (Riemann wall face with the body velocity, Coulomb cap, small-cell mixing), no penalty. `docs/research/contact-schemes.md` (9 sources), `prototype/contact2d.py`; `tests/test_proto_contact_schemes.sh` passes (about 25 s): piston stress 1.006 zv, noise 0.2 % (penalty at k = M/h: 27 %), 30 deg frictional wall acceleration within 0.3 % (penalty 11 % off), mass 1e-15. LS-DYNA-style coupling points tried after Sia's note: smoother than a plain spring but leak 2 % and half the force; decision 0042 addendum adopts its point density, leak report and controls into the cut-cell wall |
| 4.1 | Rigid body: STL input and kinematics | Done | Decision 0041: STL walls and bodies are immersed rigid surfaces on the Cartesian mesh (Sia, 2026-11-17). `src/celLagrangian` (`rigidBody`, `rigidBodies`), `constant/rigidBodies` (tutorial `rigidBodyMotion`): Function1 velocity and omega, exact path; no forces yet. `tests/test_rigid_body.sh` passes in about 1 s: pose 4e-16, STL 3e-8, 2 processors identical |
| 4.2 | Cut cells near each body | Done | Decision 0043: exact edge, face and cell coverage from the STL triangles (no sampling), body-frame octree, utility `celBodyCells`. `tests/test_body_cells.sh` passes in about 6 s: faceted sphere volume to 4e-16, plane at 30 deg to 2e-7 (single-precision STL), moving and turning cube exact, 2 processors 9e-16. Search 6e-5 s per searched cell |
| 4.3a | Fixed STL walls with friction | Done | Decision 0047: cut-cell walls in `eulerianSolid` (open areas and volumes, wall face with the body's friction, ghosts in covered cells, small cells one cell with a neighbour in the Lagrange step), `friction <deg>` per body in `constant/rigidBodies`, per-body wall force and leak report; celFoam couples fixed bodies (tutorial `stlTank`). `tests/test_stl_walls.sh` passes (about 3 min): STL side walls 0.24 % from the mesh-wall Janssen tank, tank turned 10 deg 2.1 %, leak 1e-11, force balance 0.06 %. Not with `freeSurface` (4.3c); no mixing across processors |
| 4.3c | STL walls with the free surface and outlets (split from 4.3a) | Next | Brief `docs/wp/WP-4.3c.md`; needed by 4.5 |

## Validation references

`docs/spec/validation-references.md` (2026-11-17): six papers with field data (silo
velocity and pressure fields, collapse profiles, inclined-plane mu(I) data) for later
comparison packages, after WP 3.6b to 3.6e. Some DOIs still to resolve.

## Standing instruction from Sia (2026-11-14)

Proceed without asking for permission, and merge each pull request automatically once
CI is green. Sia is away; report results, do not wait for review. Still ask before
anything that cannot be undone outside this repository.

## Waiting on Sia

- mu(I) constants (mu_2, I0) and grain size and density for the T0 material, if
  measured (decisions 0036, 0037). Placeholders now: mu_s = tan 35 deg, mu_2 - mu_s =
  0.26 and I0 = 0.279 from Jop et al. (2005) glass beads, d = 1 mm, rho_s = 2650 kg/m^3.
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
   opens and empties (`docs/results/silo-prototype.md`).
5. WP 3.3q: discharge rate. Done (decision 0027): an outlet bug held the rate at a third
   of Beverloo; fixed, the rate follows D^1.5 with C about 1.0 at E = 0.1 MPa.
6. WP 3.3r: the rate at T0's stiffness. Done (decision 0028): rules 12 and 13 make the
   rate independent of E (within 9 percent from 0.1 to 10 MPa) and of the mesh (3
   percent between 5 and 2.5 mm).
7. WP 3.2t: one rule set for the collapse and the discharge. Closed without one
   (decision 0029): in the split prototype every cap that lets the 10 MPa collapse stop
   also slows the discharge.
8. WP 3.1: C++ gravity and walls. Done (decision 0030).
9. WP 3.2a: C++ void and free surface. Done (decision 0031).
10. WP 3.2b: C++ column collapse. Done (decision 0032): shear cap and rule 9 behind
   `shearCap` / `spreadLayers`; run-outs 0.55 / 0.35 / 0.20 at 0.1 / 1 / 10 MPa, close to
   the prototype.
11. WP 3.3: C++ outlet. Done (decision 0033): closed while settling, then open (rule
   10b); 74.1 kg/s/m at E = 0.1 MPa, 9 percent below the prototype.
12. WP 3.4: C++ discharge across E. Closed without its done-when (decision 0034): with
   the collapse code the rate falls with E (5 mm: 78 at 0.1 MPa, 55 at 1 MPa; 10 mm
   hides most of it); two processors give the serial rate.
13. WP 3.4b: rule 12 in C++. Closed without its done-when (decision 0035): the collapse
   is fine with it and 10 mm passes, but at 5 mm the rate still falls 16 % from 0.1 to
   1 MPa.
14. Sia chose rate-dependent friction over a soft E (2026-11-16: parameters must be
   measurable). Survey `docs/rate-dependent-friction-survey.md`, plan decision 0036:
   mu(I) with the Barker-Gray regularisation in the return mapping, plus an impedance
   from the plastic tangent at flowing faces (mu(I) alone is estimated to be ~100 times
   too weak against the numerical viscosity at 5 mm).
15. WP 3.6a: mu(I) in the constitutive library. Done (decision 0037). The Barker-Gray
   regularisation does not exist for mu_s = tan 35 deg; T0 uses the plain Jop law.
16. WP 3.6s: shear bands. Done (decision 0038): the collapse band is set by the cell
   (1 to 4 cells at 10, 5, 2.5 mm) and run-out does not converge; mu(I) does not change
   that. The cure (compressible mu(I)/phi(I), or non-local) is a later package.
17. WP 3.6b: inclined layer. Closed without its done-when (decision 0039): the
   numerical viscosity is about 0.06 rho c_s h, 3 to 35 times the mu(I) tangent
   viscosity from 0.1 to 10 MPa; the flow is 5 to 100 times too slow.
18. WP 3.6c: flow impedance. Closed without its done-when (decision 0040): the layer
   now flows at 0.65 to 0.97 of exact at 0.1 MPa (was 0.14 to 0.20) and within 5 % at
   2.5 mm for 0.1 and 1 MPa; 10 MPa and the top cell remain.
19. WP 3.6c2: stiff layer. Done (decision 0044): wall ghost for flowing cells and a 1e-4
   void tolerance; the layer is within 5 % of Bagnold at 0.1, 1 and 10 MPa (5 mm).
20. Next: WP 3.6d (silo E sweep with `flowImpedance yes`, on two meshes per decision
   0038; the first mu(I) silo case to look at), then 3.6e. WP 3.5 waits.
   Keep runs coarse and short: Sia asked to be frugal with Claude usage (2026-11-15).

Then the C++ counterparts 3.1 to 3.4 (2.4b is done; briefs not written yet) and 3.5
(cohesion).

Phase 4 (Lagrangian bodies and STL walls). Sia asked on 2026-11-17 for wall
interactions defined by STL (decision 0041: immersed rigid surfaces on the Cartesian
mesh, as Abaqus CEL). Order:

1. WP 4.1: STL input and kinematics. Done.
2. WP 4.0: contact scheme investigation. Done (decision 0042: cut-cell moving walls,
   not penalty). Shown to Sia 2026-11-17.
3. WP 4.2: cut cells near each body. Done (decision 0043).
4. WP 4.3a: fixed STL walls with friction. Done (decision 0047; tutorial `stlTank`).
5. WP 4.3c: STL walls with the free surface and outlets (`docs/wp/WP-4.3c.md`). Next.
6. WP 4.5: silo and hopper with STL walls (`docs/wp/WP-4.5.md`).
7. WP 4.3b: moving bodies, forces both ways; then 4.4 benchmarks.

## Known issues

- Meshes (decision 0041): `freeSurface yes` (the silo physics) needs a Cartesian hex
  mesh (every face normal to an axis; grading allowed). Without void the code runs on
  any mesh but is untested there: face states assume the face midway between centres,
  and the time step uses cbrt(V). General meshes are WP 3B.4, not planned.

- Shear bands (WP 3.6s, decision 0038): the C++ column collapse localises into bands 1 to
  4 cells wide on every mesh; run-out grows with refinement (0.30 / 0.33 / 0.43 at 10 / 5 /
  2.5 mm, a = 0.5, E = 1 MPa). mu(I) as in T0 does not regularise it. Collapse run-outs
  are mesh-bound until a regularisation package is done.

- mu(I) (WP 3.6a, decision 0037): for mu_s = tan 35 deg the Jop law is ill posed (in
  Barker and Gray's incompressible sense) at every I, and their regularisation cannot be
  built. Whether our compressible elastic model is affected is not known (WP 3.6b).
  The cost of the per-point mu(I) solve in celFoam is not measured.
- Inclined layer (WP 3.6b, decision 0039): the face Riemann solver's numerical
  viscosity (about 0.06 rho c_s h) swamps mu(I) even at 0.1 MPa; the steady layer runs
  at 0.14 to 0.20 of the exact speed (5 mm). The top cell of a layer under tilted
  gravity slides over the body (2.6 times the exact speed at 10 MPa); likely the shear
  cap at the static phi. WP 3.6c (decision 0040) found three sources (elastic Hancock
  predictor, elastic shear jump, static-phi cap) and fixed most of it behind
  `flowImpedance yes`; WP 3.6c2 (decision 0044) fixed the 10 MPa floor cell and the
  sliding top cell. Left: at 0.1 MPa the layer creeps faster over time (+0.7 % of the
  surface speed per second at 5 mm, +1.05 % at 2.5 mm; cells 3 to 11 % fast at 4 s).
  `flowImpedance` counts void below 1e-4 as full and cells within 3 % of the static
  surface as flowing: watch static regions in the silo (WP 3.6d).

- C++ silo (decisions 0034, 0035): with the collapse code the discharge rate falls with E
  (-29 % from 0.1 to 1 MPa at 5 mm; -16 % with `normalCap yes`). Fix planned as WP 3.6
  (decision 0036); Sia rejected a soft E as unmeasurable.

- C++ silo (WP 3.1, decision 0030): the free top of a full tank lets material in through
  the remap as the column settles (mass +6e-4); WP 3.2 replaces it with void. Velocity,
  symmetry and traction boundaries still give the Hancock predictor the cell value
  (walls use their Riemann solution).
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
  discharge numbers of WP 3.3p (decision 0026) were mostly an outlet bug (decision
  0027); with it fixed the rate still fell with E, fixed for the discharge in WP 3.3r by
  capping the normal jump term too (rule 12, with per-side strain; decision 0028) and
  dispersing interior void (rule 13). With both rules the collapse runs too far (a = 0.5
  run-out 0.8 / 1.05 / 1.2 at 0.1 / 1 / 10 MPa), so they are off there; WP 3.2t found
  no variant that serves both (decision 0029), and the discharge's E independence leans
  partly on capping faces of cells that are still elastic (softening near the outlet). Surfaces were staircases of whole cells (fixed in
  WP 3.2r); fronts are still blunt (WP 3.2s). The C++ Lagrange step has the same face solver and
  needs the shear cap before the silo runs in celFoam.
- Void prototype (decisions 0024, 0025): partly filled cells carry no stress (void-first
  closure). The a = 0.5 collapse still creeps at about 1 cm/s at t = 2 s with 10 cells
  (stops at 0.5 s with 20): full cells near the front are at yield in their normal
  stresses (active state), so the shear cap lets their faces slip with no resistance.
  A rate-dependent friction (mu(I)) is the physical cure. The toe is a 1 to 2 cell step
  (the stopped nose of the flowing layer).

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
