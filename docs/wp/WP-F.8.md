# WP F.8: Collapse in gas (end of Stage B)

State: in progress, split (F.8b, decision 0133)
Branch: claude/project-thread-ndfemq (stacked on F.7)

## Goal

Initially fluidised dam break against Roche et al. (2010); read the paper first (the HAL page refused automated reading in F.0). Fall-back: Rondon et al. (2011).

## Read first

- `docs/decisions/0125-gas-coupling-plan.md`
- `docs/spec/gas-references.md`, `tools/gas_reference.py`
- Briefs of the earlier F packages (their logs)

## Build

- Case, mesh shown to Sia, checkpointed run

## Done when

- `tests/manual/test_fluidised_collapse.sh`: run-out and front speed within 20 % of the paper; dry collapse with gas off unchanged.

Parallel (Sia, 2026-11-24: the solver must be fully parallel): the done-when test also runs on 2 and 4 processors, with the cuts through the bed, the free surface and the outlet. Gas fields match serial to 1e-10 relative, with the pressure solve converged to 1e-12; total gas mass is conserved to 1e-12 on every processor count. Also report the cost per step with gas on against gas off.

Tests use the `ergunWenYu` drag law (Sia, 2026-11-24).

## Out of scope

- Release (F.9).

## Log

Newest first. Each entry: date, what was done, what was learned, exact next step.

- 2026-11-25 (F.8g): diagnosed the slow front from the papers (our carried p acted as a force; the papers use it only to lower friction) and added gasProperties `frictionOnly`. 10 mm: run-out 5.67 h0 (-1 %), front 1.04 sqrt(g h0) (-26 %). Learned: the carried p must also follow the dynamic total pressure (the papers carry the ratio p/rho g h); 2 processors part at the MC apex as the dry frictionless case does. Next: 5 mm run, then carry the ratio.

- 2026-11-25 (F.8f, step 1): added the Princeton packing-dependent pressure (`princetonFr`) as a limit on the effective mean stress, and `tests/test_packing_pressure.sh`. Dam break: dry 2.47 h0, compressible gas 3.50 h0 (-39 %), route B 4.47 h0. Learned: the base pore pressure now drops at release and recovers as in Breard; the flowing part still loses its pressure. Next: Breard's other settings (air through the reservoir, initial packing, wall friction) on the compressible model.

- 2026-11-25 (F.8d, F.8e): dry column calibrated (mu(I) + flowImpedance + fullBlend, Gueugneau 24/19 deg): 2.375 h0, -9 % of Roche's pressure-balanced run. Added gasProperties `diffusivity` (route B: pore pressure carried by the grains and diffusing) and `tests/test_pore_diffusion.sh`. Dam break at kappa 0.015: run-out 4.45 h0 (-23 %), front -25 %; at the theoretical 0.11: 3.30 h0. Learned: the compressible gas's suction at release costs about 1.3 h0 of run-out; with route B we are close but not yet within 20 %. Next: 5 mm mesh, then F.8f.

- 2026-11-25 (F.8c): read Roche (2008) in full and the published models (Gueugneau 2017, Aravena 2021, Breard 2019; decision 0133). Case set to the lab (h0 0.4, a = 2, 15 deg floor, air on). `overfillStress` fixes the remap densification crash. Run-out 1.25 h0. Learned: the dry collapse is short too (0.85 against about 2.6 h0); our compressible gas loses support at release like Breard's, but our elastic-plastic grains then lock. Next: F.8d dry calibration, then F.8e (route B), F.8f (route A), Sia's order.

- 2026-11-25 (F.8b): papers still unreadable. Fixed the blow-ups (gas shut in nearly empty remapped cells: eMinPacked 0.2, face density floor) and two serial/parallel jumps (the cap's 'not set' p_c mixed by the remap: p_c = max(stored, old P, capP0)). Run-out 2.20 h0 (-62 %), front -67 %; parallel 1e-10 to 0.36-0.46 s, 1e-8 at the end; v1 dry check now plain MC (identical to main). Next F.8c: the paper's floor and gas supply; parallel round-off growth.

- 2026-11-24: read the Roche (2008) abstract (full text refused by HAL, AGU, hal.science); case `tests/gas/collapse_cases.py`, mesh shown to Sia, `tests/manual/test_fluidised_collapse.sh` (checkpointed). Found and fixed (option `surfaceCarry`) the free surface sucking gas through the grains; expanded bed made soft with the cap. Result: run-out 1.95 h0 (-66 %), front -70 %: FAIL; dry identical to main, gas mass 1e-14; parallel parts within 0.002 s and blew up; cost 2.3x. Learned: the defluidised base's floor friction holds the flow (no floor friction: 7.0 h0); void flow near a thin front is not robust. Next (F.8b): full paper, then floor/gas supply, robust void flow and drag guard, parallel from an apex start, 5 mm mesh.

- 2026-11-24: brief written in WP F.0. Next: start when the package before it has merged.
- 2026-11-25: 5 mm refit kappa 0.22 (run /mnt/project-files/f8/run10_ratio5mm_k022): run-out 5.64 h0 (-2 %), stops 6.7, front 1.11 sqrt(g h0) (-21.5 %). Theory passes front, refit passes run-out; run-out mesh dependent.
- 2026-11-25 (F.8i): zero pore pressure on the material's surface (tests/test_surface_drainage.sh, 1.1 %). Dam break theoretical kappa: 10 mm 5.17 h0, 5 mm to the channel end. Learned: the surface error was real (15-28 % slow drainage) but the mesh dependence comes from the partly supported frictional flow; fully supported fronts agree within 5 %. Next: dry 5 mm check.
- 2026-11-26 (F.8i): dry check: dry run-out 3.625 h0 on 5 mm against 2.375 h0 on 10 mm (+53 %): the grain law is mesh dependent here.
- 2026-11-26 (F.8j, Sia: fix the grain law's mesh dependence first): dry dam break to 1.4 s, run-out from the gate (h0):

  | Variant | 10 mm | 5 mm |
  | --- | --- | --- |
  | As calibrated (mu(I), flowImpedance, fullBlend) | 2.47 | 3.61 |
  | No flowImpedance | 1.08 | 1.59 |
  | No fullBlend | 1.38 | 1.88 |
  | Plain Mohr-Coulomb 24 deg | 1.08 | 1.56 |
  | Plain MC, E 0.25 / 1 / 4 MPa (10 mm) | 1.40 / 1.08 / 0.93 | |
  | Calibrated, E 0.25 / 1 / 4 MPa (10 mm) | 2.90 / 2.47 / 2.20 | |

  Every variant runs 35-50 % further on 5 mm, plain MC too, and the run-out falls with stiffness: it follows rho c_s D,
  the damping of the face Riemann solver (shear jump term), not the grain law. A rate-independent collapse should not
  depend on E or D. Proposal to Sia pending.
- 2026-11-26 (F.8j, Sia approved "drop shear damping in yielding cells"): `yieldShearJump no` (no elastic shear jump beside yielding cells, no elastic predictor along the flow direction in them). Dry run-out (h0): plain MC 10 mm 1.85 (was 1.08), 5 mm 2.59 (was 1.56): +40 % from 10 to 5 mm, no better; calibrated 2.50 / 3.71 (was 2.47 / 3.61). So it is NOT the cure. Time step: at Courant 0.4 / 0.2 / 0.1 plain MC 10 mm 1.08 / 1.28 / 1.50, calibrated 10 mm 2.47 / - / 3.30 (5 mm at 0.4: 3.61). The loss scales with the time step, and dt halves with the cell; checking 5 mm at Courant 0.2 against 10 mm at 0.1 (same dt).
- 2026-11-26 (F.8j): same time step on both meshes, calibrated dry column (h0): 10 mm at Courant 0.4 / 0.2 / 0.1 / 0.05: 2.47 / 3.03 / 3.30 / 3.25; 5 mm at 0.4 / 0.2: 3.61 / 3.53. Converged in time, 10 and 5 mm differ by about 7 %. The "mesh dependence" is mostly a time-step error that the default Courant 0.4 leaves on 10 mm (dt halves with the cell). 2.5 mm mesh shown to Sia (/mnt/project-files/f8/mesh_2p5mm.png).
- 2026-11-26 (F.8j): 2.5 mm calibrated dry column (Courant 0.4, 4 procs, 1.96 h): 4.43 h0 at 1.4 s, still moving (5 mm 3.61, 10 mm 2.47; 10 mm at Courant 0.1: 3.30). Not converged. The calibrated case ends in a thin creeping film, so its run-out also depends on the height threshold that defines the front (at 1.4 s, front where the deposit is >= 5 mm / >= 10 mm: 10 mm 2.47 / 1.62, 5 mm 3.06 / 2.83, 2.5 mm 4.21 / 2.76). Plain Mohr-Coulomb has a clean front (1.07 / 1.07 on 10 mm, 1.53 / 1.49 on 5 mm), and its 43 % growth from 10 to 5 mm is a real mesh dependence. Back to Sia.
- 2026-11-26 (F.8j, Sia chose "research the core scheme"):
  - Sliding block (scratch script, not in the suite): a block at 2 m/s on a smooth floor without gravity keeps its shape and speed exactly. A cohesive block (c 5 kPa) on a 19 deg floor stops at 0.990 m (exact 0.992). So translation, the void front and wall friction are accurate.
  - Side findings, not pursued: c >= 100 kPa with E 1 MPa blows up within 0.02 s under gravity; the compaction cap combined with cohesion compacts stress-free material (alpha -5 % in free flight, c 1 kPa).
  - Plain Mohr-Coulomb with yieldShearJump no, 10 mm: secondOrder no gives 0.93 / 0.82 h0 at Courant 0.4 / 0.1 (second order 1.85 / 2.35); shearCap + normalCap 1.93 / 2.35 (no effect). The loss is in the face fluxes, not in the caps.
  - Time: 5 mm at Courant 0.4 / 0.2 gives 2.59 / 2.60, converged in time; 10 mm needs Courant <= 0.1.
  - Mesh: front where the deposit is >= 5 mm thick: 10 mm (Courant 0.1) 2.35, 5 mm 2.60, 2.5 mm 2.81 h0 (+10 %, +8 %). With the half-cell threshold, 2.5 mm gives 3.30 (thin sheet).
  - Reading: the remaining dependence is spatial dissipation of the acoustic face fluxes (about rho c h) in slow flow (Mach about 0.1), the known low-Mach weakness of Godunov schemes. Mohr-Coulomb MPM collapses overshoot experiments (Fern & Soga 2016), so the coarse-mesh run-out was short partly from numerical dissipation.
- 2026-11-27 (F.8k, Sia chose "add the slow-flow fix"): `slowFlowJump yes` (decision 0138), test tests/manual/test_collapse_convergence.sh. Plain MC dry run-out at 5 mm front thickness 2.30 / 2.46 / 2.84 h0 on 10 / 5 / 2.5 mm (FAIL, 5 % asked); at 20 mm thickness 1.58 / 1.69 / 1.62 (within 4 %; before 1.40 / 1.67 / 1.75). 10 mm at Courant 0.1: 2.42. A thin sheet runs ahead on fine meshes. Back to Sia.
- 2026-11-27 (F.8k, Sia chose "investigate the thin sheet"): the sheet is the deposit's tapering tail, and each mesh ends it where it is one cell thick. Final thickness (mm) at x = 0.90 / 1.02 / 1.10 / 1.26 / 1.42 m: 10 mm cells 19 / 10 / 9 / 0 / 0; 5 mm 20 / 20 / 12 / 5 / 0; 2.5 mm 20 / 15 / 12.5 / 6 / 2.5. Where all meshes resolve the layer (two coarse cells, 20 mm) the profiles agree within 4 %. Floor normal stress in the front 0.1 m is 0.2-1.2 of the weight on every mesh (stress waves, no trend with the mesh). Early grains thrown ahead (1e-3 alpha at 0.04 s) change between reruns of the same case (round-off) and are too thin to count. Linear extrapolation of the tail to zero thickness: 2.42 / 2.85 / 2.83 h0 (fit 10-30 mm), 2.42 / 2.59 / 2.45 (15-40 mm), not robust. So a grid cannot carry a layer thinner than a cell; the 5 mm-threshold run-out grows with refinement for that reason, not from the flux scheme.
- 2026-11-27 (F.8k, Sia chose "resolved thickness, 5 mm mesh"): convergence test now measures the front at 20 mm thickness: 1.58 / 1.69 / 1.62 h0, PASS (3.8 %). Calibrated grains (mu(I), flowImpedance, fullBlend, cap) + yieldShearJump no + slowFlowJump yes, front at 5 mm thickness:
  - Dry, floor 19 deg (Gueugneau's measured value): 10 mm 3.17, 5 mm 3.27 h0 against 2.6 (+22 / +26 %).
  - Dry sweep, 10 mm: floor 22 / 25 / 28 deg 2.85 / 2.73 / 2.67; phi 28 2.83. Floor 22 on 5 mm: 2.73 h0 (+5 %).
  - Fluidised (route B, frictionOnly, kappa 0.22), floor 22, 5 mm: run-out 4.77 h0 (-17 %, passes 20 %), front speed 2.03 m/s (-28 %, Fr 1.02; fails 20 %); fluidised / dry 1.75 (Roche 1.8-2.2). 10 mm floor 19: 4.97 h0, front -31 %. 5 mm run 2.4 h on 2 processors.
  - collapse_cases.py is not yet changed (floor 22 departs from the measured 19 deg; asked Sia).
