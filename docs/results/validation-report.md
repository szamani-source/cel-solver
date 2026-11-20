# celFoam validation report (WP 6.3b)

Date: 2026-11-19 (updated the same day by the accuracy pass, decisions 0090 and 0091). Covers every benchmark run up to the end of phase 5 and the phase 6
checks. Each line names the test that reproduces it; "CI" tests run on every push,
"manual" ones are in `tests/manual/` (too long for CI, or not passing).

Reference data come from published sources (`docs/spec/silo-reference.md`,
`docs/spec/validation-references.md`, the decision notes). Material values are
placeholders unless stated (Sia, 2026-11-18). Sia accepted about 20 % uncertainty on
discharge rates for now; accuracy work comes after the solver is complete.

Material T0: cohesionless Mohr-Coulomb, rho 1800 kg/m^3, phi 35 deg, nu 0.3, E 0.1 to
10 MPa, with the T0 rule set (mu(I) friction, flow impedance, normal cap, decision 0050).
Cells are 10 mm unless stated.

## Summary

| # | Benchmark | Reference | Result | Tolerance | Status | Test |
|---|---|---|---|---|---|---|
| 1 | Silo discharge, T0 tank, rate across stiffness | Beverloo, C 1.0 to 1.9 for mu_s 0.32 to 0.38 (decision 0019) | 65.8 / 52.0 / 60.8 kg/s/m at E 0.1 / 1 / 10 MPa (C 0.79 / 0.63 / 0.73); 5 mm 77.6 / 61.2 | E-independent within 10 % | Known issue: about 20 % spread | manual `test_discharge_stiffness.sh`, `test_discharge_scatter.sh`; CI `test_discharge.sh`, `test_outlet.sh` |
| 2 | Settled tank, wall stresses | Janssen, K 0.29 to 0.43; Jaky K0 | within 1.7 %, K 0.424 at the walls; K0 exact | 2 % | Pass | CI `test_silo_settle.sh` |
| 3 | Steady inclined layer (31 deg) | Bagnold profile, Jop mu(I) | profile error 3.6 / 0.9 / 0.8 % at E 0.1 / 1 / 10 MPa (5 mm) | 5 % | Pass | CI `test_inclined_layer.sh`, manual `test_inclined_layer.sh` |
| 4 | Sand column collapse, a = 0.5 and 2 | Lube et al. 2005, Lajeunesse et al. 2005 run-out bands | run-out 0.40 / 1.70 (1 MPa), 0.30 / 1.30 (10 MPa) | bands 0.30 to 0.90, 1.2 to 3.6 | Pass (bottom of the bands; mesh dependent) | manual `test_collapse_mui.sh`; CI `test_collapse.sh` |
| 5 | Cohesion: arching above the outlet | Jenike critical width, 109 mm (c 500 Pa), 55 mm (c 250 Pa) | stops at 80 / 100 mm, flows at 120 / 140 mm; c 250 Pa switches between 40 and 60 mm | bracket | Pass | CI `test_cohesion.sh` |
| 6 | Tank with STL walls and floor | same tank with mesh walls | 68.8 against 65.8 kg/s/m (+4.4 %); hopper 60 deg 66.2 (Rose-Tanaka predicts 1.21 x flat, got 1.01; 1.07 at 5 mm, decision 0091) | 5 % | Pass (close to the limit) | CI `test_stl_silo.sh` |
| 7 | Strip footing on Tresca clay | Prandtl, Nc = 2 + pi = 5.14 | 4.73 (-8.1 %) at 10 mm, 5.04 (-1.9 %) at 5 mm | 10 % | Pass | CI `test_footing.sh` |
| 8 | Strip footing on sand (phi 30) | Martin 2005, N_gamma 7.653 | psi 30: 20 at 0.01 m/s, 11.9 at 0.00125 m/s (10 mm); psi 0 about 11 at 10 mm, 9.7 at 5 mm | 15 % | Known issue: numerical rate effect, rho c v (decision 0095) | manual `test_sand_footing.sh`; CI `test_body_above_surface.sh`, `test_sand_dilation.sh` |
| 9 | Pile section through clay, smooth | Randolph and Houlsby 1984, 9.14 | 9.29 (+1.7 %) | 10 % | Pass | CI `test_pile.sh` |
| 10 | Pile section, rough (adhesion) | Randolph and Houlsby, 11.94 | 10.47 (-12 %) at 10 mm, 11.26 (-5.7 %) at 5 mm | 15 % | Pass | CI `test_rough_pile.sh` |
| 11 | Taylor bar, Lagrangian hex8 (3D quarter) | Abaqus Benchmarks copper rod: shortening 13.10 mm, widening 5.528 mm (C3D8R) | 13.14 (+0.3 %), 5.76 (+4.2 %) | 5 % / 10 % | Pass | CI `test_taylor_bar.sh` |
| 12 | Taylor bar, Eulerian on the axisymmetric wedge | Abaqus axisymmetric CAX4R: 13.11 mm, 6.006 mm | 13.02 (-0.7 %), 5.02 (-16 %) at 0.2 mm cells | 5 % / 20 % | Pass; widening short (known issue) | CI `test_taylor_wedge.sh` |
| 13 | Silo fed from the top while it discharges | unfed rate; mass balance | outflow -1.0 / -2.7 % with feed 36.8 / 73.7 kg/s/m; matched feed holds the tank within 1.05 %; balance 4e-14 | 5 % | Pass | CI `test_fed_silo.sh` |
| 14 | Cylindrical silo on the wedge against a 3D quarter model | 3D Beverloo; each other | C 0.678 (wedge) and 0.680 (quarter, square outlet by its hydraulic diameter), 0.2 % apart (6.8 % by an equal-area circle, decision 0090); wedge 5 times faster | 10 % / 3 % | Pass | CI `test_wedge_silo.sh` |
| 15 | Projectile: free steel slug into clay | same slug as a rigid body; Prandtl energy estimate 32.3 mm | 27.60 mm against rigid 27.52 (+0.3 %); both 15 % short of the estimate (26 % before A.3b, decision 0094) | 5 % / 30 % | Pass; 15 % short (known issue) | CI `test_projectile.sh` |
| 16 | Deformable footing push (steel hex8) | the same footing rigid | mean force -3 %, density 1.013 rho0 | 5 % / 2 % | Pass | CI `test_deformable_push.sh` |
| 17 | Free body dropped onto T0 | its weight at rest; momentum | force / weight 1.0000, momentum 1e-5 | 2 % | Pass | CI `test_free_body.sh` |
| 18 | MPI with bodies (1, 2, 4 processors) | serial run | first 100 steps 5e-15; outcomes within the round-off band | decision 0080 | Pass | CI `test_mpi_bodies.sh` |
| 19 | Restart | continuous run | identical bit for bit (serial and 2 processors, with bodies) | exact | Pass | CI `test_restart.sh` |

## Verification checks underneath

These check the numerics against exact or independent answers and run in CI:

- Constitutive library (elastic, J2, Drucker-Prager, Mohr-Coulomb with tension cut-off,
  mu(I), dilatancy): C++ matches the Python reference to 1e-12 to 1e-15 on hundreds of
  strain paths; steady shear reaches mu(I) to 2e-14 (`test_constitutive_cpp.sh`,
  `test_mu_i.sh`, `test_dilatancy.sh`).
- Eulerian cycle: elastic P and S waves at 0 and 45 deg, order 1.8 to 2.1, speeds within
  0.2 %; rigid rotation 0.11 % per turn; plastic wave and elastic limit within 0.12 %
  (`test_lagrange_step.sh`, `test_remap.sh`, `test_eulerian_cycle.sh`).
- Free surface: lake at rest 2.5e-17 m/s, block translating 50 cells with exact edges,
  free fall at g, opening gap 2 v t exact (`test_void.sh`).
- Axisymmetric wedge: Lame cylinder 0.9 % beyond two wall cells, cylindrical Janssen
  2.0 % (`test_wedge.sh`, `test_wedge_void.sh`).
- Hex8 element: patch test 1e-15, bar wave 0.3 %, energy 1.4e-7 (`test_hex8.sh`).
- Cut cells from STL: volumes exact to 4e-16 (`test_body_cells.sh`, `test_stl_walls.sh`).

## Known issues

1. **Discharge rate depends on the artificial stiffness (about 20 %).** At E = 1 MPa the
   T0 tank discharges about 20 % slower than at 0.1 and 10 MPa (bands of three runs; single
   10 MPa runs scatter up to 18 % from round-off). The cause is how loosened cells above the
   outlet count (full or partly filled, decision 0046); dilatancy (0048) and a blend of the
   rules (0049) did not cure it. C is also below the published band, partly because T0's
   friction (tan 35 deg = 0.70) is about twice that of the beads in the references.
   Accepted by Sia for now; mass flow is secondary until Sia's validation case (2026-11-19).
   A small dilation angle narrows the gap: 63.0 / 55.7 at psi 0.5 deg, 61.8 / 52.5 at
   psi 1 deg (0.1 / 1 MPa, decision 0091).
2. **Sand footing 27 % high or more** (decisions 0057, 0091). psi = 0 reads N_gamma
   about 9.7 at 5 mm against 7.653 (and a psi 0 answer should be lower still, about 4.4
   by the Davis reduction). Cause found in part: a rigid body within a cell of the free
   surface disturbs the material at rest (edge cells at 0.1 to 0.2 m/s and 20 kPa for
   psi 0; associated flow blows up before contact; cause: small-cell mixing into empty
   slivers, decision 0092). Fixing that leaves psi 0 unchanged, so the excess has
   another cause too. The hover effect is fixed (A.3b, decision 0094,
   `tests/test_body_above_surface.sh`). A.4 (decision 0095) found the rest: the reading
   scales with the push speed times rho c (numerical viscosity of the face Riemann
   solver and predictor, as decision 0040); psi 30 now runs (`voidPlasticFlow`) and
   reads 20 at 0.01 m/s, 11.9 at 0.00125 m/s, still falling. Treat bearing results on
   sand as upper estimates until A.5.
3. **Projectile 15 % short** (decisions 0074, 0094). The slug stops short of the
   quasi-static Prandtl estimate. It read 26 % short (more so at low speed: -33 % at
   0.5 m/s, -21 % at 2 m/s) until A.3b removed the early contact from small-cell
   mixing. Related to the clay footing reading 8 % low at 10 mm.
4. **Wedge Taylor bar widening short** (decision 0063). 5.02 against 6.006 mm (-16 %) at
   0.2 mm cells: the thin lip at the impact face is smeared over a cell or two; first
   order in the cell size. Use the Lagrangian body for thin flanges and lips.
5. **Shear bands are mesh dependent** (decision 0038): collapse run-outs grow with
   refinement (0.30 / 0.33 / 0.43 at 10 / 5 / 2.5 mm). mu(I) does not regularise them.
6. **STL tank +4.4 %** against mesh walls (decision 0052), within its 5 % limit. The
   60 deg hopper speeds the flow by 1.01 at 10 mm and 1.07 at 5 mm against Rose-Tanaka's
   1.21 (decision 0091): it is resolution, so silo results need 5 mm cells or finer
   (no larger than half the shear band of 10 to 20 grain diameters).
7. **Soft deformable bodies** need hourglass 0.3 and loads spread over the wall face
   (decision 0073); a body finer than the cells packs the material.
8. **Creep**: a stopped cohesive tank creeps out at under 1 % of the Beverloo rate
   (decision 0054); the inclined layer at 0.1 MPa speeds up 0.7 % per second.
9. **Parallel runs** agree with serial only to the round-off band once a switch flips
   (decision 0080); compare bands of runs.

## Cost

About 4 microseconds per cell per step on one core with the T0 rule set; 4 processors
give 2.6 (10 mm, 620 cells each) to 3.3 (5 mm) times serial
(`docs/results/scaling.md`). The 10 mm T0 tank discharge (2.3 s, E 10 MPa) takes about
11 minutes on one core.

## What is not validated

- Real material data: mu(I) constants, grain size and the tank are placeholders.
- 3D silos beyond the quarter model of item 14; non-central outlets.
- Rate and temperature effects in metals (no Johnson-Cook yet).
- Field data (velocity and pressure fields in silos, collapse profiles over time) from
  `docs/spec/validation-references.md`: kept for the accuracy work.
