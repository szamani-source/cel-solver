# 0032 Shear cap and rule 9 in C++, column collapse at three stiffnesses (WP 3.2b)

Date: 2026-11-15. Status: accepted (WP 3.2b). Follows decisions 0023 (rule 8, face shear
cap) and 0024 (rule 9, top layers spread sideways) in the unsplit C++ step of decision
0031; neither rule 12 nor 13 (decision 0029).

## Context

The C++ cycle has void since WP 3.2a (decision 0031) but the same acoustic face solver
as the prototype, whose shear jump term over-damps slow plastic flow in one-cell shear
bands (decision 0023). Decision 0029 moved the rule 12/13 question to C++: run the
collapse and the discharge with one code, and add a cap only if the C++ discharge rate
depends on E. This package is the collapse half.

## Decision

1. **Switches** in `materialProperties`, both default `no` so every earlier case runs
   bit for bit: `shearCap yes;` (rule 8, needs `phi`, reads `c`, default 0) and
   `spreadLayers yes;` (rule 9, needs `freeSurface yes`).
2. **Shear cap, face by face** (`capShear` in `eulerianSolidLagrange.C`). At each
   contact face (both sides touch and the Riemann solution is kept; with
   `freeSurface no` every interior face) the face stress is the mean of the two cells'
   stresses (alpha sigma with the free surface) with its normal row replaced by the
   Riemann traction. The shear jump term, a vector in 3D, is scaled by the largest
   fraction (20 bisection steps) that keeps the Mohr-Coulomb function (principal
   stresses from `eigenValues`) no larger than without the jump, or zero. In 2D this is
   the prototype's face stress (normal and shear traction, mean tt and zz).
3. **Capped faces slip.** Each side sees the velocity of its own characteristic at the
   capped traction, v_tL + (t - t_tL)/z_sL and v_tR - (t - t_tR)/z_sR, in its velocity
   gradient; momentum uses the one capped traction. The two slip velocities are face
   fields beside vStar; on coupled patches the owner side solves and they swap sides
   with the face values (`syncCoupled`). The normal jump term is not capped.
4. **Rule 9 per axis.** Gravity picks the vertical axis (largest component of g). A
   partly filled cell with alpha >= 0.1 is spread (lean 0) along each other axis where
   its Youngs gradient (centred alpha differences summed 1-2-1 over the neighbours
   across) is larger along gravity than along that axis. Neighbours beyond a
   non-coupled boundary are the cell itself, as the prototype pads by copying. The new
   lean goes to the contact test and the remap, as in the prototype.
5. **Test** `tests/test_collapse.sh` (`tests/collapse/collapse_checks.py`, driver
   `celVoidCheck` with a variable time step at maxCo 0.4 and snapshots): the
   prototype's set-up and measures (decision 0025), 10 mm cells, t_max 2 s. The 1 MPa
   runs are checked against decision 0019 (as the prototype's CI test); 0.1 and 10 MPa
   are recorded, not checked: their run-out falling with E is the known effect of the
   uncapped normal jump term, which is the rule 12 question for WP 3.4.

## Result (`tests/test_collapse.sh`, 10 mm cells, about 2.5 min)

| Case | C++ run-out | Prototype (0029, 0025) | Band | Stops? | Fitted slope |
| --- | --- | --- | --- | --- | --- |
| a = 0.5, E = 0.1 MPa | 0.55 (0.35 at 0.5 s, 0.45 at 1 s) | 0.55 | 0.30 to 0.90 | creeps (KE 9e-3 of peak at 2 s) | 17.4 deg |
| a = 0.5, E = 1 MPa | 0.35 (0.20, 0.30) | 0.40 | 0.30 to 0.90 | creeps (2e-3) | 22.4 deg |
| a = 0.5, E = 10 MPa | 0.20 (0.10, 0.15) | 0.25 | 0.30 to 0.90 | creeps (3e-3) | 24.1 deg |
| a = 2, E = 1 MPa | 1.40 (1.00, 1.20) | 1.8 | 1.2 to 3.6 | stops at 1.5 s | 25.5 deg |

Top kept within 0.02 cell in every a = 0.5 run; mass to 1e-14; no stress in empty cells.
Without the cap (a = 0.5, 1 MPa) the run-out is 0.25, as the prototype before decision
0023 (0.25); the cap gives 0.35 (prototype 0.40).

Two processors (a = 0.5, 1 MPa, split across the column at y = 0.06 m) match serial to
round-off at t = 0.1 s (stress 8e-13 of rho g H, velocity 5e-14 m/s; 71,551 face-steps
capped by then). The runs part later (1e-10 Pa at 0.1 s, 1e-7 Pa at
0.2 s, order one at 0.3 s with both rules on; each rule alone stays near 1e-7 Pa at
0.3 s): the order of face sums differs on a processor boundary, and the cap's bisection
and rule 9's gradient comparison are switches that round-off flips; the plastic flow
then amplifies the difference. The test compares at 0.1 s.

## Consequences

- The C++ collapse reproduces the prototype within 0.05 of run-out at a = 0.5 for all
  three E and keeps the E dependence of decision 0023 (0.55 / 0.35 / 0.20): the 10 MPa
  run-out is below the 0019 band, as the prototype's (0.25). a = 2 runs 1.40 against
  1.8, inside the band.
- The a = 0.5 columns still creep at 2 s (front 1 cm per 0.5 s at every E), as the
  prototype at 10 cells per height (decision 0025): faces held at yield slip with no
  rate dependence; a rate-dependent friction (mu(I)) is the physical cure.
- WP 3.3 (outlet, rule 10b) and WP 3.4 (discharge at 0.1, 1, 10 MPa) run on this code
  with both switches on. If the discharge rate depends on E, a cap of the normal jump
  (rule 12) is tried there with this test as its check (decision 0029).
- Serial and parallel runs of a flowing collapse agree only to round-off growth; a
  parallel check of a plastic flow should compare early, or use cases without
  switches.
