# 0028 Normal jump cap (rule 12) and dispersed interior void (rule 13) for the discharge

Date: 2026-11-15. Status: accepted for the prototype silo discharge (WP 3.3r). Not yet
for the column collapse (see Consequences).

## Context

After the outlet fix (decision 0027) the T0 discharge rate still fell with the
stiffness: 81 / 73 / 54 kg/s/m at E = 0.1 / 1 / 10 MPa (10 mm cells, D = 0.06 m). At
1 MPa it also fell with refinement (73 at 10 mm, 61 at 5 mm), while at 0.1 MPa the
2.5 mm mesh drained 14 percent faster than 5 mm.

## What was measured (10 mm, D = 0.06 m, first 0.4 to 0.6 s after opening)

- Not the time step: at 1 MPa, CFL 0.1 instead of 0.5 gave the same rate (66 against
  70 kg/s/m at 0.25 s).
- Not slower flow alone: at 1 MPa the outlet speed was higher than at 0.1 MPa (0.92
  against 0.85 m/s), but the outlet zone (|x - centre| < D, y < 6 cm) held 4 percent
  void and the bottom outlet rows 10 to 20 percent (0.1 MPa: 0.5 and 1 to 4 percent);
  at 10 MPa 7 to 11 percent, with half the cells there partly filled.
- Separation tested on the Riemann stress without its jump term: less void, but the
  rate did not rise (67 kg/s/m at 1 MPa): the void was also a release for a resistance.
- The resistance is the normal jump term z (vR - vL) / 2 of the acoustic Riemann
  solver, which rule 8 (decision 0023) left uncapped. In plastic flow at psi = 0 the
  face normal velocity jumps are mostly deviatoric extension (dvx/dx = -dvy/dy), so
  the term is a numerical viscosity on plastic flow of order rho c h, growing with
  sqrt(E). At 10 MPa, z dv for dv of 0.05 to 0.2 m/s between cells is 8 to 30 kPa,
  ten times the stress near the outlet (about 1 kPa).
- The void comes from rule 4 (decision 0021): a sweep stretches a cell by dt dv/dy,
  which at 10 MPa is about 8 times the elastic strain the material has near the outlet
  (s / M), so void opens in every sweep. A cell with even 1e-6 of void then leans its
  slab towards the fuller neighbour (rule 2) and the face on the other side becomes a
  free surface with zero traction, so the flowing bulk is cut into blocks.

## Decision 1: rule 12, the normal jump term is capped at yield

As rule 8 for the shear jump: at each interior face the normal jump term is scaled
back (bisection, 20 steps) so that the face stress (normal and shear tractions, mean of
the two cells' tt and zz stresses) is no further outside Mohr-Coulomb than without the
jump. The capped traction is used for the momentum of both cells (one traction per
face: momentum exact) and the face keeps one normal velocity for the geometry and the
remap (mass exact). For the strain, each full cell takes the velocity of its own
characteristic at the capped traction, as a capped shear face slips; the dropped
length strain dstr is split: its deviatoric part is dropped (plastic flow at psi = 0
keeps volume), its volume is kept as an in-plane isotropic strain (dstr / 2 in each
in-plane direction), so the pressure stays tied to the density.

Tried first: dropping the whole dstr (uniaxial). The rate rose but the material near the
outlet lost 11 percent of its density in 0.4 s (dilation without stress). Keeping the
volume removed the drift.

## Decision 2: rule 13, void inside the bulk is dispersed

In a sweep, a partly filled cell whose two neighbours both hold at least 0.5 of
material is spread (lean 0, rule 2 otherwise unchanged): it touches both faces, so its
faces are contact faces. Small voids opened by one sweep inside the flowing bulk are
porosity, closed first by the next compression (rule 4), instead of internal free
surfaces. 0.9 instead of 0.5 gave the same rate (86 against 86 kg/s/m at 1 MPa).

Each rule alone is not enough: rule 13 alone at 1 MPa let the void grow (to 15 percent
in the outlet zone) and the rate fell to 50 to 75 kg/s/m; rule 12 alone gave 85 at
1 MPa and 71 at 10 MPa with void still growing. Together the outlet zone holds 0.1 to
0.5 percent void at every E and the rate is the same at CFL 0.1 and 0.5.

## Result (`tests/test_proto_silo_stiffness.sh`, 10 mm, D = 0.06 m, rate 97 to 88 percent)

| E (MPa) | Before (0027) | Rate (kg/s/m) | C | 0.1 s windows | Corner / outlet speed (mean) |
| --- | --- | --- | --- | --- | --- |
| 0.1 | 81 | 91.2 | 1.10 | 5.1 % | 0.63 % |
| 1 | 73 | 87.1 | 1.05 | 4.2 % | 0.10 % |
| 10 (T0) | 54 | 84.0 | 1.01 | 4.0 % | 0.09 % |

Spread 8.7 percent (band 10 percent: thin margin). The stiff runs are also quieter
(instantaneous corner speed 0.25 percent at 10 MPa against 2.3 at 0.1 MPa).

With the rules on, the other discharge tests still pass (bands unchanged):

| Check | WP 3.3q | Now |
| --- | --- | --- |
| Rate law, 5 mm, E = 0.1 MPa: C for D = 0.04 / 0.06 / 0.08 m | 0.99 / 0.99 / 1.06 | 1.14 / 1.16 / 1.12 |
| Exponent (1.3 to 1.7), fitted C (1.0 to 1.9) | 1.59, 1.035 | 1.48, 1.133 |
| Head independence, 5 mm | 0.992 | 0.992 |
| Corner / outlet speed, 5 mm, D = 0.06 m (mean) | 0.83 % | 0.67 % |
| Outlet test, 10 mm: C at D = 0.06 m | 0.98 | 1.10 (head 0.945) |
| Mesh, D = 0.06 m, E = 0.1 MPa, rate 0.2 to 0.5 s: 5 against 2.5 mm | 81.5 against 93.3 (14 %) | 96.2 against 93.4 (3 %) |

The exponent is now closer to 1.5 and the fine-mesh rise is gone: the 3.3q mesh
dependence at 0.1 MPa was the same normal-jump viscosity (of order rho c h, smaller on
the finer mesh).

## Switches

`VoidBox.normal_cap` (rule 12) and `VoidBox.disperse` (rule 13) default to off, so the
void, collapse and settling tests are unchanged. `discharge2d.make_tank` switches both
on (`NORMAL_CAP`, `DISPERSE` = 0.5); `--old-rules` reproduces WP 3.3q.

## Consequences

- The discharge can now run at T0's own E = 10 MPa (10 times more steps than 0.1 MPa).
- Column collapse with both rules (a = 0.5, 10 cells, E = 1 MPa): run-out 1.05 against
  0.40 without (band 0.30 to 0.90), fitted slope 9.5 deg against 22.
  At 0.1 and 10 MPa: 0.8 and 1.2 (without: 0.45 and 0.25). The rules make the collapse
  flow too far and it still depends on E (other way round).
  So the rules are not yet general: the jump terms were also what stopped the collapse.
  A follow-up package must reconcile them before the C++ Lagrange step copies either.
- The C++ Lagrange step is unsplit (decision 0013): rule 13's trigger (void opened by one
  sweep's stretch) may not arise there, but rule 12's viscosity argument applies as is.
