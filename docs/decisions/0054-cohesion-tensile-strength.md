# 0054 Cohesion in the silo: tensile strength at faces and in void opening (WP 3.5)

Date: 2026-11-18. Status: accepted (WP 3.5 done at 10 mm).

## Context

T0 started cohesionless (c = 0, decision 0016); Sia asked for cohesion to be added
later. The constitutive side already had it: Mohr-Coulomb and Drucker-Prager take `c`,
Mohr-Coulomb an optional Rankine cut-off `sigmaT` (decisions 0009, 0011, 0015), the
shear and normal caps (rules 8 and 12) use c cos(phi), and the mu(I) return keeps c
with phi_m. The free-surface rules did not: a material counted either as "no tension"
(c = 0 or sigmaT = 0: contact faces separate as soon as they pull, cells open void as
soon as their normal stress reaches zero) or as having unlimited tensile strength
(faces never separate, a full cell never opens void). A cohesive soil was in the
second group, so it could never tear: the material under an arch stretched instead of
dropping off.

Done-when (`docs/work-packages.md`): flow stops below a critical outlet width.

## Decision

1. **One tensile strength sT** per material (`eulerianSolid::tensileStrength_`):
   0 without tensile strength (as before), `sigmaT` if given, else c cot(phi) (the
   Mohr-Coulomb apex, the largest tension the Rankine cut-off allows) for a cohesive
   Mohr-Coulomb or Drucker-Prager material, unlimited for elastic and J2. Printed in the
   log ("tensile strength 714.074 Pa" for c = 500 Pa, phi 35).
2. **Contact faces tear at sT:** a face where both sides touch keeps the Riemann
   solution while its normal traction is at most sT; beyond it the face is traction-free
   (rule 3 of decision 0021 with 0 replaced by sT).
3. **Void opening at sT:** expansion of a cell first unloads its normal stress to sT,
   the rest opens void (rule 5 with 0 replaced by sT). For a full cell this applies only
   to materials with a finite sT; elastic and J2 full cells still never open void, and
   partly filled cells of those materials still use 0, as before.
4. Walls and outlets keep no adhesion (material separates from a wall at zero tension).
   Wall adhesion is not in scope.

Cohesionless and elastic runs are unchanged bit for bit (sT = 0 and unlimited give the
old expressions); the whole CI suite passes unchanged.

## Findings

T0 tank, W 0.4 m, fill 0.3 m (the arch does not depend on the head), 10 mm cells, T0
rule set (decision 0050), walls 20 deg, outlet opens at 0.3 s. Reference: Jenike's
critical width of a slot, B_c = sigma_c / (rho g) with the Mohr-Coulomb unconfined
strength sigma_c = 2 c cos(phi) / (1 - sin(phi)): 109 mm at c = 500 Pa, 55 mm at 250 Pa.
"Stops": rate over the last 0.15 s below 5 % of the Beverloo rate with C = 1.

| c (Pa) | D (mm) | E (MPa) | rate / Beverloo C=1 | discharged in 0.5 s | |
| --- | --- | --- | --- | --- | --- |
| 500 | 80 | 1 | 0.007 | 0.004 | stops |
| 500 | 100 | 1 | 0.008 | 0.039 | stops (after a short drop) |
| 500 | 120 | 1 | 0.44 | 0.16 | flows |
| 500 | 140 | 1 | 0.78 | 0.34 | flows |
| 0 | 80 | 1 | 0.85 | 0.22 | flows (control) |
| 500 | 100 | 10 | about 0.01 | 0.005 | stops |
| 500 | 120 | 10 | flows | 0.17 | flows |
| 250 | 40 | 1 | 0.045 (creeps) | 0.005 | creeps |
| 250 | 60 / 80 | 1 | 0.27 / flows | 0.04 / 0.15 | flows |

The critical width at c = 500 Pa is between 100 and 120 mm (Jenike 109 mm), the same
at E = 1 and 10 MPa: the arch is set by strength, not stiffness. At c = 250 Pa it is
between 40 and 60 mm (Jenike 55 mm), so it scales with c.

Without items 2 and 3 (old code, same c = 500 Pa): the 100 mm tank creeps at 0.041 of
Beverloo (borderline, not stopped) and the 120 mm tank flows at 0.30; with them 0.008 and
0.44. Tearing makes the transition sharp: material under the arch drops off instead of
stretching.

## Consequences

- `tests/test_cohesion.sh` (CI, about 80 s on 4 cores) is the done-when: c = 500 Pa stops
  at 80 and 100 mm, flows at 120 and 140 mm; c = 0 flows at 80 mm; mass, empty cells.
  `tests/manual/test_cohesion.sh` repeats the bracket at E = 10 MPa and at c = 250 Pa.
- Tutorial `tutorials/siloCohesion`: T0 tank (0.6 m fill, E 10 MPa) with c = 500 Pa,
  outlets 100 mm (arches) and 140 mm (flows).
- A stopped tank still creeps at 0.5 to 2 kg/s/m (under 1 % of Beverloo at c = 500 Pa,
  4.5 % at c = 250 Pa and 40 mm, steady over 1 s). Same family as the creep of the
  inclined layer (decision 0044). Recorded as a known issue.
- Only 10 mm checked. The arch is a few cells wide, so the critical width may shift with
  the mesh; a 5 mm bracket is the obvious next check.
- Drucker-Prager uses the same c cot(phi), which is not its own apex for every fit
  (the cone has no Rankine cut-off). A `sigmaT` entry there sets only the face and
  void-opening strength. Only Mohr-Coulomb is checked.
