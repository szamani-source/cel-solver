# 0095: Sand footing: surface cells now dilate; the N_gamma excess is a rate effect

Date: 2026-11-19 (accuracy package A.4, thread "Fix sand footing bearing capacity")
Status: accepted (switch `voidPlasticFlow`, default off; finding recorded).

## Context

Decision 0057 found the smooth footing on phi 30 sand reading N_gamma about 11 at 10 mm
and 9.7 at 5 mm (psi 0) against Martin's exact 7.653 (psi = phi), and associated flow
blowing up. A.3b (decision 0094) removed the hovering-body cause; the associated push
still blew up at about 0.6 s. Sia (2026-11-19): find where the error comes from; cells
must resolve the shear band; mesh dependence is secondary.

## Findings (10 mm cells unless named; implied N_gamma over s = 0.01 to 0.05 B)

**1. Domain size is not the cause.** A layer 1.2 m wide (12 B) reads as the 0.6 m one
(10.6 to 10.9 against 10.6 to 11.2).

**2. The associated blow-up.** The cell that blows up is a partly filled surface cell
in the heave beside the footing (x = 0.145, alpha 0.8, 22 kPa where the geostatic stress
is about 100 Pa, 113 kPa 0.02 s later). In a partly filled cell the void closure
(decision 0021 rule 5, `lagrangeStep` step 5) gives the material only the elastic
unloading of an expanding axis and no strain at all while void closes; the plastic
dilation of the return mapping then has to be taken up by elastic compression, so the
pressure, and with it the strength, grows exponentially with shear. `voidPlasticFlow
yes` lets the positive plastic normal strain of each axis go into the cell's void (a
few fixed-point passes of the model call; the axes share the void, so the material
volume stays within the Lagrangian cell volume: without that share the surface cells
packed to 1.35 rho0 on slow pushes). The associated push now runs to 0.06 B (and 0.1 B),
density 1.001 to 1.002 rho0, mass to 1e-14. psi 0 is unchanged (9.2 to 9.5 at half
speed with and without the switch).

**3. The reading depends on the push speed times the wave impedance.**

| Push (Z v: impedance times speed) | psi 0 | psi 30 |
| --- | --- | --- |
| 0.01 m/s, E 10 MPa (Z v) | 10.6 to 11.2 | 20.1 to 20.4 |
| 0.005 m/s (Z v / 2) | 9.2 to 9.5 | 16.4 to 17.1 |
| 0.0025 m/s (Z v / 4) | (7.3 to 7.5, packed run) | 13.4 to 14.7 |
| 0.00125 m/s (Z v / 8) | (5.4 to 5.8, packed run) | 11.2 to 12.6 |
| 0.005 m/s with E 2.5 MPa (Z v / 4) | (7.2 to 7.3, packed run) | |
| T0 rule set (flowImpedance, caps), 0.01 m/s | 8.6 to 9.3 | 14.2 to 14.8 |
| T0 rule set, 0.005 m/s | 8.3 to 8.7 | 11.9 to 13.0 |

(Runs marked packed used the first version of the switch, which packed surface cells;
they show the trend only.) Halving the speed or quartering E (halving rho c) moves the
reading alike, so it depends on rho c v, not on the step count. Material velocities
scale with the push (no noise). This is the numerical viscosity of the face Riemann
solver and the Hancock predictor that decision 0040 found in flowing mu(I) layers
(about 0.06 rho c_s h per unit strain rate): at 0.01 m/s, rho c v = 1.6 kPa against
0.5 gamma B N_gamma = 6.8 kPa. A frictional material amplifies it (the extra normal
stress raises the strength) and associated flow more (dilation needs normal velocity
jumps, priced at the P impedance), which is why the clay footing (decision 0055) is
right and the sand one is not. With psi 30 the push also stiffens past 0.06 B (q / q_ref
2.2 at 0.08 B, 3.9 at 0.1 B at every speed): dilatancy locking under the footing.

The readings were still falling at the slowest affordable push (25 min at Z v / 8), so
the rate-free value of this mesh is not known.

## Decision

1. Add `voidPlasticFlow` (materialProperties, default off, earlier results unchanged bit
   for bit); the sand footing check uses it. A.4's done-when is
   `tests/test_sand_dilation.sh` (CI): the associated push runs through with a bounded,
   steady load, mass kept, density at most 1.01 rho0.
2. The N_gamma check stays a failing manual test (`tests/manual/test_sand_footing.sh`,
   new options `--speed`, `--E`, `--smax`, `--width`, `--depth`, `--novoidflow`).
3. The excess needs a rate-free treatment of slowly yielding Mohr-Coulomb cells (the
   flow impedance of decision 0040 for a rate-independent law), a package of its own
   (A.5), with the psi 30 locking past 0.06 B.
4. Until then, quasi-static loads on frictional material from celFoam are upper
   estimates; keep rho c v well below the load scale (here gamma B), by a slower push or
   a softer E where the elastic response does not matter.

## Consequences

- Sand bearing results: about +40 % (psi 0) and +160 % (psi 30) at 0.01 m/s, falling
  with speed. Silo discharge is a fast flow and is not affected in the same way.
- `voidPlasticFlow` may help other free-surface cases with dilation (the psi 0.5 to
  1 deg silos Sia asked for); turning it on by default needs the full suite.
