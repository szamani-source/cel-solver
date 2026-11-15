# 0029 No single rule set for collapse and discharge in the split prototype (WP 3.2t)

Date: 2026-11-15. Status: accepted. Rules 12 and 13 stay as decision 0028 left them: on
in the discharge, off in the column collapse. The question moves to the C++ Lagrange
step (unsplit, decision 0013).

## Context

Decision 0028 made the T0 discharge rate independent of E with rule 12 (normal Riemann
jump capped at yield) and rule 13 (void inside the bulk spread as porosity). With both
on, the a = 0.5 column collapse (10 cells) runs out 0.8 / 1.05 / 1.2 at E = 0.1 / 1 /
10 MPa against a band of 0.30 to 0.90. WP 3.2t looked for one rule set that passes
both `tests/test_proto_collapse.sh` (with run-outs at 1 and 10 MPa within 0.1) and
`tests/test_proto_silo_stiffness.sh`.

## What was measured

Collapse: a = 0.5, 10 cells, `void2d.column_collapse(0.5, cells_short=10, E=...,
normal_cap=..., disperse=...)` with the study knobs below (module constants in
`prototype/void2d.py`). "Runs away" = still at peak kinetic energy at t = 2 s and stopped
by the far wall (run-out 1.2). Discharge: 10 mm, D = 0.06 m, rate from 97 to 88 percent
of the mass, `discharge2d.run_case(cell=0.01, E=..., window=CHECK_WINDOW)`. With
decision 0028's rules the discharge gives 91.2 / 87.1 / 84.0 kg/s/m at 0.1 / 1 / 10 MPa.

### Each rule alone (collapse run-out at 0.1 / 1 / 10 MPa)

| Rules | 0.1 MPa | 1 MPa | 10 MPa |
| --- | --- | --- | --- |
| Neither (CI setting) | 0.55 | 0.40 | 0.25 (STATUS, measured before WP 3.2r) |
| Rule 13 only | 0.80 | 0.85 (top lost) | runs away |
| Rule 12 only | 0.95 | runs away | (not run; 1 MPa already fails) |

Both rules break the collapse on their own, for different reasons.

### Rule 13: porosity diffuses through the bulk

Counted partly filled cells with all four neighbours at least half full (1 MPa): without
rule 13 they stay below 7 and hold no void after 0.3 s; with it they grow to 90 cells
holding half a cell of void, and the kinetic energy keeps rising. A spread cell gives
material through both faces in the remap, so its void spreads into its full neighbours,
and every partly filled cell has had its normal stress unloaded by rule 4. The bulk
becomes a frictionless porous layer. The void itself comes from the split sweeps: in
isochoric plastic flow one sweep stretches a cell (rule 4 unloads it and opens void)
and the other compresses it back.

Knob `STRETCH_ISO`: in a full cell the part of a stretch matched by the transverse
contraction rate (from the transverse velocity gradient) goes to the material; only the
net in-plane expansion opens void. Rule 13 plus `STRETCH_ISO` (rule 12 off) gives a
normal collapse: 0.40, top kept, fitted slope 26.3 deg at 1 MPa. But it slows the
discharge: rules 12 and 13 plus `STRETCH_ISO` give 76.1 kg/s/m at 1 MPa (87.1 without),
and without rule 12 the 10 MPa discharge is 67.6 kg/s/m with 0.1 s windows 28 percent
apart.

`DISPERSE_4` (rule 13 also needs the cross-sweep neighbours) and `DISPERSE_CONTACT`
(rule 13 changes the contact faces only, not the remap slab) did not stop it (0.65 and
1.0 at 1 MPa); `DISPERSE_SELF` = 0.95 (only nearly full cells) made it worse (1.15).

### Rule 12: no capped variant stops the collapse at 10 MPa

| Variant (knob) | Collapse 1 MPa | Collapse 10 MPa | Discharge 1 MPa |
| --- | --- | --- | --- |
| Rule 12 as in 0028 | runs away | runs away | 87.1 |
| Viscous floor eta = 200 Pa s (`NCAP_ETA`) | 0.40 | runs away | 73.8 (no rule 13) |
| Volumetric part of the jump kept (`NCAP_VOL`) | runs away | | 57.4 at 10 MPa |
| Only above 2 or 5 1/s (`NCAP_GAMMA`), with 13 and ISO | | runs away | |
| Only between full cells (`NCAP_FULL`), with 13 and ISO | 0.85 | | |
| Only where both cells are at yield (`NCAP_YIELD` 0.05), with 13 and ISO | 0.35 | runs away | 48.3 |
| Capped cells strain with the geometry (`NCAP_STRAIN` common) | runs away | | |
| Capped faces keep an odd-even jump (`NCAP_HYPER`) | runs away | runs away | 90.3 |
| Rule 12 as in 0028, CFL 0.15 instead of 0.5 | 0.90, stops | | |

(eta = 200 Pa s is about 80 percent of the jump viscosity z h / 2 at 1 MPa, 26 percent
at 10 MPa.) `STRETCH_CONFINED` (any stretch of a full cell pressed across goes to the
material) stopped the collapse with `NCAP_VOL` (0.50 / 0.35 at 1 / 10 MPa) but the
discharge fell to 31 kg/s/m: material at the outlet must open void to fall.

Findings:

1. With the plain cap, the collapsing column turns hydrostatic (1 MPa, median over full
   cells under more than 200 Pa: horizontal / vertical stress 0.8 to 1.1 against 0.30
   without the cap, q / p 0.37 against 0.64) and the kinetic energy reaches 20 times the
   peak without the cap. Capped compressive faces drop the deviatoric part of the strain
   while the cells are still elastic, so the deviator ratchets down.
2. Restricting the cap to cells at yield stops that (q / p stays 0.62 at 10 MPa) but the
   10 MPa column still never stops: the kinetic energy stays at 0.05 to 0.1 J instead
   of decaying. The result depends on the time step (1 MPa: runs away at CFL 0.5,
   stops at 0.90 at CFL 0.15), so the normal jump is also what keeps the split
   MUSCL-Hancock step stable at faces at yield; the brief's "the jump terms were the
   only damping of the slow front" is right.
3. The discharge needs the opposite: every variant that keeps the collapse frictional
   (cap only at yield, isochoric stretch kept, a viscous floor) brings the rate back
   down to 48 to 76 kg/s/m at 1 MPa. The E independence of decision 0028 therefore
   comes partly from capping faces whose cells are still elastic, which softens the
   material near the outlet. The Beverloo constant and exponent there are still within
   their bands, but that rate should not be read as a pure removal of numerical
   viscosity.
4. A fixed viscosity (eta) or a rate threshold cannot serve both: the collapse at
   10 MPa needs at least about a third of the full jump (an E-proportional damping),
   and the discharge needs much less.

## Decision

- No change to the defaults: rules 12 and 13 stay on in the discharge
  (`discharge2d.make_tank`) and off elsewhere, as decision 0028 left them. WP 3.2t is
  closed without meeting its done-when; its test was not added.
- The study knobs stay in `prototype/void2d.py`, all off, so the tables above can be
  reproduced: `NCAP_ETA`, `NCAP_GAMMA`, `NCAP_FULL`, `NCAP_YIELD`, `NCAP_VOL`,
  `NCAP_STRAIN`, `NCAP_HYPER`, `STRETCH_ISO`, `STRETCH_CONFINED`, `DISPERSE_4`,
  `DISPERSE_SELF`, `DISPERSE_CONTACT`. `column_collapse` takes `normal_cap` and
  `disperse`.
- The C++ Lagrange step is unsplit (decision 0013), so neither the one-sweep void
  (rule 13's trigger) nor the split-step instability of finding 2 need arise there. The
  C++ silo packages copy the jump terms in full (no rule 12 or 13) first, and measure
  the discharge rate against E and the collapse with the same code. Only if the C++
  discharge rate then depends on E is a cap added there, with the collapse as its
  check from the start.

## Consequences

- The prototype keeps two rule sets. Its discharge numbers (decisions 0027, 0028) and
  collapse numbers (0023 to 0025) stay valid for their own rule sets; finding 3 is
  recorded against the 0028 rates.
- The stiffness margin of `tests/test_proto_silo_stiffness.sh` (8.7 percent against
  10) was not improved. `NCAP_HYPER` gave 90.3 kg/s/m at 1 MPa (87.1 with 0028), not
  run at 0.1 or 10 MPa.
- mu(I) friction was not built: with T0's phi = 35 deg and a grain of about 1 mm the
  inertial number in the collapse is about 0.02, which adds about 0.02 to the friction
  coefficient, far too little to stop a run-away of this kind.
