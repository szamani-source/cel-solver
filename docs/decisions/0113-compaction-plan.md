# 0113 Compaction before the Linux install: a volumetric cap (plan, packages K.1-K.3)

Date: 2026-11-20. Status: accepted (Sia chose compaction over the collapse
toe, 2026-11-21 03:29 UTC); K.1 and K.2 done.

## Context

Sia (2026-11-20 16:04 UTC): compaction (a cap, or density-dependent
behaviour) is the third item before the Linux install. Decision 0108 showed
why: filling the Couto silo from the top doubles the wall pressures because
Mohr-Coulomb and Drucker-Prager have no volumetric yield, so an impact
locks stress in and the bulk never densifies. Sia's rules: parameters must
be lab-measurable (decision 0036) and come from published sources; develop
on simple cases first.

## Proposal

A hardening volumetric cap on both frictional models, with the standard
oedometer parameters (lambda, kappa, initial void ratio e0 or density):

- elastic bulk modulus from the swelling line, K = (1 + e) p / kappa
  (floored by the given E at low p), shear modulus from nu;
- an elliptical cap closing the existing cone (modified Cam clay shape in
  p-q), with its size p_c hardening with plastic volumetric strain:
  d eps_v^p = (lambda - kappa)/(1 + e) d p_c / p_c;
- p_c is a new history variable, remapped mass-weighted like epsPEq (the
  History struct, the Eulerian field list and restart need it);
- the cone, mu(I), the nonlocal law and cohesion stay as they are; without
  `lambda` the models are unchanged bit for bit.

| WP | Title | Done when |
| --- | --- | --- |
| K.1 | Cap in the kernel | Single-point oedometer (1D confined compression, load and unload) gives e = e0 - lambda ln(p/p0) on loading and kappa on unloading, within 1 %; MC and DP; isotropic compression likewise |
| K.2 | Oedometer in the Eulerian solver | A column compressed by a lid (confined, frictionless walls): load-unload curve against K.1 within 2 %; density and mass conserved; restart and 2 processors bit for bit |
| K.3 | Settling column (simple case for filling) | A column filled under gravity compacts: density profile against the 1D closed form within 5 %; at-rest K and wall pressure no longer depend on dropping height (decision 0108 effect gone) |

## Consequences

- One more history variable per cell (remap cost, restart files).
- The silo comes back to filling-from-the-top only after K.3 passes.

## What was built (K.1, K.2, 2026-11-21)

The proposal above changed in two ways while building it:

- **Linear elasticity, not K = (1 + e) p / kappa.** The cap's hardening is
  written so that loading along p = p_c gives the total compression line
  exactly whatever K is: eps_v = lambda* ln(p/p0) with lambda* = lambda/(1 + e0)
  (entry `capLambda`), p_c,new = p_c exp((P' - p_c)/(K lambda*)) in closed
  form. Unloading is elastic with the model's K (E and nu: the oedometer
  unloading modulus), so no kappa entry and the frictional models keep their
  constant elasticity. Needs K lambda* > P (checked: the step stays elastic
  otherwise). Stress scaled towards the apex by p_c,new/P' (stays inside a
  cohesionless cone). `capP0` (default 100 Pa): the smallest p_c, so
  stress-free material is elastic up to P0 and then on the line.
- **K.2 is a settling column, not a lid oedometer.** A lid pushing into a
  fixed Eulerian mesh would feed material through the boundary; a column
  compacting under its own weight into the void above tests the same
  coupling (density, stress and p_c through the remap). It needs
  `fullBlend yes`: with the plain partly filled rule (decision 0024 rule 4)
  the stress-free surface cell jolts the column each time the settling
  surface crosses a cell face, and the cap keeps every overshoot (cells 5 to
  18 % above the line near the top and the base; with fullBlend 0.0 %).

Results:

- `tests/test_cap.sh` (seconds, kernel): isotropic and oedometer loading
  to 100 kPa on the line to round-off, unloading K to 4e-14, MC and DP.
- `tests/test_settle_column.sh` (about 15 s): 47 cells on
  ln(rho/rho0) = P0/K + lambda* ln(P/P0) within 1.4 % (P 0.4 to 10.8 kPa),
  mass 2e-15, at rest; the top settles 36 mm (0.7 mm without the cap).

Fields: `pc` [Pa] per cell (written, read on restart, remapped mass-weighted
as epsPEq; fed material arrives with pc 0). Not yet with immersed bodies
(fatal error). K.3 next: filling from the top (decision 0108).

## K.3 study (2026-11-21): the cap does not remove the top-fill doubling

`tests/compaction/fill_box.py`: the WP 3.7a inlet box, 0.6 m tall, fed at
2 m/s through a 40 mm lid inlet for 1 s, settled with damping, against the
same mass placed at once under a 0.5 s gravity ramp (T0 rules, fullBlend;
cap lambda* = 0.02):

| Start | Cap | K mid-deposit | Wall s_xx below mid-height | Mean density |
| --- | --- | --- | --- | --- |
| gravity on | no | 0.44 | 1.4 to 2.4 kPa | 1800 |
| gravity on | yes | 0.42 | 1.4 to 2.4 kPa | 1896 |
| filled from the top | no | 0.86 | 3.4 to 4.3 kPa | 1859 |
| filled from the top | yes | 0.87 | 3.3 to 3.9 kPa | 2015 |

The cap compacts the landing material (12 % denser than placed, 6 % denser
than the gravity-on cap run) but the locked-in lateral stress stays: K is
still about twice the gravity-on 0.43. The impact builds lateral stress in
the elastic range of each landing cell (the stream lands at about 3 m/s,
rho c v of order 100 kPa), and nothing in MC or the cap relaxes a lateral
stress below the cone afterwards. Decision 0108's doubling is therefore not
a missing-compaction effect alone; it needs either a slower landing (a
realistic filling rate: the Couto silo is filled slowly) or a stress
relaxation of material at rest. K.3's done-when is not met; open.
