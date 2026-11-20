# 0113 Compaction before the Linux install: a volumetric cap (plan, packages K.1-K.3)

Date: 2026-11-20. Status: proposed (starts when Sia picks compaction over
the collapse toe, decision 0111).

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
