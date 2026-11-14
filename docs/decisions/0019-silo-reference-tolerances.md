# 0019 Silo reference values and test tolerances

Date: 2026-11-14. Status: accepted (WP 3.0).

## Context

The silo packages (3.1p, 3.2p, 3.3p, then C++ 3.1 to 3.5) need numbers to pass or fail
against. The literature values (`docs/spec/silo-reference.md`) come from materials and
constitutive laws that are not ours: continuum discharge constants were all computed
with mu(I) at a static friction of 0.32 to 0.38, experiments used glass beads (repose
22 deg), while T0 is rate-independent Mohr-Coulomb with phi = 35 deg. Some prefactors
(Lajeunesse et al. 2005 column collapse, Jaky's formula) could not be read in an opened
source.

## Decision

1. Settled state: the smooth-wall stress ratio to expect is nu/(1 - nu) = 0.4286
   (elastic, no yield), not Jaky; with frictional walls the fitted Janssen K must lie
   between the wall-yield value (0.29) and the elastic value (0.43), tested as 0.27 to
   0.45. Janssen is applied to the 2D bin with wall length per area 2/W.
2. Discharge: the reference constant is C = 1.45 in `Wdot = C rho sqrt(g) D^(3/2)` with
   k = 0 (centre of continuum values 1.40, 1.42, 1.55). The test band is 1.0 to 1.9
   (30 percent), wide because the effect of phi = 35 deg on C is not known. The exponent
   band (1.3 to 1.7) is the real scaling test. Results outside the band are recorded,
   not tuned away.
3. Column collapse: the reference run-out is the quasi-2D experimental fit
   1.2 a (Lube et al., as quoted by Staron and Hinch and by Dunatunga and Kamrin). The
   band is 0.5 to 1.5 times the fit, because continuum mu(I) codes give about 2.2 a,
   DEM 2.5 a, and higher friction shortens the run-out. The sharp checks are geometric:
   H_f = H_i at a = 0.5, and deposit slopes not steeper than phi.
4. Values not found in an opened source stay marked "not verified" in the spec and are
   not used as test limits.

## Consequences

- Wide bands on C and run-out mean those checks catch gross errors (wrong exponent,
  missing friction, mass loss), not 10 percent ones. A tighter validation would need a
  mu(I) model (a later constitutive package) or reference runs with our own material in
  another code.
- 3.1p's tall-fill case (H = 3 W) separates the two K limits by 14 percent in bottom
  stress, more than its 10 percent tolerance, so the fitted K is meaningful.
- The Beverloo steady window is short for D = 0.08 m at H = 0.6 m (0.8 s); 3.3p measures
  the rate between 95 and 70 percent of the initial mass or also uses H = 0.9 m.
