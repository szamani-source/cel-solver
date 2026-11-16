# 0036 mu(I) friction plus a plastic-wave impedance, instead of a soft E (WP 3.6 plan)

Date: 2026-11-16. Status: accepted (Sia chose rate-dependent friction over a soft
numerical E on 2026-11-16: "we cannot use parameter that we cannot measure
experimentally").

## Context

WP 3.4b (decision 0035): with rule 12 the C++ discharge rate still falls 16 % from
E = 0.1 to 1 MPa at 5 mm. The fallback of running T0 at E = 0.1 MPa was rejected
because that E is not the material's. Survey: `docs/rate-dependent-friction-survey.md`.

## Decision

1. The friction law becomes mu(I) (Jop, Forterre and Pouliquen 2006) with the
   Barker-Gray (2017) low-I regularisation, inside the Mohr-Coulomb and Drucker-Prager
   returns (rate-dependent return, solved per point). Parameters: d, rho_s, mu_s, mu_2,
   I0, all from standard lab tests. With mu_2 = mu_s it reduces to the present models.
2. The Riemann solver gets an impedance from the elastoplastic tangent at faces where
   the material flows, replacing the after-the-fact caps (rules 8 and 12), which stay
   as switches until the new impedance matches or beats them.
3. Nonlocal fluidity is deferred.
4. A rough estimate (survey) says mu(I) adds about 3 Pa s near the outlet against a
   numerical rho c h of about 200 Pa s at 1 MPa and 5 mm, so step 1 alone is not
   expected to remove the E dependence. WP 3.6b measures this before more is built.

## Consequences

- New packages 3.6a to 3.6e in `docs/work-packages.md`; WP 3.5 (cohesion) waits.
- T0 needs grain data (d, rho_s) and mu(I) constants for Sia's material. Until Sia has
  measurements, glass-bead or sand values from the literature are placeholders.
