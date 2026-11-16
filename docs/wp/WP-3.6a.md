# WP 3.6a: mu(I) friction in the constitutive library

State: done (2026-11-16)
Branch: `claude/project-thread-uprzcj`

## Goal

Drucker-Prager and Mohr-Coulomb take an optional rate-dependent friction mu(I) (Jop law,
optional Barker-Gray regularisation), Python reference first, then C++ and celFoam.

## Read first

- Decisions 0036 and 0037, `docs/rate-dependent-friction-survey.md`
- `prototype/constitutive/mu_i.py`, `src/celConstitutive/kernel/celMuI.H`

## Done when

`tests/test_mu_i.sh` passes: steady simple shear at fixed pressure gives mu(I) to 1e-6
over I = 1e-4 to 1; mu2 = mu_s reproduces the present paths bit for bit (Python and C++);
C++ matches Python.

## Out of scope

- The inclined layer in celFoam (WP 3.6b), the plastic impedance (3.6c), the silo (3.6d).
- The mu_inf high-I term of Barker and Gray (eq. 6.3); compressible mu(I)/phi(I).

## Log

- 2026-11-16: One session. Re-read Barker and Gray (2017) and Jop et al. (2005).
  Finding: for mu_s = tan 35 deg the Jop law has no well-posed range of I, so the
  Barker-Gray branch cannot be built for T0; it is optional (`muIAlpha`) and refused
  where impossible (decision 0037). `mu_i.py` and `check_mu_i.py` (30 s): friction read
  back from the stress equals mu(I) to 2e-14 for I = 1e-4 to 1 (sand Jop: DP psi 0 and
  10, MC psi 0; glass Barker-Gray: DP and MC); mu2 = tan(phi) bit for bit on 25 paths.
  C++ `celMuI.{H,C}`, `Model::updateRate`, DP and MC refactored to a friction-dependent
  return (`returnMap`, `returnAt`); matches Python on 30 mu(I) paths to 1e-12 (an apex
  flag may flip where the cone return lands at the apex to round-off); mu2 = tan(phi)
  equals the plain C++ paths bit for bit (34 paths). celFoam passes dt to the model and
  prints the law. `tests/test_mu_i.sh` about 55 s. Outlet, void and constitutive tests
  unchanged (outlet rate 74.1 kg/s/m as before). Next: WP 3.6b, the inclined layer.
