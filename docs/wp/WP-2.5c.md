# WP 2.5c: Mohr-Coulomb tension cut-off (Python reference, then C++)

State: not started
Branch: `wp-2.5c-tension-cut-off`

## Goal

Mohr-Coulomb takes an optional tensile strength `sigma_t` (Rankine cut-off): no principal
stress may exceed it. Without it a cohesive material carries tension up to
`c cot(phi)` at the apex, which is unrealistic for soils and bulk solids. Python first
(the reference), then the C++ kernel matching it.

## Read first

- `docs/decisions/0004`, `0006`, `0011`
- `prototype/constitutive/models.py` (`mohr_coulomb`, `_mc_principal_planes`),
  `check_models.py`
- `src/celConstitutive/kernel/celMohrCoulomb.{H,C}`, `tests/gen_constitutive_reference.py`,
  `tests/test_constitutive_cpp.sh`

## Build

- Python: parameter `sigma_t` (Pa, default absent = no cut-off, results bit-identical to
  now). `0 <= sigma_t <= c cot(phi)`; larger values are rejected with a clear error.
  Add the three Rankine planes `sigma_i - sigma_t <= 0` (associated flow) to the
  principal-space return as further active planes, combined with the Mohr-Coulomb planes
  by the same active-set linear solve (tension plane, tension edge, tension apex, and the
  mixed Mohr-Coulomb/tension corners). History gains nothing new; plastic strain from
  tension goes into `eps_p` and `eps_p_eq` like the rest.
- Checks in `check_models.py`: uniaxial tension strength equals
  `min(sigma_t, 2 c cos(phi)/(1 + sin(phi)))`; hydrostatic tension is capped at `sigma_t`;
  with `sigma_t = 0` no principal stress is ever positive; 10,000 random increments end
  with every yield function (Mohr-Coulomb and Rankine) at or below tolerance; unloading
  from the cut-off is elastic; compression paths are unchanged by the cut-off.
- C++: the same in `celMohrCoulomb`, wrapper keyword `sigmaT` (optional), new reference
  paths (uniaxial tension, hydrostatic tension, tension then shear, random increments
  biased to tension, for `sigma_t = 0` and for `0.5 c cot(phi)`).
- Update `docs/constitutive-models.md` and the known issue in `docs/STATUS.md`.

## Done when

- `tests/test_proto_constitutive.sh` passes with the new checks.
- `tests/test_constitutive_cpp.sh` passes: C++ matches Python to a relative 1e-10 on the
  new paths; all existing paths unchanged.

## Out of scope

- Tension cut-off for Drucker-Prager (add a row to `docs/work-packages.md` if wanted),
  softening after tensile failure, void opening (WP 3.2).

## Log

(none yet)
