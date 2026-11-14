# WP 2.5c: Mohr-Coulomb tension cut-off (Python reference, then C++)

State: done
Branch: `claude/project-thread-m6f254` (thread sessions push only to their assigned branch)

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

### 2026-11-14

- Python: `mohr_coulomb` takes optional `sigma_t`; `mc_tension_cutoff` validates it,
  `_mc_tension_return` does the active-set return over {main, two edges, T1, T2, T3},
  `rankine_yield` is the cut-off yield function. The old Mohr-Coulomb return runs first
  and stands if it satisfies the cut-off, so models without `sigma_t` are bit-identical
  (checked on 20,000 random increments against the previous code). Decision 0015.
- `check_models.py` section 7: uniaxial tensile strength min(sigma_t, 2c cos/(1+sin)) for
  sigma_t = 0, 0.5 and 0.9 c cot(phi); hydrostatic tension capped at sigma_t (apex flag);
  unloading elastic; compression paths unchanged; out-of-range sigma_t rejected;
  4 x 10,000 tension-biased random increments inside both surfaces (max s1 with
  sigma_t = 0 is 2e-10 Pa, round-off); associated return equals an independent
  brute-force closest point to 1e-15. Whole script takes about 45 s.
- C++: `MohrCoulomb` reads optional `sigmaT`, `tensionReturn` copies the Python order;
  new `Return` values tensionPlane, tensionEdge, tensionApex, tensionCorner. Reference
  generator writes 10 new paths (sigma_t = 0 and 0.5 c cot(phi): triaxial tension at
  2 kPa confinement, hydrostatic tension, tension then shear, uniaxial strain,
  1000 tension-biased random increments). All match to 2e-14 and each requires the
  tension return it is there for; the random paths hit all four. Wrapper reads
  `sigmaT`; wrapperCheck and celFoam (smokeBox) check it, and that sigmaT > c cot(phi)
  stops celFoam.
- Learned: the sigma_t = 0 uniaxial tension path from zero confinement is all round-off
  (stress ~1e-12 Pa), so a relative comparison is meaningless; the reference path uses
  2 kPa confinement instead.
- Not done (out of scope): cut-off in the vectorised `batch.py`, Drucker-Prager cut-off,
  softening after tensile failure.
