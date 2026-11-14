# 0009 Drucker-Prager in C++: word options beside numeric parameters (WP 2.5a)

Date: 2026-11-14. Status: proposed (WP 2.5a), confirm by merging the PR.

## Context

Drucker-Prager needs a non-numeric choice, the match to Mohr-Coulomb (`fit`), but the
kernel's `Parameters` map (decision 0006) holds numbers only. The model also needs
reference paths that actually reach the cone and the apex.

## Decision

- The kernel gets a second map, `Options` (`std::map<std::string, std::string>`), and
  `makeModel(name, params, options = {})`. Models that need no options ignore it, so
  existing calls are unchanged. Unknown option values throw `std::invalid_argument`.
- The OpenFOAM wrapper constructor takes a third list of optional word keys; the
  `druckerPrager` wrapper lists `fit`, so dictionaries say `fit compression;`
  (default `compression`, also `extension`, `plane_strain`).
- The reference file gets `option <key> <word>` lines after `case`; the generator
  writes any string-valued Python parameter that way. The generator now names the C++
  model on each case and maps it to the Python function (`PY_MODELS`).
- `plasticCorrector` (the C++ `_finish`) is a shared kernel function, so Mohr-Coulomb
  (WP 2.5b) can reuse it.
- Constructor checks: `0 <= psi <= phi < 90`, `c >= 0`, and `c > 0` when `phi = 0`.
- Drucker-Prager uniaxial-strain path starts from the isotropic 100 kPa state and goes
  from 1 % compression on into 1 % extension, then reloads: uniaxial compression alone
  stays inside the cone for these friction angles (q/p ratio too low), so the J2 path
  would never yield. Random Drucker-Prager paths must contain apex returns.

## Consequences

- One more argument to `makeModel`; other non-numeric choices (e.g. a future hardening
  law name) use the same mechanism.
- `compareReference` now fails a Drucker-Prager random path with no apex return, so a
  change to the random sequence that misses the apex is caught.
