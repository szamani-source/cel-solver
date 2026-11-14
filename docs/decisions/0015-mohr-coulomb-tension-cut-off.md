# 0015 Mohr-Coulomb tension cut-off: Rankine planes, active-set return (WP 2.5c)

Date: 2026-11-14. Status: proposed (WP 2.5c), confirm by merging the PR.

## Context

Without a cut-off, cohesive Mohr-Coulomb carries tension up to the apex c cot(phi)
(17.3 kPa for c = 10 kPa, phi = 30), far more than soils and bulk solids can hold
(decisions 0004, 0011). The brief asks for an optional tensile strength `sigma_t`, a
Python reference first and a C++ kernel that matches it.

## Decision

- **Surface.** Three Rankine planes `s_i - sigma_t <= 0` with associated flow, added to
  the Mohr-Coulomb planes. Python parameter `sigma_t`, C++ kernel parameter and
  dictionary keyword `sigmaT` (Pa). Optional: absent means no cut-off.
  `0 <= sigma_t <= c cot(phi)` is required (any value is allowed for phi = 0); larger
  values would be inactive and are rejected with a clear error, as are negative ones.
- **Return.** The existing Mohr-Coulomb return (main plane, edge, apex) runs first. If
  its result satisfies the cut-off it stands. So a path that never reaches the cut-off
  is unchanged bit for bit, and a model without `sigma_t` is exactly the WP 1.4 / 2.5b
  model (checked: 20,000 random increments identical to the old code; compression
  paths identical with `sigma_t = 0`).
- Otherwise (the trial violates only the cut-off, or the Mohr-Coulomb result does) the
  return tries every active set of one, two, then three planes out of
  {main, s1 = s2 edge, s2 = s3 edge, T1, T2, T3} (41 sets) in a fixed lexicographic
  order and takes the first admissible one: multipliers >= 0, principal order kept,
  Mohr-Coulomb and Rankine yield functions <= 0 (tolerance 1e-10 of the stress scale).
  This is the Kuhn-Tucker system of the multi-surface return. It covers the tension
  plane, tension edge, tension apex and the mixed Mohr-Coulomb/tension corners without
  case-by-case formulas, and costs at most 41 tiny solves, only on tension steps.
  Planes with dependent normals are skipped by `|det A| <= 1e-12 prod |A_ii|`.
- **Check of the algorithm.** For associated flow the result is the closest point in the
  energy norm. `check_models.py` compares it with an independent brute-force projection
  onto every face of the full nine-plane surface in unsorted principal space: they agree
  to 1e-15 over 1000 trial states, including repeated eigenvalues.
- **History.** Nothing new: tensile plastic strain goes into `eps_p` and `eps_p_eq` via
  the plastic corrector. The `apex` flag is set at the tension apex (hydrostatic
  `sigma_t`), the vertex of the combined surface; the Mohr-Coulomb apex lies outside the
  cut-off whenever `sigma_t < c cot(phi)`.
- **C++.** `MohrCoulomb::tensionReturn` copies the Python order and tests; 3x3 solves by
  Gaussian elimination with partial pivoting. `updateReturn` reports four new returns:
  tension plane, tension edge, tension apex, tension corner (mixed).

## Consequences

- Cohesive soils now fail in tension at `min(sigma_t, 2 c cos(phi)/(1 + sin(phi)))` in
  uniaxial tension, and hydrostatic tension is capped at `sigma_t`.
- With non-associated flow the first admissible set in the fixed order is taken; the
  Kuhn-Tucker solution is not proven unique there. No failure was found in 40,000
  random tension-biased increments; if no set is admissible the model throws.
- The vectorised 2D prototype (`prototype/constitutive/batch.py`) has no cut-off; it is
  only used by the WP 1.3 prototypes.
- No softening after tensile failure and no void opening: WP 3.2 and later.
