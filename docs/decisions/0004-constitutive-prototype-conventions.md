# 0004 Constitutive prototype conventions (WP 1.4)

Date: 2026-11-14. Status: accepted for the Python prototype; the C++ port (WP 2.2, 2.5a,
2.5b) should keep them unless a new note says otherwise.

## Context

The Python single-point driver is the reference the C++ library will be tested against,
so its conventions and algorithm choices must be fixed and written down.

## Decision

- Voigt order `xx, yy, zz, xy, yz, zx`; strain uses engineering shear (gamma = 2 eps),
  stress uses tensor components. Stress positive in tension, angles in degrees, SI units.
- One signature for every model: `(deps, sigma, history, params) -> (sigma, history)`.
  History is a dict (`eps_p`, `eps_p_eq`, `apex`) and is never changed in place.
- Elastic predictor, plastic corrector. Drucker-Prager and Mohr-Coulomb are perfectly
  plastic (no hardening or softening yet); J2 has linear isotropic hardening.
- Drucker-Prager: `f = sqrt(J2) + eta p - xi c`, potential with `eta_bar` from psi.
  Three matches to Mohr-Coulomb, chosen per problem by `fit`: `compression` (outer cone,
  default), `extension` (inner cone), `plane_strain`. Cone return, else apex return.
- Mohr-Coulomb: return mapping in sorted principal stress space. Main plane, then the
  edge the main return crosses (triaxial compression or extension corner), then the apex.
  Every return is a small linear solve with active planes, so it is exact for perfect
  plasticity; it does not copy textbook closed-form coefficients.
- Mixed stress/strain control in the driver uses Newton iteration with a finite-difference
  Jacobian, so models need no consistent tangent yet.

## Consequences

- No separate Rankine tension cut-off yet. Cohesionless material already returns to zero
  stress at the apex; a cut-off below c cot(phi) is left for WP 2.5b if a problem needs it.
- A consistent tangent is not needed by the explicit solver; implicit use would need one.
