# Constitutive models for bulk solids

The solver's purpose is bulk and particulate solids (sand, granular material, clay), so
pressure-dependent plasticity is a core feature, not an add-on. Every model implements one
interface (`celConstitutive`, WP 2.2): given the strain increment, the old stress and the
history variables at one point, return the new stress and history. The Eulerian and
Lagrangian sides call the same code.

Each model is first written and verified in a small Python single-point driver (WP 1.4),
then ported to C++ with the same tests (WP 2.2, 2.5a, 2.5b).

| Model | Use | Parameters | Verification test | Package |
| --- | --- | --- | --- | --- |
| Linear elastic | Baseline | E, nu | Uniaxial strain against closed form | 1.4, 2.2 |
| von Mises / Tresca (J2) | Undrained clay, metals | E, nu, yield stress or su | Uniaxial yield, radial return exactness | 1.4, 2.2 |
| Drucker-Prager | Sand and granular bulk, smooth cone | E, nu, friction angle, cohesion, dilation angle | Triaxial compression peak strength, non-associated flow volume change | 1.4, 2.5a |
| Mohr-Coulomb | Sand and granular bulk, standard soil parameters | E, nu, phi, c, psi, tension cut-off | Triaxial compression and extension strengths, corner and apex return | 1.4, 2.5b; optional tension cut-off sigma_t (Python and C++) 2.5c, decision 0015 |

## Points that matter for granular material

- **Non-associated flow.** Dilation angle psi is separate from friction angle phi;
  associated flow (psi = phi) over-predicts volume increase.
- **Apex and tension.** Cohesionless material cannot carry tension. The return mapping
  must handle the cone apex, and material that separates must drop to zero stress
  (links to void handling, WP 3.2). Mohr-Coulomb takes an optional tensile strength
  `sigma_t` (keyword `sigmaT`, WP 2.5c): no principal stress exceeds it, so cohesive
  material no longer carries tension up to c cot(phi). Drucker-Prager has no cut-off.
- **Corners.** Mohr-Coulomb has corners in the deviatoric plane; the return mapping in
  principal stress space handles them exactly. Drucker-Prager is matched to it
  (compression, extension or plane-strain fit, chosen per problem).
- **Large deformation.** Stress is carried with an objective rate (WP 2.3) and remapped
  with its history variables (WP 2.4, 3.4).

- **Rate-dependent friction (WP 3.6a, decision 0037).** Both frictional models take an
  optional mu(I) law (mu2, I0, d, rhoS; optional Barker-Gray muIAlpha): the return uses
  phi_m = atan(mu(I)) with I from the plastic shear rate and pressure. Without mu2 (or
  with mu2 = tan(phi)) the models are unchanged.

## Later candidates (not scheduled)

- Strain softening / critical state (state-dependent friction and dilation).
- Hypoplasticity for sand under cyclic or complex loading.
- Cap models for compaction.
