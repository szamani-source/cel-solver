# 0070 Hex8 Lagrangian element: Flanagan-Belytschko with stiffness hourglass (WP 5.1)

Date: 2026-11-18. Status: accepted.

## Context

Phase 5 adds deformable Lagrangian bodies (decision 0001: explicit hex8 FE with hourglass
control). The element must pass the patch test on distorted meshes, stay stable under
hourglass modes, share the constitutive models of the Eulerian side, and be cheap.

## Decision

1. **Element.** One-point hexahedron with the uniform gradient of Flanagan and
   Belytschko (1981): B_I is the mean of grad N_I over the element, taken as the sum of
   the face integrals of N_I n dA, each exact by 2x2 Gauss points on the bilinear face.
   So sum_I x_I B_I = I and the element passes the patch test for any hexahedron
   (warped faces included). Internal force V sigma . B_I.
2. **Hourglass control.** Flanagan-Belytschko base vectors gamma (orthogonal to linear
   fields), stiffness form as LS-DYNA type 4: the hourglass forces are the time
   integral of k q, k = hourglass M V sum |B_I|^2 / 8, `hourglass` default 0.05. The
   stiffness form is conservative (energy kept to 1e-7 in the test); the viscous form
   would damp physical motion.
3. **Time integration.** Central difference with staggered velocities, forces at x_n,
   strain increment and spin from the velocity gradient on the mid-step configuration,
   Hughes-Winget rotation of the stress, then the celConstitutive update (with dt for
   rate-dependent models). Lumped mass. Stable step V / (largest face area) / c.
4. **Constraints and loads.** Per-node, per-direction prescribed velocities, external
   nodal forces, gravity, and `planeStrain yes` (no z velocity) for the 2D runs.
5. **Standalone first.** The body is built from nodes and elements; reading it from a
   case and coupling it to the Eulerian material are WP 5.3.

## Consequences

- `tests/test_hex8.sh` (seconds): patch test to 1e-15, 1D strain wave within 0.3 %,
  energy to 1.4e-7.
- Hourglass forces are not rotated with the element (fine for the small rotations per
  step of explicit runs; a body spinning through large angles would need it).
- One integration point: stress and plastic state are per element.
