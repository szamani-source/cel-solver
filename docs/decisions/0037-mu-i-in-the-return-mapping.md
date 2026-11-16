# 0037 mu(I) friction in the Drucker-Prager and Mohr-Coulomb returns (WP 3.6a)

Date: 2026-11-16. Status: accepted (follows decision 0036, chosen by Sia).

## Context

Decision 0036 plans the Jop, Forterre and Pouliquen mu(I) friction with the Barker and
Gray (2017) low-I regularisation inside the existing returns. Sia requires every
parameter to be measurable. Both papers were re-read for this package (Jop et al. 2005,
JFM 541, eq. 4.3 and section 4.1; Barker and Gray 2017, JFM 828, eqs. 3.9, 5.2, 5.3).

## Decision

1. **How mu enters.** mu sets a mobilised friction angle phi_m = atan(mu), used in
   place of phi by the unchanged Drucker-Prager (same fit) or Mohr-Coulomb return; psi,
   c and the tension cut-off are kept (sigmaT is lowered to c cot(phi_m) if it would lie
   above that apex). So mu_s = tan(phi): T0's phi = 35 deg gives mu_s = 0.700.
2. **Rate.** I = dgamma/dt d sqrt(rho_s/P), with dgamma = sqrt(2 e:e) of the deviatoric
   plastic strain of the return (tensor components; gamma_xy in simple shear) and P the
   mean pressure after the return. P <= 0 with flow gives I = infinity, mu = mu_2.
3. **Solve.** phi_m is the root of R(mu) = mu - mu(I(mu)), which increases with mu, so
   the root is unique and bracketed (by mu_s and mu_2, or below mu_s with the
   regularisation). Illinois method to 1e-14 mu_2 (several base returns per flowing
   point; not counted). The kernel gets `Model::updateRate(deps, dt, ...)`; models without
   mu(I) call `update`, so their results are unchanged bit for bit, and so is a model
   with mu2 = tan(phi).
4. **Regularisation is optional, and off for T0.** The Barker-Gray branch
   mu = sqrt(alpha/ln(A/I)) replaces the Jop law below the first neutral point I_1^N of
   their well-posedness function C(I). For mu_s = tan 35 deg the Jop law has **no**
   well-posed range at all (C > 0 for every I, tried mu_2 - mu_s = 0.26 to 0.8), so the
   branch cannot be built; it exists for low-friction materials (Jop's glass beads:
   I_1^N = 0.0194). The code builds it when `muIAlpha` is given (Barker and Gray use
   1.9; it must be below 2) and refuses it when there is no well-posed range. T0 runs
   the plain Jop law: elastic below mu_s, as the present model.
5. **Placeholders for T0** (not Sia's material): mu_s = tan 35 deg; mu_2 - mu_s = 0.2615
   and I0 = 0.279 from Jop et al. (2005) glass beads (tan 32.76 - tan 20.9, d = 0.53 mm,
   rho_s = 2450); d = 1 mm and rho_s = 2650 kg/m^3 (quartz). No sand mu(I) values were
   found and verified in this session.

## Consequences

- Parameters: mu2, I0, d, rhoS (all measurable: grain size and density directly; mu_2
  and I0 from inclined-plane flow rules, drum or controlled-pressure shear tests).
  alpha is not a material parameter; it is only used where the regularisation exists.
- The incompressible ill-posedness that Barker and Gray cure is not cured for T0. Our
  model is compressible and elastic, which may regularise it (Heyman et al. 2017,
  Barker et al. 2017; not checked). WP 3.6b should look for mesh-dependent shear bands;
  if they appear, the compressible mu(I)/phi(I) route (survey section 2) is next.
- With the regularisation there is no yield stress (mu -> 0 as I -> 0): every stressed
  point creeps, negligibly far below mu_s.
- Mohr-Coulomb in steady isochoric simple shear gives tau/sigma_n = sin(phi_m), not
  tan(phi_m) (the usual Davis result); Drucker-Prager gives the fit's eta(phi_m). Lab
  values of mu measured as tau/sigma_n on the shear plane therefore map to phi_m by the
  chosen model, as phi does today. The WP 3.6a test reads phi_m back from the stress.
- Cost: each flowing point runs several returns per step. Not measured in celFoam yet
  (WP 3.6b).
