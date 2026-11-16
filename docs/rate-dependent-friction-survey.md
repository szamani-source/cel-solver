# Rate-dependent friction for the silo: literature survey (WP 3.6 planning)

Date: 2026-11-16. Asked by Sia after WP 3.4b (decision 0035): drop the soft numerical
stiffness (E = 0.1 MPa cannot be measured on the real material), go for rate-dependent
friction, and survey better approaches first. Recommendation: decision 0036.

Sources were found by title on 2026-11-16; the formulas and values below are as these
papers are widely cited and were **not re-read in the original text in this session**
(marked "not verified" where a number matters). Each WP that uses one re-reads it first.

## The problem in one paragraph

After rule 12 the C++ discharge rate still falls 16 % from E = 0.1 to 1 MPa at 5 mm
(decision 0035). The cause is numerical: the acoustic Riemann solver adds a viscosity of
order rho c h to plastic flow (decision 0028). A physical rate law adds its own
viscosity, but it only removes the E dependence if that viscosity is much larger than the
numerical one. A rough estimate for T0 (not measured, to be checked in WP 3.6b): the
mu(I) viscosity at small I is about (mu_2 - mu_s) d sqrt(rho_s p) / I0, which is
0.26 x 1e-3 m x sqrt(2650 x 5 kPa) / 0.28, about 3 Pa s near the outlet, while rho c h at
1 MPa and 5 mm is about 1800 x 24 m/s x 0.005 m, about 200 Pa s. **So mu(I) alone will
not make the rate independent of E.** The survey therefore covers both the rate law
(physics) and how other codes handle the elastic-to-flowing transition (numerics).

## 1. Local mu(I) rheology

- GDR MiDi (2004); Jop, Forterre and Pouliquen (2006, Nature 441): friction
  mu(I) = mu_s + (mu_2 - mu_s) / (I0 / I + 1), inertial number
  I = gamma_dot d / sqrt(p / rho_s). Glass beads: mu_s = tan 20.9 deg, mu_2 = tan 32.8,
  I0 = 0.279 (not verified).
- Parameters: grain diameter d and grain density rho_s (direct), mu_s, mu_2, I0 from
  standard lab tests: inclined-plane flows (h_stop curve and the Pouliquen flow rule
  Fr = beta h / h_stop), rotating drum, or annular shear / rheometer at controlled
  pressure (e.g. the controlled-pressure viscometric tests in JFM 2020, "Viscometric
  flow of dense granular materials under controlled pressure and shear stress").
  **All measurable.** mu_s is close to tan(phi) at critical state, so T0's phi = 35 deg
  becomes mu_s.
- Weakness: ill-posed (unbounded short-wave growth) at small and large I in the
  incompressible form (Barker, Schaeffer, Bohorquez and Gray 2015, JFM).
- Silo use: Staron, Lagrée and Popinet (2012, Phys. Fluids, "The granular silo as a
  continuum plastic flow: the hour-glass vs the clepsydra", arXiv 1211.5916) and Staron,
  Lagrée and Popinet (2014, EPJE, "Continuum simulation of the discharge of the granular
  silo: a validation test for the mu(I) visco-plastic flow law", arXiv 1306.6279):
  incompressible Navier-Stokes with mu(I) as a regularised viscosity (Gerris), Beverloo
  scaling recovered; they report that a constant friction already gives the D^1.5 law and
  that mu(I) mainly sets the prefactor and the velocity profiles (not verified). This is
  where silo reference constants 1.40 to 1.55 in decision 0019 came from.

## 2. Regularised mu(I)

- Barker and Gray (2017, JFM 828, "Partial regularisation of the incompressible mu(I)
  rheology for granular flow"): changes the low-I branch so the equations are well-posed
  for all I, keeping the experimentally fitted curve over the range where it was
  measured. Adds no parameter that cannot be fitted from the same tests.
- Heyman, Delannay, Tabuteau and Valance (2017, JFM, "Compressibility regularizes the
  mu(I)-rheology for dense granular flows") and Barker, Schaeffer, Shearer and Gray
  (2017, Proc. R. Soc. A, "Well-posed continuum equations for granular flow with
  compressibility and mu(I)-rheology"): adding a volume fraction law phi(I)
  (dilatancy) restores well-posedness. Schaeffer, Barker, Tsuji, Gremaud, Shearer and
  Gray (2019, JFM 874, "Constitutive relations for compressible granular flow in the
  inertial regime"): full compressible mu(I)/phi(I) model. phi(I) is measurable (volume
  fraction against I in the same shear tests), but it brings a second field law.
- Relevance: our model is compressible already (elastic pressure tied to density), so
  the low-I ill-posedness may be less severe than in the incompressible form, but a
  rate law with mu'(I) large at I -> 0 still needs the Barker-Gray branch to avoid
  mesh-dependent shear bands.

## 3. Nonlocal rheologies

- Kamrin and Koval (2012, PRL) and Henann and Kamrin (2013, PNAS, "A predictive,
  size-dependent continuum model for dense granular flows"): nonlocal granular fluidity
  (NGF). A fluidity field g obeys a screened Poisson equation with a cooperativity length
  xi = A d / sqrt(|mu - mu_s|); one extra parameter A (order 0.5 to 1), fitted from one
  flow geometry and then predictive in others (split-bottom cells, silos). Captures
  creep below mu_s and size effects (h_stop), which local mu(I) cannot.
- Cost: one more PDE per step (a scalar diffusion, explicit or a small implicit solve)
  and a parameter A that needs one calibration test beyond the mu(I) ones. Kamrin's later
  work extends it to the elastic regime and to MPM.
- Relevance: the a = 0.5 collapse creep (decision 0025, 0032) and the stagnant zones of
  the silo are where local mu(I) is weakest. A later package, not the first step.

## 4. The elastic-to-flowing transition in solid codes

- Dunatunga and Kamrin (2015, JFM, "Continuum modelling and simulation of granular
  flows through their many phases", arXiv 1411.5447): elasto-viscoplastic model in MPM:
  elastic below yield, mu(I) plastic flow above, a disconnected (stress-free) state below
  a critical density. Run on silos and column collapses. The structure is the same as
  ours (hypoelastic, Mohr-Coulomb/Drucker-Prager return, void). The rate law enters the
  return mapping: the plastic shear rate and the friction are solved together in one
  local (per point) equation. Their stiffness is a numerical choice as well (not
  verified), but MPM has no Riemann solver, so no rho c h face viscosity.
- Incompressible viscoplastic codes (Gerris/Basilisk: Lagrée, Staron, Popinet 2011
  column collapse; Chauchat and Médale 2014 granular rheology in an FV solver) have no
  elasticity at all: the regularised viscosity (capped at a large value) stands in for
  the solid regime. E then plays no role, but nor does elastic stress (Janssen at rest,
  arching), which T0 needs while it settles.
- Godunov solvers for elastoplastic solids (Miller and Colella 2001 and the many that
  followed): at a yielding face the waves are the plastic ones, with the tangent
  stiffness, not the elastic one. Using the elastic impedance in plastic flow is what
  gives our rho c h viscosity. The consistent fix is an impedance from the elastoplastic
  tangent at faces where both cells are at yield: much lower for the deviatoric part,
  the same for the volumetric part. Rules 8 and 12 are crude versions of this (cap at
  yield after the fact); an impedance from the tangent would be derived, not tuned, and
  has no new parameter.

## Summary

| Approach | Removes numerical E dependence? | New parameters | Measurable? | Cost |
| --- | --- | --- | --- | --- |
| Soft E (rejected by Sia) | Hides it | numerical E | No | none |
| Local mu(I) (Jop 2006) | No, adds ~3 Pa s against ~200 (estimate) | mu_s, mu_2, I0, d, rho_s | Yes (inclined plane, drum, shear cell) | small: rate law in the return map |
| Barker-Gray regularised mu(I) | Same as mu(I); well-posed | same (low-I branch from the same fit) | Yes | small |
| Compressible mu(I)/phi(I) | Same | + phi(I) law | Yes | medium |
| Nonlocal NGF (Henann-Kamrin) | No (physics); fixes creep and size effects | + A | Yes, one more calibration | medium: one PDE |
| Elastoplastic impedance at yielding faces | Yes, in principle (removes rho c h from plastic flow) | none | n/a | medium: Riemann solver |

## Recommendation (decision 0036)

1. Build mu(I) with the Barker-Gray low-I branch as the friction of the Mohr-Coulomb and
   Drucker-Prager models (rate-dependent return, Dunatunga-Kamrin style). It is the
   standard for silo continua, its parameters are lab-measurable, and the silo
   reference constants (decision 0019) were obtained with it.
2. In the same plan, replace rules 8 and 12 by an impedance from the elastoplastic
   tangent at faces where the material is flowing. mu(I) alone is not expected to remove
   the E dependence (estimate above); this is the numerical fix, without parameters.
3. Nonlocal fluidity (NGF) later, if creep or stagnant zones matter for Sia's cases.

The first package measures, on the simplest flow that has an exact answer, how large the
numerical viscosity is against the mu(I) one, so the estimate above is replaced by a
number before anything larger is built.
