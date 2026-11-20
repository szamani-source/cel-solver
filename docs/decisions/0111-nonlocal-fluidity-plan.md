# 0111: Mesh dependence cure: nonlocal granular fluidity (plan, packages N.1 to N.4)

Date: 2026-11-20. Status: accepted (plan).

## Context

Sia (2026-11-20 16:04 UTC): before the Linux install the solver needs (1) a
friction law with a length scale so that shear bands and rates stop depending
on the mesh, (2) cohesion, (3) compaction. Decision 0038 recorded the
mesh-dependent collapse band (1 to 4 cells wide on every mesh, run-out 0.30 /
0.33 / 0.43 on 10 / 5 / 2.5 mm) and listed the nonlocal (gradient) fluidity
of Kamrin and Koval (2012) as a cure with a measurable length (the grain size
d). Sia asked for the rule with a length scale in it, so this is the cure.

## Decision

Model: the nonlocal granular fluidity (NGF) model of Kamrin and Koval (2012)
and Henann and Kamrin (2013), in the elasto-viscoplastic form Henann and
Kamrin (2014) used in Abaqus/Explicit:

- fluidity g [1/s] with mu = sqrt(J2)/P and the plastic shear rate
  gamma_dot_p = g mu;
- local fluidity g_loc = sqrt(P/(rhoS d^2)) (mu - mu_s)/(b mu) for
  mu > mu_s, else 0, with b = (mu_2 - mu_s)/I0 (the steady homogeneous
  limit is the linear local law mu = mu_s + b I, the small-I slope of the
  solver's mu(I) entries);
- nonlocal equation, solved each step on the Eulerian mesh for the material
  cells (void and walls: zero normal gradient):
  g - xi^2 lap(g) = g_loc, xi = A d / sqrt(|mu - mu_s|) (capped);
- stress update per point: elastic predictor, then the deviatoric stress
  relaxed along its own direction at fixed P with gamma_dot_p = g mu, where
  g = g_loc(mu_new) + (g_field - g_loc(mu_old)): the local part implicit
  (one scalar root), the nonlocal part lagged one step. A lagged g alone is
  unstable once G dt mu g'(mu)/P > 2, which explicit steps reach at once.
  No volume change in the flow (dilatancy with the nonlocal law is later).
  Tension (P <= 0) falls back to the base model's return.

Entries: `nonlocal yes; nonlocalA 0.48;` in a mohrCoulomb or druckerPrager
model with the mu(I) set (mu2, I0, d, rhoS). Off by default: existing cases
unchanged bit for bit. Glass-bead constants from Henann and Kamrin (2013):
mu_s 0.3819, mu_2 0.6435, I0 0.278, A 0.48 (to be checked against the paper
in N.1, as for every copied value).

Sia (18:19 UTC): develop and test on simple cases, not the silo (the silo comes back only after the law is proven).

Packages (one done-when test each, written first):

| WP | Title | Done when |
| --- | --- | --- |
| N.1 | Fluidity field and the point update | `tests/test_nonlocal.sh`: (a) homogeneous simple shear gives mu(I) to 1e-6; (b) the steady fluidity across a 1D layer with a stress gradient matches an independent Python solution of g - xi^2 g'' = g_loc to 1 % |
| N.2 | Coupled into the Eulerian step: sheared layer | A granular layer sheared between a moving lid and a rough base, under gravity, on 3 meshes (cells 4, 2, 1 d): the shear band width converges (in d) and the velocity profile matches the 1D nonlocal solution within 5 % |
| N.3 | Mesh convergence of the column collapse | Collapse a = 0.5 on 10 / 5 / 2.5 mm with the T0 set plus nonlocal: band width and run-out converge (5 and 2.5 mm within 5 %; now 0.30 / 0.33 / 0.43, decision 0038) |

## Consequences

- d becomes a physical length in the results; cells must still resolve the
  band (Sia: 10 to 20 d), so very coarse meshes get no benefit.
- One Helmholtz solve per step (OpenFOAM's linear solvers, host only).
- Cohesion and compaction follow as their own package series.

## N.1 result (2026-11-20)

`tests/test_nonlocal.sh` (under a minute) passes: kernel homogeneous shear
(MC and DP, dt 1e-6 and 1e-4 s, rates 1 to 50 1/s) gives mu_s + b I exactly;
the fluidity across a 40-cell layer with mu from 0.30 to 0.50 matches a
4000-cell independent solution within 0.1 % of g_max, with and without void
on top, and spreads fluidity into the part below mu_s (12 to 19 % of g_max).
Code: kernel `celNonlocal` (`nonlocalA` entry in mohrCoulomb and
druckerPrager), `nonlocalFluidity()` in celEulerian (OpenFOAM PCG, Gauss
linear corrected), field `g` written, driver `celFluidityCheck`.
Adding virtuals to the kernel Model changes its layout: rebuild every
application (./Allwmake) after pulling.

## N.2 first runs (2026-11-20 evening): the face shear jump dominates

`tests/nonlocal/shear_layer.py` (layer of 40 d, lid at 0.02 m/s, P0 200 Pa
under the lid, gravity, started at rest; 10 cells of 4 d, 2 s, steady):

| Run | Shear stress across the layer | Flow |
| --- | --- | --- |
| local mu(I) | 87 Pa in the body, 44 Pa in the top cell | top two cells only |
| nonlocal | 70 to 75 Pa, 32 Pa in the top cell | top three cells |
| local + flowImpedance | 60 Pa, top cells 50 to 69 | the lid slips over the top cell |
| nonlocal + flowImpedance | 46 Pa | the lid slips over the top cell |

In steady simple shear the shear stress must be uniform; it is not,
because the faces carry the Riemann shear jump z_s (v_R - v_L)/2 on top of
the cell stresses (z_s = sqrt(G rho), about 24 kPa s/m here: 240 Pa for a
0.01 m/s jump). That numerical viscosity (decision 0039: about
0.06 rho c_s h) is larger than the granular stress, so the band is set by
the face rule and the cell size, not by the law. Next (N.2): between
material cells with fluidity, take the face shear from the nonlocal
viscosity eta = P/g (the tangent of gamma_dot = g mu) in place of z_s, as
flowImpedance does for mu(I), and at fixed-velocity walls likewise; then
run 10 / 20 / 40 cells.
