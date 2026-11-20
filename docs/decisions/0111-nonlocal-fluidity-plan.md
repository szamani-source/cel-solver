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

- fluidity g [1/s] with mu = |tau|/P and the plastic shear rate
  gamma_dot_p = g mu;
- local fluidity g_loc = sqrt(P/(rhoS d^2)) (mu - mu_s)/(b mu) for
  mu > mu_s, else 0, with b = (mu_2 - mu_s)/I0 (the steady homogeneous
  limit is the Jop mu(I) law the solver already has);
- nonlocal equation, solved each step on the Eulerian mesh for the material
  cells (void and walls: zero normal gradient):
  g - xi^2 lap(g) = g_loc, xi = A d / sqrt(|mu - mu_s|) (capped);
- stress update per point with g fixed: elastic predictor, then the
  deviatoric stress relaxed along its own direction with
  gamma_dot_p = g mu (closed form for fixed g and P); the pressure part and
  dilatancy as in the current MC/DP return.

Entries: `nonlocal yes; nonlocalA 0.48;` in a mohrCoulomb or druckerPrager
model with the mu(I) set (mu2, I0, d, rhoS). Off by default: existing cases
unchanged bit for bit. Glass-bead constants from Henann and Kamrin (2013):
mu_s 0.3819, mu_2 0.6435, I0 0.278, A 0.48 (to be checked against the paper
in N.1, as for every copied value).

Packages (one done-when test each, written first):

| WP | Title | Done when |
| --- | --- | --- |
| N.1 | Fluidity field and the point update | `tests/test_nonlocal.sh`: (a) homogeneous simple shear gives mu(I) to 1e-6; (b) the steady fluidity across a 1D layer with a stress gradient matches an independent Python solution of g - xi^2 g'' = g_loc to 1 % |
| N.2 | Coupled into the Eulerian step | Steady flow down an incline (layer of 20 d): velocity profile within 5 % of the 1D solution, and a thin layer stops (h_stop) as the 1D model predicts within 10 % |
| N.3 | Mesh convergence of the collapse | Column collapse a = 0.5 on 10 / 5 / 2.5 mm with the T0 set plus nonlocal: band width in mm and run-out converge (5 and 2.5 mm within 5 %) |
| N.4 | Couto silo on two meshes | Rate and wall pressures on the aligned mesh at two resolutions within 10 % of each other; compared with the papers |

## Consequences

- d becomes a physical length in the results; cells must still resolve the
  band (Sia: 10 to 20 d), so very coarse meshes get no benefit.
- One Helmholtz solve per step (OpenFOAM's linear solvers, host only).
- Cohesion and compaction follow as their own package series.
