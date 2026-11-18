# 0038 Collapse shear bands are mesh dependent; cure deferred (WP 3.6s)

Date: 2026-11-17. Status: accepted.

## Context

Decision 0037 put the plain Jop mu(I) law into the returns for T0 (phi = 35 deg) because
the Barker-Gray regularisation does not exist there, and asked whether shear bands become
mesh dependent. WP 3.6s ran the C++ column collapse (a = 0.5, E = 1 MPa) on 10, 5 and 2.5 mm
cells with plain Mohr-Coulomb and with mu(I) (`docs/wp/WP-3.6s.md`).

## Decision

1. Recorded as a known issue: the band is 1 to 4 cells wide on every mesh, its peak shear
   rate grows with refinement, and run-out does not converge (0.30, 0.33, 0.43 at 1 s).
   mu(I) with T0's placeholder constants changes nothing at 10 and 5 mm (I about 0.05).
2. No cure in WP 3.6s: each candidate is a package of its own. In order of preference,
   for a later package (after WP 3.6b to 3.6d, so the silo rate work is not held up):
   a. compressible mu(I)/phi(I) (dilatancy, survey section 2), which Barker et al. (2017)
      and Heyman et al. (2017) are cited as well posed (decision 0037; not checked);
   b. a non-local (gradient) fluidity with a grain-size length scale (Kamrin and Koval
      2012), measurable through d;
   c. a viscous regularisation scaled to a fixed physical length, the cheapest, but its
      length is not a lab-measured parameter (Sia's rule, decision 0036).
3. Mesh studies of collapse and discharge report results on at least two meshes and say
   which quantities converge (here: top height and the deposit body do; run-out and band
   width do not).
4. `inertialNumber` treats P below 1e-300 rho_s as P = 0 (I infinite), to avoid an
   overflow trap at free-surface cells. Results for any P above that are unchanged.

## Consequences

- The 10 mm collapse numbers in `tests/test_collapse.sh` and the 3.6e validation against
  Lube and Lajeunesse are mesh-bound until one of 2a to 2c is in.
- The silo discharge rate (WP 3.6d) may also depend on the band at the stagnant-zone edge;
  WP 3.6d should run two meshes, as WP 3.4b did.
- `tests/manual/test_shear_band.sh` repeats the study (about 9 min); the cure package
  turns its band width into a checked number.
