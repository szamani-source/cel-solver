# Validation references (journal papers with data)

Status: list agreed in kind (Sia asked for "a couple of good journals" for later
validation, 2026-11-17). Not yet used by any test. Scalar laws already fixed for the
tests (Janssen, Beverloo, collapse run-out prefactors) live in `silo-reference.md`;
this list is for the **field data** we will compare against next: velocity profiles,
free-surface shapes over time, wall stresses, and runs of a continuum mu(I) code on
the same geometry.

DOI check: Crossref and doi.org are blocked from the cloud container. "Search" means
the DOI was confirmed by a web search on 2026-11-17; "memory" means it was not
resolved and must be checked before a value is copied from the paper.

## The papers

| # | Paper | Journal | DOI / open copy | DOI check | Data it gives | Serves |
| --- | --- | --- | --- | --- | --- | --- |
| V1 | L. Staron, P.-Y. Lagree, S. Popinet, "Continuum simulation of the discharge of the granular silo: a validation test for the mu(I) visco-plastic flow law" (2014) | Eur. Phys. J. E 37, 5 | 10.1140/epje/i2014-14005-6; arXiv 1306.6279 | Search | 2D flat-bottomed silo, continuum mu(I) against DEM on the same geometry: discharge rate against outlet width (Beverloo form), velocity and pressure fields, rate independent of the head | T0 discharge, WP 3.6 (mu(I)) |
| V2 | Y. Zhou, P.-Y. Lagree, S. Popinet, P. Ruyer, P. Aussillous, "Experiments on, and discrete and continuum simulations of, the discharge of granular media from silos with a lateral orifice" (2017) | J. Fluid Mech. 829, 459 | [Cambridge](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/experiments-on-and-discrete-and-continuum-simulations-of-the-discharge-of-granular-media-from-silos-with-a-lateral-orifice/0B924467BC7DF119B7E03B190EDDA999) | Memory (10.1017/jfm.2017.643) | Quasi-2D **experiments** plus DEM plus continuum mu(I): measured rate against orifice size, velocity profiles at the orifice, free-surface shape | T0 (side-outlet variant), the only one with experiments and a continuum code side by side |
| V3 | J. Choi, A. Kudrolli, M. Z. Bazant, "Velocity profile of granular flows inside silos and hoppers" (2005) | J. Phys.: Condens. Matter 17, S2533 | arXiv cond-mat/0501568; [ADS](https://ui.adsabs.harvard.edu/abs/2005JPCM...17S2533C/abstract) | Memory (10.1088/0953-8984/17/24/008) | Measured velocity fields in a quasi-2D silo and hopper (glass beads, imaging): vertical velocity against x at several heights, width of the flowing channel growing with height | T0 flow pattern (channel, stagnant corners); hopper variant later |
| V4 | G. Lube, H. E. Huppert, R. S. J. Sparks, A. Freundt, "Collapses of two-dimensional granular columns" (2005) | Phys. Rev. E 72, 041301 | 10.1103/PhysRevE.72.041301 | Search | Quasi-2D **experiments**, several materials: run-out and final height against aspect ratio, profiles over time, failure surface | Column collapse (WP 3.2 tests) |
| V5 | P.-Y. Lagree, L. Staron, S. Popinet, "The granular column collapse as a continuum: validity of a two-dimensional Navier-Stokes model with a mu(I)-rheology" (2011) | J. Fluid Mech. 686, 378 | Memory (10.1017/jfm.2011.335) | Memory | Continuum mu(I) against DEM for 2D collapse: free-surface profiles at fixed times, run-out against aspect ratio, deposit shape; same group and law as V1 | Column collapse with mu(I) (WP 3.6) |
| V6 | GDR MiDi, "On dense granular flows" (2004) | Eur. Phys. J. E 14, 341 | Memory (10.1140/epje/i2003-10153-0) | Memory | Compiled experiments and DEM in six geometries; inclined plane: h_stop(theta), Pouliquen flow rule, Bagnold velocity profile; source of the mu(I) form | WP 3.6b inclined layer, mu(I) constants |
| V7 | N. J. Balmforth, R. R. Kerswell, "Granular collapse in two dimensions" (2005) | J. Fluid Mech. 538, 399 | [author PDF](https://secure.math.ubc.ca/~njb/Research/gran.pdf) | Search (journal and title; added in WP 3.6s, numbers not read) | Quasi-2D **experiments** in narrow and wide channels with several grains: run-out and final height against aspect ratio, sensitivity to channel width and grain friction | Column collapse (WP 3.6e), with V4 |
| V8 | Q.J. Zheng, A.B. Yu, "Finite element investigation of the flow and stress patterns in conical hopper during discharge" (2015) | Chem. Eng. Sci. 129, 49 | 10.1016/j.ces.2015.02.022 | Paper (Sia's copy) | Abaqus Eulerian Mohr-Coulomb hopper: outlet velocity, flow mode, flat-bottom Beverloo C 0.59, wall pressure | Conical hopper, `conical-hopper.md` |
| V9 | A. Couto, A. Ruiz, P.J. Aguado, "Experimental study of the pressures exerted by wheat stored in slender cylindrical silos..." (2013) | Powder Technol. 237, 450 | 10.1016/j.powtec.2012.12.030 | Paper (Sia's copy) | **Experiments**: wheat silo wall pressures at rest and in discharge, 7 tests | Conical hopper, `conical-hopper.md` |

For the later targets (T1 footing, T2 pile), the first cross-check is
S. Qiu, S. Henke, J. Grabe (2011), "Applying the coupled Eulerian-Lagrangian method
to geotechnical problems", Computers and Geotechnics 38, 30 (Abaqus CEL; DOI from
memory 10.1016/j.compgeo.2010.09.002, not resolved).

## Journals to watch

Journal of Fluid Mechanics, Physical Review E / Letters, European Physical Journal E,
Granular Matter (already the source of S1 and S6 in `silo-reference.md`), Powder
Technology (silo experiments), Computers and Geotechnics (CEL in geotechnics).

## How these will be used

Each comparison becomes its own work package with a "done when" test, after the
mu(I) packages (WP 3.6b to 3.6e). Order: V1 (closest to T0 and to our model), then
V3 and V2 for velocity fields, then V4 and V5 for collapse. Before a number is copied
into a test, open the paper, record the page or figure, and resolve the DOI.
