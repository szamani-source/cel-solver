# How a rigid body or STL wall pushes on the Eulerian material (WP 4.0)

Date: 2026-11-17. Recommendation: decision 0042. Experiments: `prototype/contact2d.py`,
checked by `tests/test_proto_contact_schemes.sh`.

## Sources (opened 2026-11-17; text read through a fetch tool, so quotes are as extracted)

| Id | Source | What it gave |
| --- | --- | --- |
| S1 | H. Chen (Ansys/LST), "*ALE_STRUCTURED_FSI: the new S-ALE FSI solver", 16th Int. LS-DYNA Users Conference, 2020. https://www.dynalook.com/conferences/16th-international-ls-dyna-conference/fsi-ale-computing-technology-t8-1/t8-1-b-fsi-ale-081.pdf | `*CONSTRAINED_LAGRANGE_IN_SOLID` (CLIS) has 38 parameters over 22 manual pages and "leakage is not always prevented" (blast example). The new S-ALE coupling is penalty only, 9 parameters, only PFAC to tune; leakage control, normals, edge coupling and coupling points automatic. A compressed bottle leaked with CLIS, not with S-ALE |
| S2 | LSTC manual excerpt, `*ALE_STRUCTURED_FSI` (undated). https://ftp.lstc.com/anonymous/outgoing/hao/sale/manuals/ALE_STRUCTURED_FSI.pdf | One penalty method, "similar to a combination of CTYPE = 4 and 5". PFAC scales the estimated critical stiffness, or is a curve of coupling pressure against penetration. Coupling to all materials not recommended (no interface to track leakage) |
| S3 | Ansys Innovation Space forum, "ALE-Lagrange coupled simulation leakage issue", Oct 2023, answer by an Ansys engineer. https://innovationspace.ansys.com/forum/forums/topic/ale-lagrange-coupled-simulation-leakage-issue | CLIS practice: CTYPE 4 (penalty) in most cases; NQUAD at least 2 per ALE element, too many "gives stiff response and instability"; PFAC curve is the main lever (penetration ~0.1 ALE width); ILEAK 2 if leakage persists; Lagrangian to ALE mesh about 1:1; leakage often from bad meshes and sharp corners |
| S4 | Abaqus 2024 documentation, "Contact constraint enforcement methods in Abaqus/Explicit". https://docs.software.vt.edu/abaqusv2024/English/SIMACAEITNRefMap/simaitn-c-expcontactconstraints.htm | General contact always uses penalty enforcement |
| S5 | Abaqus 2023 documentation, "Eulerian analysis". https://docs.software.vt.edu/abaqusv2023/English/SIMACAEANLRefMap/simaanl-c-euleriananalysis.htm | Eulerian-Lagrangian contact is penalty, stiffness "automatically maximized subject to stability limits"; "an enhanced immersed boundary method"; material boundaries rebuilt from volume fractions with planar facets each increment; default "may exhibit leakage when modeling liquids or gases" (an enhanced, costlier option prevents it); penetration at corners of about one Eulerian element (fillet radius of one element advised); by default no tension and zero friction across the interface; "a simple regular grid of Eulerian elements often yields the best accuracy" |
| S6 | Abaqus 2024 documentation, "General contact". https://docs.software.vt.edu/abaqusv2024/English/SIMACAEITNRefMap/simaitn-c-contactgeneral.htm | E-L contact only on Lagrangian solid and shell surfaces; compensates mesh size differences |
| S7 | Monasse, Daru, Mariotti, Piperno, Tenaud, "A conservative coupling algorithm between a compressible flow and a rigid body using an embedded boundary method", arXiv 1012.5479 (2010). https://arxiv.org/abs/1012.5479 | Cut-cell coupling with exact mass conservation and momentum and energy balance claimed. Method sections not read (not verified) |
| S8 | Arienti, Hung, Morano, Shepherd, "A level set approach to Eulerian-Lagrangian coupling", J. Comput. Phys. 185 (2003) 213-251. https://shepherd.caltech.edu/EDL/publications/reprints/vtf_jcp.pdf | Ghost-fluid coupling (no cut cells): no small-cell limit, but conservation not guaranteed; force on the body = boundary integral of the traction |
| S9 | Natarajan, Grout, Zhang, Day, "A moving embedded boundary approach for the compressible Navier-Stokes equations in a block-structured adaptive refinement framework", arXiv 2108.00126 (2021). https://arxiv.org/abs/2108.00126 | Unsplit cut-cell Godunov scheme: wall flux from the Riemann solver with the body's normal velocity; small cells by hybrid update plus conservative redistribution; uniform motion preserved exactly |
| R1 | Staubach and Martinelli 2023 (in `docs/reuse-survey.md`, read earlier) | Abaqus CEL tip stresses oscillate strongly in penetration |

Not found or not read: the CLIS constraint coupling (CTYPE 1 to 3) in a primary source, CLIS
damping and friction options, a published comparison of LS-DYNA and Abaqus coupling on the
same problem.

## The schemes

| | Penalty (LS-DYNA CLIS/S-ALE, Abaqus CEL) | Cut-cell moving-wall Riemann (embedded boundary) | Ghost fluid / immersed forcing |
| --- | --- | --- | --- |
| Leakage | Depends on stiffness; both vendors document it and give options against it (S1, S3, S5) | None through the wall by construction (no flux through covered face parts) | Possible; not conservative (S8) |
| Conservation | Momentum yes (equal and opposite forces), energy no | Mass exact, momentum exact with the body force | No |
| Force noise | Grows with stiffness (R1; our experiments) | Small; wall traction from the same Riemann solver as interior faces | Moderate (not measured) |
| Friction | Regularised penalty or damper on slip | Coulomb cap on the face shear, as the box walls (decision 0030) | Hard |
| Time step | Stiffness must stay under its own limit (S3: instability) | Small cut cells need mixing or redistribution (S9) | No limit |
| Body edge in a cell | Overlap measured from volume fractions | Exact cut geometry: face apertures and the wall segment | Level set |
| Fit with our scheme | Adds a force term; the Riemann solver never sees the body | The wall is one more face with prescribed velocity, which our Lagrange step already has (walls, outlets, prescribed motion) | Needs ghost states for a cell-centred scheme |
| Parallel | Overlap is cell-local | Cut geometry cell-local; mixing reaches one neighbour (may cross a processor) | Ghost layer |

## Experiments (prototype/contact2d.py, 5 mm cells, E = 10 MPa, nu = 0.3, rho = 1800)

First-order unsplit Lagrange-plus-remap, linear elastic with void, the same for both
couplings except at the body. Penalty stiffness k = PFAC M / h (M the P-wave modulus),
scaled by each cell's material fraction, with an impedance damper on approach.

1. **Piston**: full-width rigid plate pushed into a 0.2 m layer at 2 m/s. Exact contact
   stress z v until the wave returns. Plate edge at 0, 0.25, 0.5 of a cell.
2. **Translating block**: 4 x 4 cell rigid block moving 6 cells at 2 m/s through material
   at rest in a closed box; edges on and off cell faces.
3. **Inclined wall**: elastic block (0.12 x 0.05 m) on a fixed wall at 30 deg to the mesh
   with friction delta_w = 20 deg, under gravity. Exact sliding acceleration
   g (sin 30 - tan 20 cos 30) = 1.81 m/s^2.

Leak = material beyond the open volume of cut and covered cells (overfill above 2 % elastic
compression), as a fraction of all material. Noise = standard deviation of the contact
stress over z v (piston), mean step-to-step jump of the drag over z v A (block).

| Scheme | Piston stress / zv (offsets 0, 0.25, 0.5) | Piston noise | Piston leak | Block noise (0, 0.5) | Block drag ratio between offsets | Block leak | Incline acc. / exact | Incline leak |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cut-cell Riemann | 1.006, 1.006, 1.006 | 0.002 | 0.7e-4 to 5e-4 | 0.064, 0.017 | 0.93 | 3e-4, 5e-4 | 1.003 | 4e-6 |
| Penalty PFAC 0.1 | 0.84, 0.88, 0.86 | 0.14 to 0.21 | 1.1e-2 | 0.044, 0.065 | 0.83 | 5e-3 | 0.97 | 1.9e-3 |
| Penalty PFAC 1 | 0.99, 0.99, 1.00 | 0.27 | 3.8e-3 | 0.32, 0.13 | 7.7 | 1.6e-3 | 1.11 | 5.7e-4 |
| Penalty PFAC 10 | 0.97, 0.96, 1.14 | 1.1 to 1.2 | 1.8e-3 | 5.0, 1.0 | 3.2 | 8e-4 | 1.95 | 1.8e-4 |

Mass is conserved to 1e-15 in all runs (after limiting each donor's outflow to its mass;
before, small cut cells lost 1e-4). Riemann leak values are compressed material in small
cut cells, not material inside the body: no mass crosses a covered face. Penalty shows
the textbook trade: soft springs leak and under-predict the force, stiff springs ring
(noise 27 % at PFAC 1, over 100 % at 10) and get friction wrong (the stiff normal force
oscillates, so the Coulomb limit mu |Fn| is applied to a noisy Fn).

## Is LS-DYNA's coupling better than Abaqus's?

Both are penalty couplings (S2, S3, S5). In what we could read, LS-DYNA gives the user
more control (PFAC curves, coupling point density, leakage control options, and the newer
S-ALE coupling that automates them), and its documentation addresses leakage head on;
Abaqus picks the stiffness automatically and documents leakage for fluids and about one
element of penetration at corners. So "better" holds in the sense of more control and
an actively developed leakage treatment, not of a different principle; no side-by-side
published comparison was found (not verified either way). For us the point is that a
cell-centred Godunov code can avoid the penalty trade altogether.

## Follow-up: LS-DYNA-style coupling points (Sia, 2026-11-17)

Sia noted that LS-DYNA handles contact points better than Abaqus and asked us to use
its method or a similar one. LS-DYNA's coupling differs from Abaqus's in three ways
we could act on: coupling points spread over each Lagrangian segment (NQUAD, at least
two per ALE element, S3), penetration tracked at each point rather than from the
volume fractions of a cell, and controls on top of the spring (PFAC curve, damping,
leakage control ILEAK, S1 to S3). The keyword card and theory text for
`*CONSTRAINED_LAGRANGE_IN_SOLID` could not be read online (the R14 theory manual and
the 960 keyword and 2005 theory manuals returned only their contents pages), so the
details below are our reconstruction, **not verified** against the manual.

`prototype/contact2d.py` scheme `points`: points every h / NQUAD along the body
surface; at each point the material velocity is interpolated from the cells holding
material, penetration d grows by the relative normal velocity into the body (never
below 0, compression only), pressure k d (or a pressure-penetration curve) plus a
damper DAMP z (relative normal velocity), Coulomb friction; the force goes to the
four nearest cells by weight times mass, the reaction to the body. Leak control
(`ileak`) adds the cell-overlap penalty on top.

| Variant (k in M/h) | Piston stress / zv | Piston noise | Piston leak | Block noise | Incline acc. / exact | Incline leak |
| --- | --- | --- | --- | --- | --- | --- |
| points, k 0.1 | 0.48, 0.46 | 0.03, 0.07 | 2.1e-2 | 0.007 | 1.03 | 3.9e-2 |
| points, k 1 | 0.46, 0.48 | 0.06, 0.09 | 2.0e-2 | 0.015 | 1.03 | 3.9e-2 |
| points, k 1, DAMP 0.5 | 0.46, 0.47 | 0.05, 0.08 | 2.0e-2 | 0.014 | 1.03 | 3.9e-2 |
| points, k 10, DAMP 0.5 | 0.46, 0.47 | 0.20, 0.13 | 2.0e-2 | 0.02 | 1.03 | 4.0e-2 |
| points, k 1, NQUAD 4 | 0.46, 0.48 | 0.06, 0.09 | 2.0e-2 | 0.013 | 1.03 | 3.9e-2 |
| points, k 1, leak control | 0.98, 1.00 | 0.27, 0.26 | 3.8e-3 | 0.30, 0.16 | 1.11 | 5.7e-4 |
| cut-cell Riemann (decision 0042) | 1.006, 1.006 | 0.002 | 5e-4 | 0.06, 0.02 | 1.003 | 4e-6 |

Coupling points give a smooth force (noise 3 to 9 % instead of 27 %) and friction close
to right, which is the improvement over a cell-overlap spring that LS-DYNA users report.
But the material in a cell the surface cuts is spread over the whole cell, so it sits
partly inside the body (leak 2 to 4 %) and the plate pushes on it at about half the
right stress; more points or damping change nothing. Leak control fixes the leak and
the force but brings back the noise of the plain spring. In this test the point idea
and the cut-cell wall agree on what matters: the force must act at the true surface
position, and the material must not be spread across it. The cut-cell wall does both
exactly, with no stiffness to tune.
