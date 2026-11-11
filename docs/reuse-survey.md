# Reuse survey: existing open-source codes

Date: 2026-11-11. Pages checked on that date; licences as stated on each project page.

**Finding:** no open-source code provides CEL inside OpenFOAM, so there is nothing to adopt
whole. Several codes cover single pieces well and are worth reusing as code (GPL-compatible)
or as algorithm references (any licence).

| Code | What it is | Licence | How we use it |
| --- | --- | --- | --- |
| OpenFOAM v2512 isoAdvector | Geometric VOF interface advection, in OpenFOAM itself | GPL-3.0 | **Use directly** for volume fractions (WP 3.1) |
| [solids4foam](https://github.com/solids4foam/solids4foam) | Solid mechanics and FSI toolbox for OpenFOAM; v2 adds multi-material and solid-to-solid contact ([v2 notes](https://solids4foam.github.io/about/solids4foam-v2-vs-v1.html), [JOSS 2025](https://joss.theoj.org/papers/10.21105/joss.07407)) | GPL-3.0 | Reference and possible code source for constitutive laws and contact; Lagrangian, not Eulerian |
| [explicitSolidDynamics (vertex-centred)](https://github.com/jibranhaider/explicitsoliddynamics_vc) | Explicit FV for large-deformation solids in OpenFOAM: first-order hyperbolic laws for momentum and deformation gradient | GPL-3.0 | Closest FV reference for explicit solid dynamics in OpenFOAM; informs the WP 1.3 formulation choice |
| [Remhos](https://github.com/CEED/Remhos) | Monotone, conservative high-order remap for the Eulerian phase of ALE (MFEM) | BSD-2-Clause | Algorithm reference for the remap (WP 2.4) |
| [Fierro](https://github.com/lanl/Fierro) | LANL C++ code for transient material dynamics, Lagrangian and Eulerian, multi-material volume fractions | BSD-3-Clause | Algorithm reference for multi-material cells (WP 3.3) |
| [NimbleSM](https://github.com/sandialabs/nimblesm) | Sandia explicit Lagrangian FE with contact; archived June 2025 | Licence file present, type not read | Reference for explicit hex8 FE and contact (WP 5.1, 4.3) |
| [WeldFormFEM](https://github.com/mfkiwl/WeldFormFEM) | Explicit FEM solver for forming; lists a CEL solver based on Benson's work, contact still in progress | Not stated | Read-only reference; do not copy code (no licence) |
| [Anura3D](https://github.com/Anura3D/Anura3D_OpenSource) | MPM for large deformation and soil-water-structure interaction | LGPL-3.0 (licence file in repo) | Independent reference results for penetration benchmarks (WP 4.4) |

## A known weakness to target

[Staubach and Martinelli (2023), "MPM vs. CEL: Numerical modelling of penetration processes"](https://proceedings.open.tudelft.nl/seg23/article/view/600)
found that MPM (Anura3D) and Abaqus CEL agree on cone factors for undrained penetration,
but CEL tip stresses oscillate strongly and refinement barely helps. A smooth contact
force is therefore a concrete goal for WP 4.3 and 4.4, and a way our code can improve on
the commercial one.
