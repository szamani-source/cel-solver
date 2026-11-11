# 0001 Architecture: explicit FV Eulerian solver on OpenFOAM, own FE library

Date: 2026-11-11. Status: proposed (recommended in the plan doc; Sia asked to proceed).

## Context

Options considered: (A) modify a standard OpenFOAM flow solver; (B) a new explicit FV
solver using OpenFOAM only as infrastructure, plus our own FE library; (C) a standalone
C++ FE-based CEL code; (D) the Material Point Method instead of CEL.

## Decision

Option B. OpenFOAM provides mesh, fields, MPI, I/O and isoAdvector. Physics lives in our
libraries (`celConstitutive`, `celEulerian`, `celLagrangian`, `celCoupling`). Start with
Lagrange-plus-remap so results compare directly with Abaqus CEL; revisit in WP 1.3.

## Consequences

- OpenFOAM's standard solvers are not modified. Upgrades touch only thin adapters.
- FE/FV coupling is our own work, the main research risk.
- MPM (Anura3D) remains a useful validation reference, not the method.
