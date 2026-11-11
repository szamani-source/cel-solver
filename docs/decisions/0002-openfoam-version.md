# 0002 Pin OpenFOAM v2512 (openfoam.com)

Date: 2026-11-11. Status: accepted.

## Context

OpenFOAM has two main lines (openfoam.com by OpenCFD/ESI, openfoam.org by the
Foundation) with incompatible APIs. Supporting both would double maintenance.

## Decision

Pin openfoam.com v2512, the newest release on Docker Hub (`opencfd/openfoam-dev:2512`,
updated 2026-01-23). It ships isoAdvector geometric VOF
(`src/transportModels/geometricVoF/advectionSchemes/isoAdvection`), which WP 3.1 builds on.

## Consequences

- CI uses the `opencfd/openfoam-dev:2512` image.
- Moving to a newer release is a deliberate work package, not a side effect.
