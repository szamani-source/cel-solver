# 0003 Licence: GPL-3.0

Date: 2026-11-11. Status: proposed, needs Sia's confirmation.

## Context

The solver links against OpenFOAM, which is GPL-3.0. Anything distributed that links
OpenFOAM must be GPL-compatible. Using it privately (not distributed) has no such duty.

## Decision

GPL-3.0 for the repository. This also lets us reuse code from solids4foam and
explicitSolidDynamics, which are GPL-3.0.

## Consequences

- If Sia wants a closed or commercial product, OpenFOAM cannot be the base when it is
  distributed; option C (standalone) would be needed. Confirm before WP 2.x.
