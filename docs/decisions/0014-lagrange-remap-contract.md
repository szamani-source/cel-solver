# 0014 Contract between the Lagrange step and the remap

Date: 2026-11-14. Status: accepted (written with the WP 2.3 and 2.4 briefs so the two
packages can be built in parallel).

## Context

Decision 0013 fixes the scheme: unsplit cell-centred Lagrange step, then one remap.
WP 2.3 fills `eulerianSolid::lagrangeStep(dt)` and WP 2.4 fills `eulerianSolid::remap(dt)`.
They are built by different sessions at the same time, so what one hands the other must
be fixed in advance.

## Decision

- `eulerianSolid` has a face field `phiL_` [m^3/s]: the face velocity of the Lagrange
  step dotted with the face area vector, positive out of the owner cell (OpenFOAM flux
  convention), on internal and boundary faces.
- After `lagrangeStep(dt)`: the cell fields describe the moved (Lagrangian) cells.
  Cell mass is unchanged; the cell volume is `V_L = V + dt * sum(phiL out of the cell)`;
  `rho_` is the Lagrangian density `m / V_L`; `U_`, `sigma_`, `epsP_`, `epsPEq_` are the
  updated values of that material. The mesh itself never moves.
- `remap(dt)`: moves the volume `dt * phiL` back across each face, taking material from
  the cell the flux leaves (owner if `phiL > 0`), and returns every field to the fixed
  mesh: density by volume flux, velocity, stress and history mass-weighted, limiters per
  field as decision 0013 item 5. After it, `rho_` is mass over the fixed cell volume.
- With `remap` empty (before WP 2.4 merges) the state is consistent only to first order
  in the strain per step; WP 2.3's own tests therefore use small amplitudes or
  volume-preserving prescribed motion.
- Each package keeps its code in its own source file
  (`eulerianSolidLagrange.C`, `eulerianSolidRemap.C`) to avoid merge conflicts.
- Test programs may set a state by hand through the `...Ref()` accessors.

## Consequences

- Later packages reuse `phiL`: volume fractions (WP 3.1, isoAdvector) and contact.
- A corner-transport or swept-region remap would need more than a face flux; add it to
  this contract only if WP 2.4 or the integration checks (WP 2.4b) show the need.
