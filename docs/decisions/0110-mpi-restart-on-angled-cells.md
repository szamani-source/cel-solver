# 0110: MPI and restart on angled cells (WP G.6), and the wedge on the aligned mesh (WP G.4)

## Context

WP G.6 asks that 4-processor and restarted runs on angled meshes match serial
to round-off. WP G.4 asks that the Couto silo at rest on the wall-aligned
wedge mesh (fanned trapezoid blocks, free surface) gives Janssen within 5 %.

## Decision

- Coupled (processor, cyclic) faces now read the neighbour cell's mesh-line
  frame for the neighbour side of the free-surface contact test, as an
  internal face does (`dirsNbr` in the Lagrange step). Before, the owner's
  frame was used for both sides. Cartesian and single-processor runs are
  unchanged (the frame is only built on angled meshes).
- `tests/test_angled_mpi_restart.sh` (about 1 minute):
  - column collapse (G.2 case) on leaning cells, serial against 4
    processors split 2 x 2, without the shear cap: largest scaled
    difference 9e-13 (limit 1e-10);
  - the siloDischarge tank with leaning blocks, discharging, restarted at
    0.4 s against the continuous run to 0.5 s, serial and on 2 processors:
    fields, state and mass counters identical bit for bit.
- With the shear cap on, the leaning-cell collapse on 4 processors parts
  from serial at about 1e-8 (stress, scaled by rho g H). The capped-face
  count differs from step 6 (5 against 6 faces): the cap is a switch that
  flips on a round-off difference, which the plastic flow then carries
  (as for the Cartesian collapse later in its run, decision 0032). On
  leaning cells the processor faces already differ at 1e-16 at the first
  steps (face geometry is not exact in binary), so the flip comes earlier.
  Without the cap the same runs match at 1e-13; the Cartesian capped run
  matches at 1e-13 to 0.1 s. The test reports the capped difference and
  does not check it.
- `tests/test_angled_wedge.sh` (G.4, about 13 minutes): Couto silo at rest
  on the 25 mm aligned wedge mesh, T0 rules and fullBlend, 0.5 s. Mass
  drift 6e-15, closed outlet holds, largest speed 0.0018 m/s, interior
  stress ratio 0.425 against nu/(1-nu) 0.429, GA1 to GA7 within -2.6 to
  +0.8 % of Janssen. Hopper wall at rest 11.6 kPa (GA0 8.6 to 9.6, as on
  the cut-cell wedge, H.3).

## Consequences

- Parallel runs on angled meshes with the shear cap give the same physics
  but not identical numbers to serial; compare them by results (rates,
  pressures), not bit for bit.
- G.4 and G.6 are done. G.3 keeps its open items (top full row under the
  void on turned cells, sheared mesh).
