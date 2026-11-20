# 0108: Silo at rest by switching gravity on, not by filling from the top

## Context

The Couto silo at rest reads GA0 (top of the hopper) 25-30 % above the
experiment. Couto filled the silo from the top with the outlet closed, while
our runs place the material in the silo and ramp gravity on (0.1 s, damping
100 1/s, 0.5 s). Sia (2026-11-20) asked for a fill through an inlet narrower
than the silo, outlet closed, then settle and read the pressures.

## What was done

- `controlDict dampingStart` (default 0): velocity damping only from this
  time, so material free-falls while the silo fills and settles afterwards.
- `tests/hopper/top_fill.py`: aligned mesh with the STL shell inside it
  (decision 0107), 0.40 m headroom, inlet disc r < 0.25 m (silo R 0.5 m) at
  3 m/s until the Couto fill mass is in (3.0 s), damping from 3.3 s, end 4.0 s.
  About 3.5 h on 4 cores.

## Result (MC, T0 rules)

| | top fill | gravity on | Couto |
|---|---|---|---|
| GA0 kPa | 16.8 | 12.4 | 8.6-9.6 |
| GA1 | 9.4 | 4.6 | 4.0-4.5 |
| GA5 | 8.0 | 2.8 | 2.8-3.2 |
| GA7 | 4.3 | 1.3 | 1.4-1.7 |
| K = sxx/syy (core, z 1.6 m) | 0.90 | 0.43 | ~0.4 implied |

No leak; pressures steady over the last 0.8 s; mass 1541 kg against 1503
(packed denser). An earlier full-width top inlet gave the same doubling, so
an inlet pressing on the column is not the cause.

## Decision

Keep the gravity-on start for at-rest states. The grain lands at about 7 m/s
and is squeezed from all sides; Mohr-Coulomb (and Drucker-Prager) have no
volumetric yield, so this near-isotropic stress (any K between active and
passive is admissible) stays locked in. Realistic filling needs a compaction
law (cap or density-dependent stiffness), a separate package if wanted.
Zheng and Yu also start from a geostatic state.

## Consequences

The GA0 at-rest gap stays open with other causes to look at. `dampingStart`
and `top_fill.py` remain for later filling studies.
