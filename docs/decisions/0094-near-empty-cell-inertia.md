# 0094: Near-empty cells get a floor on their inertia; the sliver fix is merged

Date: 2026-11-19 (accuracy package A.3b, thread "Fix body contact above surface")
Status: accepted.

## Context

Decision 0093 fixed the geometry bug behind the steel packing, but the parked sliver
fix of decision 0092 still could not merge: in the full suite `test_wedge_silo.sh`
crashed (3D quarter, STL tube) and the steel push in `test_deformable_push.sh` spiked to
1.032 rho0, each on round-off-level changes of the mixing arithmetic.

## Cause of the silo crash

At t = 3.5 s, with the tank nearly empty, cells at the tube wall near the top held a
round-off residue of material whose mass and fill disagree: 1e-18 kg at alpha 3e-6,
and one cell at alpha 1 with 7e-21 kg. Their face forces (1e-8 N) come from the
stress and the impedances, which use the material density rho/alpha bounded to 0.5 to
2 rho0, but the velocity update divided by the true mass: accelerations of 1e9 m/s2,
speeds of 1e8 m/s in a step, then a singular least-squares matrix in the remap. The
sliver rule never fires there (the tube's wall faces are vertical); any change in the
last bits decided whether such a residue formed.

## Decision

1. The Lagrange step's inertia (predictor and corrector, before the group sums) is
   the cell's mass; but with free surface a residue, a cell holding less than 1e-6 of
   a full cell's mass (m < 1e-6 rho0 Vo), gets at least 0.5 rho0 alpha Vo, the bound
   the impedances already use. Every other cell is unchanged to the bit. Mass is not
   changed, so conservation stays exact; the residue's momentum is not exactly
   conserved, which is negligible at 1e-18 kg. A floor for every cell (not only
   residue) broke `test_wedge_void.sh`, where expanding elastic material legitimately
   reaches rho/alpha below 0.5 rho0; a floor of 1e-3 rho0 alpha Vo was too weak to
   stop the silo run-away.
2. Merge the sliver fix of decision 0092 (empty slivers above a full target stay
   empty, and the target keeps all the mass). `test_body_above_surface.sh` moves to CI.

## Results (10 mm cells)

| Check | Main | Now | Limit |
| --- | --- | --- | --- |
| Footing 2.5 mm above sand at rest, psi 0 / 30 | 0.17 m/s / blows up | 8e-17 m/s both | 1e-3 m/s |
| Projectile rigid / steel: penetration | 23.84 / 23.89 mm (-26 %) | 27.52 / 27.60 mm (-15 %) | 30 % |
| Projectile steel: density | 1.021 | 1.015 | 1.03 |
| Deformable push steel: density | 1.016 | 1.007 | 1.02 |
| Wedge silo quarter C at maxCo 0.399 / 0.4 / 0.401 | crashes with the sliver fix | pass all three | 0.45 to 0.75 |
| Steel push density at maxCo 0.399 / 0.4 / 0.401 | 1.014 to 1.032 (variant) | 1.008 / 1.007 / 1.007 | 1.02 |

## What remains

The steel push's corner group beside the footing's left edge (cells x 0.09 to 0.11 m)
still reaches 1.032 rho0 at t = 0.12 s in one variant of the mixing arithmetic
(ga/gV for the group's fill instead of r ga/gm, equal in exact arithmetic): the
neighbour on its left starts emptying and moving right while the group is driven left
to 3 m/s. The shipped form keeps r ga/gm (the form decision 0092 tested). This is a
free-surface corner effect, not residue; it is a known issue for a later package.
