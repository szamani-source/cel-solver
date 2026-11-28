# Solid and fluid tutorials

Version 2 of celFoam solves the gas in the pores of the material and in the void with
the granular solid (decisions 0125-0135). `constant/gasProperties` switches it on;
without that file a case runs exactly as version 1. The gas model holds for packed
beds that the gas does not lift or expand (decision 0135); the lift check in
`log.celFoam` warns when a run leaves that range.

| Tutorial | What it shows | Run time (one core) |
| --- | --- | --- |
| [`airImpededHopper`](airImpededHopper/README.md) | Fine sand discharging from a slot hopper, slowed by the air it must draw in (WP F.5); with gas and dry | about 2 minutes |
| [`fluidisedDamBreak`](fluidisedDamBreak/README.md) | A fluidised column of glass beads collapsing along a channel against Roche et al. (2008) (WP F.8); fluidised and dry | about 20 minutes (10 mm) |

Each folder has `Allrun`, `Allclean`, a `README.md` with the sections Setup, Numerics,
How to run and Expected results, and `plot.py` for its figures.
`tests/test_gas_tutorials.sh` checks the layout and runs a short stretch of each.
