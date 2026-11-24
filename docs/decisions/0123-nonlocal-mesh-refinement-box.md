# 0123: Mesh refinement of the nonlocal law on a box discharge; collapse toe tries

Date: 2026-11-24. Thread "Add nonlocal law for mesh dependence" (Sia's
pre-gas list item 1, 2026-11-23).

## Context

The nonlocal fluidity law (decision 0111, `nonlocalA`) was already in the
solver and used in the final Couto validation (decision 0117). What was
left of the mesh-dependence item: (a) the column-collapse toe still grows
with refinement (decision 0111, N.3: a thin film of partly filled cells,
which carry no stress, decision 0024); (b) no refinement check of a
discharge with and without the law. A Couto hopper refinement (8 / 5.7 /
4 mm, 1 to 8 h per run) was prepared (`tests/hopper/run_refine.sh`), but
Sia (2026-11-23 23:23 UTC) asked for a small 2D box instead, with a local
and a nonlocal run on every mesh.

## N.3b: two tries at friction for partly filled cells (no effect, removed)

`CEL_LAWS=NL CEL_NL_D=0.01 CEL_NO_SPREAD=1 tests/collapse/shear_band.py`,
run-out at 1 s (front where the layer is 5 mm), 5 / 2.5 mm:

| Variant | Run-out 5 / 2.5 mm |
| --- | --- |
| nonlocal, as main (decision 0111) | 0.225 / 0.250 |
| wall friction of a partly filled cell at least mu times its weight | 0.225 / 0.250 |
| plus shear-cap normal stress at least the film's weight on faces | 0.225 / 0.312 (stray mass reached the far wall) |

Neither helped, so the code was not kept. What the runs show instead: at
1 s the deposit is not at rest. Its surface layer (the partly filled cells
and one or two full rows under them, P 10 to 30 Pa, 349 full cells above
mu_s) creeps down the slope at 1 to 2 cm/s, with the kinetic energy at
about 1e-3 of its peak, and the creep is slower on the finer mesh. The late
run-out is that creep, not a lack of floor friction. The cure belongs with
the free-surface stress of partly filled cells (decision 0024's open item)
and the rate effect (decision 0096, parked), not with the law.

## N.4: box discharge, local against nonlocal on three meshes

`tests/nonlocal/box_refine.py` (checkpoints every 0.05 s, resumes; the
runs survived two container restarts): the T0 tank of WP 3.3 (W 0.4 m, fill
0.6 m, outlet 60 mm in the middle of the floor, walls 20 deg, MC phi 35,
E 10 MPa as `tutorials/siloDischarge`, opens at 0.3 s), mu(I) set of
decision 0037 with grain size d = 5 mm (outlet 12 d), flowImpedance and
normalCap on. Local = mu(I); nonlocal = the same plus `nonlocalA 0.48`.
Cells 10 / 5 / 2.5 mm (2 / 1 / 0.5 d).

Rate: least squares over the tank mass from 97 to 88 % of the mass at
opening (WP 3.3). Speed: vertical speed on the centreline at y = 0.15 m
(2.5 D above the outlet), 0.3 and 0.5 s after opening. Change: from the
next coarser mesh.

| Law | Cells | Rate kg/s/m (C) | Change | Speed 0.3 / 0.5 s, m/s | Change |
| --- | --- | --- | --- | --- | --- |
| local | 10 mm | 62.5 (0.754) | | 0.123 / 0.136 | |
| local | 5 mm | 65.2 (0.787) | +4.3 % | 0.139 / 0.154 | +13 / +13 % |
| local | 2.5 mm | 62.4 (0.753) | -4.3 % | 0.155 / 0.159 | +12 / +3 % |
| nonlocal | 10 mm | 62.4 (0.753) | | 0.162 / 0.210 | |
| nonlocal | 5 mm | 70.1 (0.847) | +12 % | 0.237 / 0.243 | +46 / +16 % |
| nonlocal | 2.5 mm | 66.8 (0.806) | -4.7 % | 0.232 / 0.243 | -2 / 0 % |

Mass balance 4e-14 on runs without a restart, 2e-7 on runs resumed after
a container restart (the discharged-mass counter, as decision 0117).

Findings:

- The flow field over the outlet converges with the nonlocal law (5 to
  2.5 mm: 2 % and 0 %) and does not with the local law (12 % and 3 %, still
  rising). The nonlocal law also concentrates the flow: the centreline is
  about 55 % faster than with the local law on the same mesh.
- The rate is not strongly mesh dependent here with either law: 5 to
  2.5 mm moves both by 4 to 5 %, within the run-to-run scatter of this tank
  (decision 0046, about 15 % at 10 mm). The local rates go up and down
  (62.5, 65.2, 62.4). The nonlocal rate is 5 to 7 % above the local one on
  the meshes that resolve d.
- On 10 mm cells (2 d) the two laws give the same rate: the law's length
  is not resolved and it acts as the local law, as N.2 found (cells of
  about 1 d near the band).

## Consequences

- Item 1 of Sia's pre-gas list (mesh dependence) closes with the law as it
  is: use `nonlocalA` with cells of about one grain diameter or finer where
  the material shears, and check a second mesh. Manual FAQ and STATUS say so.
- The collapse toe (run-out growing with refinement) stays a known issue
  under the free-surface stress of partly filled cells and the parked rate
  effect (decisions 0024, 0096), not under the friction law.
- The Couto hopper refinement script (`tests/hopper/run_refine.sh`) stays in
  the repository for a later check on the hopper; not run.
