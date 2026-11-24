# 0124 Compaction cap with STL (immersed) walls, and the plan K.4-K.7

Date: 2026-11-23. Status: accepted (Sia's pre-gas list item 2, "proceed",
2026-11-23 22:50 UTC).

## Context

The compaction cap (`capLambda`, `capP0`, decision 0113) refused to run with
immersed bodies, and every validated Couto run uses an STL shell inside the
wall-aligned mesh (decision 0107). So the cap was never tried on the hopper,
where GA0 at rest reads 36 % high (decision 0117).

## Decision

Plan (one PR each where code changes):

| WP | Title | Done when |
| --- | --- | --- |
| K.4 | Cap with STL walls | The K.2 settling column between frictionless STL walls (a quarter-open small cell mixed into its neighbour) lies on the line within 2 % and settles as the mesh-wall column within 1 % (`tests/test_cap_stl.sh`) |
| K.5 | Wheat compaction values; Couto at rest with the cap (coarse mesh) | lambda* for wheat from a published compression test; table GA0-GA7 cap on and off against Couto |
| K.6 | Pressure-dependent stiffness (only if K.5 shows the cap alone barely moves GA0) | Kernel and settling-column checks |
| K.7 | 4 mm Couto with the cap, MC and DP | At rest GA0-GA5 against Couto; discharge rate inside Beverloo 110-131 kg/s |

K.4 as built: `pc` is carried through the three immersed-wall steps like
`epsPEq` (ghosts of covered cells, small-cell mixing, hand-over when a body
moves), with one difference: pc = 0 means "not yet loaded", so only loaded
material enters the averages. Averaging that sentinel in mass-weighted
gave a cell pc = 25 Pa (a quarter of P0) in the first test run; the cap
then compacted it below P0 and left a 0.5 % density bump (11 % off the
line). `celVoidCheck` writes pc as the last cells.dat column when the cap is
on.

## Consequences

- `tests/test_cap_stl.sh` (about 30 s): 47 cells on the line within 0.53 %,
  p_c = P, mass 6e-13, settlement 36.34 mm against 36.31 mm (0.07 %).
- Runs without the cap are unchanged (the pc arrays are not touched).
- Known: the ordinary remap still mass-weights pc with the sentinel, so
  inlet-fed material (pc 0) mixing into loaded material lowers the mix's pc.
  Harmless for the gravity-on hopper; revisit if a fed case uses the cap.

## K.5 (2026-11-24): the cap does not move GA0 with the gravity-on start

Wheat value: lambda* = Cc/(ln 10 (1 + e0)) = 0.08/(2.303 x 1.53) = 0.023,
with the mean compression index Cc 0.08 of the confined compression tests in
Zeng and Wang (2019), Int. Agrophys. 33, 347-354, doi 10.31545/intagr/110809
(wheat, 11.5 % moisture, 12.5-600 kPa), and e0 = 1310/856 - 1 = 0.53.

`tests/hopper/run_cap_study.sh mc cap|nocap coarse` (12 mm rows, nonlocal,
gravity-on start, 1 process each), wall pressure at 0.5 s, kPa:

| | GA0 | GA1 | GA2 | GA3 | GA4 | GA5 | GA6 | GA7 |
|---|---|---|---|---|---|---|---|---|
| Couto at rest | 8.6-9.6 | 4.0-4.5 | 4.1-4.5 | 3.7-4.2 | 3.2-3.6 | 2.8-3.2 | 2.5-2.7 | 1.4-1.7 |
| no cap | 12.41 | 4.56 | 4.58 | 4.03 | 3.44 | 2.79 | 2.08 | 1.30 |
| cap | 3.30 | 0.78 | 0.30 | 0.08 | 0.06 | 0.05 | 0.03 | 0.02 |

The cap run was not at rest: the wheat starts at Couto's 856 kg/m^3, which
is already the filled (compacted) density, and the cap compacts it again
from P0 = 100 Pa (9 % denser at the base), so the upper column was still
falling at the damping-limited 0.1 m/s at 0.5 s. That counts compaction
twice. Started consistently (normally consolidated under its own weight)
the cap is inactive at rest and GA0 is the no-cap 12.4 kPa: with the
gravity-on start the cap cannot change the at-rest pressures. It acts only
where the loading path matters: filling. Both runs were stopped at 0.5 s.

Closed here (Sia, 2026-11-24 02:28 UTC: do not pursue the GA0 gap; Sia
expects a finer mesh at that point to resolve it). K.6 and K.7 are not
needed. Ready for a later filling study, not run: top fill with the cap, wheat
entering loose at 772 kg/m^3 (standard test weight, Bhadra et al. 2018,
Trans. ASABE 61(2), doi 10.13031/trans.12645) until the Couto mass 1503 kg
is in, restRelax, 2 s settle: `tests/hopper/run_cap_topfill.sh`
(checkpointed; `top_fill.py` gained CEL_RHO, CEL_CHECKPOINT, CEL_RESUME).
Reference without the cap: decision 0115 run 2.

## Parallel fix found on the way

Top fill on 2 processes stopped at the first step with MPI_ERR_TRUNCATE, on
main before K.4 as well. In the free-surface remap the "turned cells" flag
(WP G.3) was local: a processor without turned cells skipped the turn back
and its correctBoundaryConditions, so the processor exchanges no longer
matched. The flag is now reduced over all processors (returnReduceOr).
