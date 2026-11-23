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
