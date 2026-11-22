# 0115: Couto silo with 12 mm cylinder rows, and top fill at the Couto mass

Date: 2026-11-22. Thread "Close Couto silo gaps on finer mesh".

## Context

Two Couto silo gaps were left open (decisions 0109, 0114):

- In discharge GA1 is 25-30 % low and GA2 11-14 % low. The nonlocal law's
  cooperative length (about 3 d = 12 mm) was below the 25 mm rows of the
  upper cylinder, so it could not widen the converging zone.
- After a top fill with `restRelax` the upper cylinder reads 40-80 % high.
  That run put 1541 kg in for 1503.

## Changes

- `tests/hopper/gap_case.py`: `CEL_H_TOP` (m, default 0.025) sets the row
  height in the upper cylinder; the headroom stays 0.3 m. 0.012 gives
  12,336 cells (was 7,584); the time step (5.85e-6 s) is set by the outlet
  cells and does not change. The start picture has an upper-cylinder panel.
- `tests/hopper/top_fill.py`: the overfill was the inlet's size on the mesh.
  topoSet takes the top faces whose centre is inside R_IN, so the inlet ends
  at the outer edge of a column, 0.2531 m instead of 0.25 m: 2.5 % more area,
  exactly 1541/1503. The script now meshes once, reads the inlet's real
  radius and times the inlet from it. `CEL_FILL_MASS` (kg, default the Couto
  mass) sets the mass fed.
- numpy < 2 fallback for `np.trapezoid`.

## Results (MC, T0 rules plus fullBlend, 4 cores)

Run 1: 2-cell wall-gap start on the rounded STL mesh, nonlocal A 0.48,
spreadLayers yes, upper rows 12 mm (3.7 h). Discharge mean 0.25-1.0 s after
opening, kPa:

| | rate kg/s | GA0 | GA1 | GA2 | GA3 | GA4 | GA5 | GA6 | GA7 |
|---|---|---|---|---|---|---|---|---|---|
| 25 mm rows (0109) | 119.4 | 24.1 | 5.63 | 5.03 | 4.51 | | | | |
| 12 mm rows | 119.7 | 23.9 | 5.25 | 5.39 | 4.61 | 4.15 | 2.88 | 2.05 | 1.18 |
| Couto | (Beverloo 110-131) | 18.7-24.8 | 7.1-8.2 | 5.8-6.5 | 4.4-5.6 | 3.2-4.0 | 2.9-3.7 | 3.4-3.6 | 1.4-1.7 |

No leak; mass at opening 1503 kg. GA2 rises 0.4 kPa (now 7 % below the
range), GA1 does not rise (26 % low), the rate and GA0 do not change.
GA6/GA7 are low in discharge. At rest (t 1.0 s) this start reads GA0 17.3,
GA1 5.1, GA2 6.9, GA3 6.9: well above the experiment, as the wall-gap start
was before (decision 0108); the gravity-on start stays the default for the
at-rest pressures.

Run 2: top fill through the r 0.25 m inlet, restRelax, 2 s settle, now
1503 kg in (4.9 h). At rest, kPa:

| | GA0 | GA1 | GA2 | GA3 | GA4 | GA5 | GA6 | GA7 |
|---|---|---|---|---|---|---|---|---|
| 1541 kg (0114) | 9.3 | 6.8 | 6.3 | | | 4.2 | | 2.6 |
| 1503 kg | 9.27 | 6.80 | 6.16 | 5.37 | 4.63 | 3.99 | 3.23 | 2.33 |
| Couto | 8.6-9.6 | 4.0-4.5 | 4.1-4.5 | 3.7-4.2 | 3.2-3.6 | 2.8-3.2 | 2.5-2.7 | 1.4-1.7 |

No leak, largest speed at the end 0.004 m/s, heap top 2.664 m (silo top
2.476 m). The right mass lowers the upper cylinder by 5-10 %, but it stays
25-60 % high, and the hopper and lower cylinder do not move.

## Decision

- Keep both options. The inlet fix stays on (it is a bug fix: the asked mass
  now goes in).
- The 12 mm rows do not close GA1. With 12 mm rows the cells are still
  about the cooperative length, not several cells per length; resolving
  it needs cells near d (4 mm) in the cylinder, about 9 times the cells and
  a smaller time step there is not needed (outlet sets it), so about 9x the
  run time (about 30 h on 4 cores). Not worth it before the Linux install;
  GA1 stays a known gap (Abaqus in Zheng and Yu is also low there).
- The top-fill excess in the upper cylinder is not the fill mass. What is
  left is the heap (0.19 m above the silo top, where Couto's fill is level)
  and the denser, impact-placed packing near the walls. The gravity-on start
  stays the default for validation.

## Consequences

STATUS rows G.5 and K.3 updated. Run directories were scratch; the logs are
copied above.
