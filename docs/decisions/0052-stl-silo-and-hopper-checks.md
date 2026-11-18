# 0052 How the STL silo and hopper are checked (WP 4.5)

Date: 2026-11-18. Status: accepted.

## Context

WP 4.5 asks for the T0 tank discharge with STL walls within 5 % of the mesh-wall tank,
and a 60 deg hopper run to empty with its rate recorded. A single 10 mm silo rate
scatters by about 15 % from round-off (decision 0046), so one run against one run cannot
resolve 5 %. A 10 MPa run takes about 9 minutes; 0.1 MPa about 25 s, and the rate does
not depend on E beyond that scatter (decision 0045).

## Decision

1. **Bands of three runs** per tank, maxCo 0.40, 0.38 and 0.36 (as decision 0046), at
   E = 0.1 MPa, everything else T0 with mu(I) and the flow impedance. The mean STL rate
   must be within 5 % of the mean mesh-wall rate.
2. **The STL tank is all STL:** side walls at x = 0 and 0.4 cutting cells at a quarter
   and three quarters, floor at y = 0.0125 cutting the second row, floor gap 0.06 m over
   the mesh outlet patch (the patch is 0.07 m; the floor closes its end faces in part).
   This is `tutorials/stlSilo`. The fill is 0.6 m above each floor.
3. **Hopper:** two STL walls at 60 deg to the horizontal (half angle 30 deg) from the
   gap edges up to the side walls; same fill top. Pass only on mass balance and
   emptying; its rate is printed beside the flat-tank mean and the Rose-Tanaka hopper
   factor (tan 30 deg)^-0.35 = 1.21, with no pass band until Sia agrees one.

## Consequences

- `tests/test_stl_silo.sh` (about 6 minutes, in CI): mesh walls 65.5 / 66.3 / 65.7,
  STL 68.2 / 68.0 / 70.1 kg/s/m, mean +4.4 % (all three STL runs above all three mesh
  runs, so the offset is real, not scatter; within 5 %). Hopper 66.2 kg/s/m, ratio 1.01
  to the flat tank (Rose-Tanaka 1.21), empty 5.0 s after opening.
- The +4.4 % sits close to the limit; a cause (the cut floor cells beside the gap, the
  gap region below the floor line, corner cells) is not isolated.
- The hopper does not show the faster flow the correlations predict for a 30 deg half
  angle. Not judged here: the correlations are for 3D orifices and dry grains with
  measured wall friction; this is a 2D slot at 10 mm cells with placeholder mu(I).
