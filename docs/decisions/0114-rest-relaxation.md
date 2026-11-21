# 0114: Top fill: relax locked-in lateral stress at rest

Date: 2026-11-21. Status: accepted (Sia approved the three-step plan,
2026-11-21 05:16 UTC: slower landing first, rest relaxation if K stays above
about 0.5, cap-coupled K0 only if both fall short).

## Context

Decision 0113 (K.3): a box filled from the top keeps K = s_xx/s_yy about
0.87 in the middle of the deposit against 0.42 when the same mass is placed
with gravity on, with or without the compaction cap. Grains landing at
speed push sideways and the frictional models keep that push once they stop
(the stress is inside the cone, so nothing yields it away). This is the
doubling of the Couto silo wall pressures in decision 0108.

## Step 1: slower landing (study, `tests/compaction/fill_box.py`, env S_IN)

Same mass, cap on (lambda* 0.02), box 0.2 x 0.6 m, 40 mm lid inlet:

| inlet speed | K mid | wall s_xx below mid-height |
|---|---|---|
| 2 m/s | 0.87 | 3.3 to 3.9 kPa |
| 0.5 m/s | 0.86 | 2.8 to 3.7 kPa |
| 0.25 m/s | 0.80 | 2.7 to 3.1 kPa |
| gravity on | 0.42 | 1.4 to 2.4 kPa |

The landing speed is set by the free fall from the lid (about 3 m/s), not
the inlet speed, so step 1 does not help. Filling in thin lifts would need
an inlet that rises with the surface; not pursued.

## Step 2: rest relaxation (built)

`materialProperties` entry

    restRelax { time 0.05; speed 0.01; K0 0.43; }   // K0 optional

In a full cell slower than `speed` [m/s], the stress block normal to
gravity relaxes towards K0 times the stress along gravity, with time
constant `time` [s]:

    s += min(dt/time, 1) (K0 s_v P_h - P_h s P_h),  P_h = I - n n, n = g/|g|

applied after the constitutive update, only where s_v < 0. Shear between
the vertical and the horizontal (wall friction) is untouched. K0 defaults to
Jaky's 1 - sin(phi). Off by default; needs gravity.

| run (cap on) | K mid | wall s_xx below mid-height |
|---|---|---|
| top fill, no relaxation | 0.87 | 3.3 to 3.9 kPa |
| top fill, time 0.05 s | 0.42 | 1.6 to 3.0 kPa |
| top fill, time 0.2 s | 0.42 | 1.6 to 3.0 kPa |
| top fill, 0.5 m/s inlet, time 0.2 s | 0.42 | 1.4 to 2.9 kPa |
| top fill, no cap, time 0.05 s | 0.42 | 1.6 to 3.0 kPa |
| gravity on, time 0.05 s | 0.43 | 1.4 to 2.4 kPa |

Test: `tests/test_rest_relax.sh` (about 3 minutes).

## Consequences

- It is a numerical aid, not a material law: it stands in for the slow
  creep and vibration that let a real fill settle to K0. Use it only for
  the at-rest phase; with `speed` small it is inactive in flowing cells.
- The gravity-on state is unchanged (K 0.42 to 0.43).
- Step 3 (cap-coupled K0) is not needed for the box.
- The wall stress near the base of the top-filled box stays about 20 %
  above gravity-on (denser, lower deposit).

## Couto silo top fill with relaxation (decision 0108 case)

`tests/hopper/top_fill.py` with `restRelax { time 0.05; speed 0.01; }`
added to the T0 rules plus fullBlend, MC, 4 cores, 3.5 h. Wall pressure
at rest, kPa:

| | top fill | top fill + relax (1 s settle) | gravity on | Couto |
|---|---|---|---|---|
| GA0 | 16.8 | 10.2 | 12.4 | 8.6-9.6 |
| GA1 | 9.4 | 7.4 | 4.6 | 4.0-4.5 |
| GA5 | 8.0 | 4.4 | 2.8 | 2.8-3.2 |
| GA7 | 4.3 | 2.6 | 1.3 | 1.4-1.7 |

Relaxation removes about half of the excess in the cylinder and brings GA0
closest to the experiment of any start, but the pressures were still
falling at the end (GA1 8.0, 7.6, 7.4 at 0.2, 0.8, 1.0 s after the inlet
closed): the silo cells only rest once damping has stopped the flow. A
2 s settle is running.
