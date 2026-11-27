# 0138 Mach-scaled normal jump between yielding cells (WP F.8k)

Date: 2026-11-27. Status: opt-in, under test (Sia chose "add the slow-flow fix").

## Context

After decision 0137 the plain Mohr-Coulomb dry dam break still ran 1.85 / 2.57 / 2.81 h0
on 10 / 5 / 2.5 mm cells (front at 5 mm deposit thickness), and 10 mm changed with the
time step (2.35 h0 at Courant 0.1). Tests in WP F.8j: a block translates exactly and
stops on a rough floor at 0.990 m (exact 0.992); first-order face states halve the
run-out; the face stress caps change nothing. So the loss is the face flux damping,
which works at the elastic wave speed (16-30 m/s) while the grains move at 1-2 m/s.
Godunov schemes are known to over-damp such low-Mach flow; the published cure scales
the jump term by the local Mach number (Guillard & Viozat 1999, Thornber et al. 2008,
Rieper 2011).

## Decision

`materialProperties slowFlowJump yes` (needs `yieldShearJump no`): at a face between
two yielding cells the normal jump coefficient zpL zpR/(zpL + zpR) of the acoustic
Riemann solver is multiplied by max(|vL|, |vR|)/c_p, kept within [slowFlowMin, 1]
(slowFlowMin default 0.01), on internal and coupled faces. Elastic material and fast
impacts keep the full jump. Default off.

## Consequences

Dry dam break, plain Mohr-Coulomb (tests/manual/test_collapse_convergence.sh), run-out
at 1.4 s in h0, by the deposit thickness that defines the front:

| Front at thickness | 10 mm | 5 mm | 2.5 mm |
| --- | --- | --- | --- |
| 5 mm, before / after | 1.85 / 2.30 | 2.57 / 2.46 | 2.81 / 2.84 |
| 10 mm, before / after | 1.82 / 2.05 | 2.14 / 2.16 | 2.04 / 2.29 |
| 20 mm, before / after | 1.40 / 1.58 | 1.67 / 1.69 | 1.75 / 1.62 |
| 30 mm, before / after | 1.05 / 1.27 | 1.34 / 1.38 | 1.31 / 1.38 |

- The body of the deposit (front at 20 mm) now agrees within 4 % of 2.5 mm on all
  three meshes (before: 10 mm 20 % short). 10 mm at Courant 0.4 gives 2.30, at 0.1
  2.42 (before 1.85 / 2.35): the time-step dependence is mostly gone.
- At 5 mm thickness the run-outs are 2.30 / 2.46 / 2.84: each mesh ends the deposit's
  thin tail where it is one cell thick (WP F.8k log), so that measure follows the
  mesh. Sia (2026-11-27) chose to measure the front at a thickness every mesh resolves:
  the done-when uses 20 mm (two coarse cells) and passes, 3.8 %; the dam break is
  calibrated on 5 mm cells.
- Default unchanged, so existing tests are unchanged (the kn < 0 path is the old
  formula).
