# 0071 Viscous hourglass control for large plastic flow (WP 5.2)

Date: 2026-11-18. Status: accepted.

## Context

WP 5.1 chose the stiffness hourglass control (decision 0070): conservative, energy kept
to 1e-7 for an elastic body. In the Taylor bar (WP 5.2) the stiffness form is wrong: its
forces grow with the integrated hourglass motion at the elastic modulus while the
material itself flows at yield, so the hourglass forces end up far above the stresses
and hold the impact face back. 6 x 6 x 45 quarter model, reference (Abaqus C3D8R)
shortening 13.10 mm, widening 5.528 mm:

| hourglass control | shortening (mm) | widening (mm) |
|---|---|---|
| stiffness 0.2 | 15.26 | 2.39 |
| stiffness 0.05 | 13.97 | 3.74 |
| stiffness 0.05, 10 x 10 x 60 | 13.60 | 4.35 |
| stiffness 0.005 | 13.14 | 5.09 |
| viscous 0.1 | 13.12 | 5.67 |
| viscous 0.05 | 13.14 | 5.76 |

Volume is kept to 0.4 % in every run: the stiff runs bulge behind the face instead.

## Decision

1. `hourglassType viscous` (LS-DYNA type 1): hourglass forces from the current
   hourglass velocities, coefficient `hourglass` rho c V^(2/3) / 4; nothing is stored.
2. The default stays `stiffness` (elastic bodies keep their energy); bodies in large
   plastic flow (impact, penetration) use `viscous`. The Taylor bar test runs viscous
   0.05.

## Consequences

- `tests/test_taylor_bar.sh` (about 1.5 min): +0.3 % and +4.2 % against the Abaqus 3D
  model.
- The viscous form damps hourglass-shaped physical motion a little; not measured.
- A stiffness form scaled by the material's current tangent (Belytschko-Bindeman) would
  serve both cases; not done.
