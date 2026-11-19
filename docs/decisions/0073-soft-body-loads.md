# 0073 Material loads on deformable bodies: element faces, spread over the wall face (WP 5.3b)

Date: 2026-11-18. Status: accepted.

## Context

Decision 0072 coupled hex8 bodies to the material. A steel footing worked; a soft one
(E 20 MPa, the order of the material's own stiffness) packed the material under it to
1.62 rho0 and shortened 1.44 times the elastic column. Its contact force chattered from
0 to 50 N step to step, where the rigid footing's force is smooth.

What was tried (soft footing, 10 x 10 x 1 elements, 1 cm cells, maximum density):

| Change | Density |
|---|---|
| as in 0072 | 1.62 |
| implicit dashpot of the material's impedance on the loaded nodes | 1.51 |
| body damping 1, 2, 4 (subcycle shortened to suit) | 1.48 to 1.50 |
| hourglass 0.3 / 1 (triangle loads) | 1.32 / 1.17 |
| loads to the element face's 4 nodes, bilinear | 1.06 |
| ... and spread over 2 x 2 points of each cell's wall face | 1.043 |
| ... with hourglass 0.3 | 1.016 |
| ... 20 x 20 x 1 elements, hourglass 0.05 (point loads at the centre: 1.68) | 1.039 |

The cause was the element's hourglass modes. A wall cell's force went to the 3 nodes
of the nearest surface triangle. A face is two triangles, so one triangle's weights
load the face's hourglass mode, which the one-point element barely resists. The bottom
surface zigzagged, and the material, seeing that as wall velocity, packed. A body finer
than the cells is worse: nodes between the cell centres get no load at all.

## Decision

1. **Loads to the element face.** A cell's force goes to the 4 nodes of the nearest
   element face with the face's bilinear shape functions, at the nearest point (Gauss-
   Newton from the triangle's weights). This is the load the element itself sees.
2. **Over the wall face.** Each cut cell's force is split over 4 points, the 2 x 2
   Gauss points across the cell (in the plane normal to its wall face, the cell's
   extent there), and the wall velocity of a deformable body's cut cell is the mean of
   the surface velocity at the same points. The material loads the surface over the
   cell's width and sees its mean motion there, which is also the volume the cell's
   wall face sweeps. Rigid bodies are unchanged (velocity at the cell centre).
3. **Soft bodies use `hourglass 0.3`** (stiffness form). The test does; the default
   stays 0.05 (decision 0070), which serves stiff bodies.
4. The hex8 stable step accounts for `damping` (critical fraction z of the shortest
   mode: dt times sqrt(1 + z^2) - z). Without it damping 1 diverged.
5. Not kept: the impedance dashpot. It halved the error with triangle loads but made
   no difference (1.0155 against 1.016) once the loads were consistent.

## Consequences

- `tests/test_soft_push.sh` (about 1.5 min) passes: density 1.0155 rho0, shortening
  0.90 of the elastic column (band 20 %; side friction carries part of the load).
- The steel footing of WP 5.3a also improves: mean force -0.5 % of the rigid one (was
  -3 %), density 1.008 rho0 (was 1.013).
- A body much finer than the cells still loads its nodes unevenly (1.039 at 20 x 20
  with hourglass 0.05). Keep body elements near the cell size, or raise `hourglass`.
