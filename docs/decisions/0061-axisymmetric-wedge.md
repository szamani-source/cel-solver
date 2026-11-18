# 0061 Axisymmetric runs on an OpenFOAM wedge mesh (WP 3.8a)

Date: 2026-11-18. Status: accepted (WP 3.8a done; free surface on the wedge is WP 3.8b).

## Context

Sia (2026-11-18): Abaqus CEL has no axisymmetric option, so a cylindrical silo needs a
quarter model. celFoam should run a wedge slice instead. OpenFOAM already has the
mesh: one cell thick in the hoop direction, two `wedge` patches a few degrees apart,
and the cells beside the axis collapse to prisms (blockMesh). Before this package
`kindOf` stopped with "wedge is not supported".

Inside the wedge every face is along r or y, so the Cartesian rules hold unchanged;
only the two wedge sides are tilted.

## Decision

1. **Wedge side faces** (`boundaryKind::wedge`, Lagrange step): the state beyond the
   face is the cell's own predicted state (the Hancock half step, in the predictor and
   the corrector alike), turned onto the face plane with OpenFOAM's `faceT`, as the
   `wedge` patch field does. With the old cell state instead, the wedge faces lagged
   the cell's other faces by the half-step velocity, so a plug falling freely sheared
   at the axis (shear stress growing as 1/r) and blew up. The face carries the rotated cell velocity and traction
   `(T sigma T^T) n`, with no Riemann jump. The two sides together give the hoop
   strain rate u_r/r and the hoop force -sigma_tt V/r. No volume flux crosses a
   wedge side (phiL = 0). An earlier try removed the normal velocity there; that let
   a swirl grow in the axis cells, so the full rotated velocity is kept.
2. **Wall stress extrapolated on wedge meshes** (`extrapolatedWalls`, second order
   only): on a wedge the hoop stress of a cell beside a wall balances the radial face
   forces over one cell, so the first-order wall traction (the cell value) left an
   O(1) hoop error at the wall (Lame: 40 % in the wall cell, not falling with h). On
   non-coupled faces normal to gravity (side walls) the face stress is the cell value
   plus the Gauss gradient made consistent with those face values,
   g_a = g0_a/(1 - K_aa), K = sum S d/V over the wall faces (skipped where
   1 - K_aa <= 0.25, one cell across). Faces along gravity keep the hydrostatic
   mirror (decision 0030). Cartesian meshes are untouched.
3. **Time step** (`cellSize()`): cbrt(V) goes to zero at the axis (V ~ r). On a
   wedge the length is the smallest in-plane V/A_a, combined with the cell's radius
   r = V/|S1 + S2| (its two wedge faces) as 1/d^2 = 1/h^2 + 1/r^2, because the hoop
   strain adds a spring of length r. With h alone the axis prisms (r = 2h/3) grew an
   exponential radial mode at maxCo 0.4 (stable below about 0.27). The combined
   length costs about 1.8 times the steps of h alone.

## Consequences

- `tests/test_wedge.sh` (about 1.5 min): Lame thick cylinder within 2.1 % (5 mm) and
  0.9 % (2.5 mm) of max |s_rr| beyond two wall cells, wall cells first order (error
  ratio 0.40 when h halves); a K0 lake with axis cells at rest to 5e-16 m/s;
  cylindrical Janssen tank settles, K = 0.43 at the wall and mean s_v within 2.0 % of
  circular Janssen; 2 processors match serial to 3e-17.
- The wall cells of a wedge stay first order (hoop 9 % at 2.5 mm in Lame).
- The axis of the wedge needs no patch (it has no faces); an `axis` patch type with
  faces is not supported.
- Not yet: `freeSurface` on a wedge (slabs are Cartesian, WP 3.8b), outlets and inlets
  on a wedge (3.8c), STL walls on a wedge (3.8d).
