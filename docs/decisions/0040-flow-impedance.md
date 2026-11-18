# 0040 Flow impedance: less numerical resistance in flowing material (WP 3.6c)

Date: 2026-11-17. Status: accepted (switch, default off).

## Context

WP 3.6b (decision 0039) found the mu(I) inclined layer flowing at 0.14 to 0.20 of the
exact speed at E = 0.1 MPa and 0.01 to 0.02 at 10 MPa (5 mm cells), and read it as a
numerical viscosity of about 0.06 rho c_s h from the face Riemann solver. Decision 0036
planned an impedance from the plastic tangent at flowing faces.

Experiments in WP 3.6c (inclined layer, 0.1 MPa, 5 mm, scaling terms one at a time)
found three sources, not one:

1. **The Hancock half-step predictor** adds dt/2 C : sym(grad U) to the face stress
   states with the elastic C. In a flowing cell the true stress rate is almost zero, so
   the faces see a shear stress G dt gamma_dot / 2 too high. G dt scales with
   rho c_s h, which is the 3.6b scaling. With the shear cap on, this was the largest term:
   dropping the deviatoric predictor alone took the layer from 0.15 to about 1.15 of
   exact.
2. **The shear jump term** z_s (v_R - v_L)/2 of the Riemann solver. Without the shear
   cap it dominates (0.18 of exact with the predictor fixed). Removing it entirely lets
   odd-even velocity modes grow (checkerboard): some jump term is needed.
3. **The shear cap at the static friction angle** (rule 8) caps faces of flowing cells at
   phi, below phi_m, and the capped faces slip, so the cells see little strain rate and
   stay near static yield. With 1 and 2 fixed and the cap kept as it was, the layer ran
   too fast and kept accelerating.

The mu(I) law is flat: its tangent viscosity at the 3.6b layer is 0.5 to 0.8 Pa s, so a
face stress error of 1 Pa (in 200 to 400 Pa) changes the shear rate by about 5 %. The
steady layer is a sensitive test.

## Decision

`flowImpedance yes` in materialProperties (Mohr-Coulomb with mu(I) only; default off,
every earlier result unchanged bit for bit):

1. **Flowing cell**: a full cell whose stress is within 3 % (in sin phi) of the static
   Mohr-Coulomb surface or beyond it. Its mobilised angle sin(phi_m) = (s1 - s3)/(2 p_m)
   gives mu = tan(phi_m), I = I0 (mu - mu_s)/(mu_2 - mu) and the tangent viscosity
   eta = d tau/d gamma_dot = p_m (1 + mu^2)^-1.5 (mu_2 - mu_s) I0/(I0 + I)^2 d sqrt(rho_s/P),
   finite at I = 0 (about 1 Pa s at the layer's base). The 3 % band stops cells that hover
   at the static surface from switching every step; with a band of 1e-6 the 10 MPa layer
   locked near its base (body 0.4 to 0.65 of exact against 0.45 to 0.74).
2. **Face shear jump**: each side contributes z_s, or 2 eta/delta if flowing (the smaller),
   in series; delta is the distance between cell centres. Walls: z_s or eta/(distance to
   the wall). The face velocity keeps the elastic impedances. The jump now acts as the
   physical viscosity on the compact stencil (it damps odd-even modes with eta, not
   rho c_s h).
3. **Hancock predictor** in a flowing cell: the deviatoric elastic rate along the plastic
   flow direction m = (e1 e1 - e3 e3)/sqrt(2) (principal directions of s1, s3) is dropped
   while loading (m : D > 0). Dropping the whole deviatoric rate was better at 1 MPa but
   unstable at 0.1 MPa, 5 mm.
4. **Shear cap**: skipped between two flowing cells (the impedance replaces it); elsewhere
   it caps at max(phi, phi_m of the flowing side), not at phi.

## Results (inclined layer, 4 s from the exact state; u / u_exact of the cell averages)

| E (MPa) | cell | WP 3.6b body | WP 3.6c body (below the top cell) | WP 3.6c top cell | profile error |
| --- | --- | --- | --- | --- | --- |
| 0.1 | 5 mm | 0.14 to 0.20 | 0.65 to 0.97 | 0.93 | 7.8 % |
| 1 | 5 mm | 0.05 to 0.09 | 0.41 to 0.92 | 1.22 | 22 % |
| 10 | 5 mm | 0.01 to 0.02 | 0.03 to 0.74 | 2.37, accelerating | 136 % |
| 0.1 | 2.5 mm | 0.22 to 0.33 | 0.53 to 1.03 | 1.03 | 4.3 % (passes 5 %) |
| 1 | 2.5 mm | | 0.33 to 1.01 | 0.99 | 4.4 % (passes 5 %) |
| 10 | 2.5 mm | | 0.01 to 0.88 | 2.94, accelerating | 194 % |

With a smooth lid instead of void the top cell slides too (0.1 MPa: 1.23; 10 MPa: 3.16)
and the body is a little slower (0.61 to 0.88 at 0.1 MPa).

The largest profile errors at 0.1 and 1 MPa are now the two cells next to the floor,
which converge with the mesh. The done-when (5 % at 5 mm for all three E) is **not met**:

- At 10 MPa the bottom cell stays nearly stuck and the body is 25 to 35 % slow, still
  slowing at 4 s. Not found in this package; candidates are the wall face and the cells
  that hover at the static surface near the floor.
- The top cell, next to void, still slides over the body at 1 and 10 MPa (decision 0039
  item 3). At 10 MPa it carries no stress at all: a tiny expansion unloads it (void-first
  closure, decision 0021) and it then has no friction. Counting it as flowing (alpha >=
  0.5) made it worse (it flew off and the body locked).

## Consequences

- Silo, collapse and outlet runs do not use mu(I) and keep the switch off: unchanged.
  WP 3.6d turns it on for the silo with mu(I).
- WP 3.6c2 (split): the 10 MPa layer (stuck floor cell, slow body) and the sliding top
  cell; done-when unchanged (5 % at 5 mm, all three E).
- CI guard in `tests/test_inclined_layer.sh`: at 0.1 MPa, 5 mm the surface is within 10 %
  and the mean speed at least 0.85 of exact (3.6b code: 0.14 and 0.15).
- The infinite slope at 31 deg is statically stable for Mohr-Coulomb with phi = 35 deg
  when s_xx differs from s_yy; only the coaxial simple-shear state flows. Some failed
  variants locked the whole body at rest this way. Any disturbance of the layer (a
  sliding top cell) can tip it into the static state.
- The 3 % band treats cells slightly below the static surface as flowing with a low face
  impedance (eta finite). Not seen to harm static regions in the at-rest layer (largest
  speed 6.4e-5 m/s, was 2.4e-5, limit 1e-4); to be watched in the silo (WP 3.6d).
