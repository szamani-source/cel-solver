# 0039 Inclined layer: numerical viscosity measured, tilted gravity fixed (WP 3.6b)

Date: 2026-11-17. Status: accepted.

## Context

Decision 0036 estimated that the face Riemann solver adds a viscosity of order rho c h to
plastic flow, about 100 times the mu(I) viscosity at T0 (survey, "The problem in one
paragraph"). WP 3.6b measures it with a steady inclined layer (Bagnold profile) in
celFoam's cycle, with T0 Mohr-Coulomb and the placeholder Jop law of decision 0037.

Set-up (`tests/layer/inclined_layer.py`): layer h = 0.05 m on a 31 deg slope, x periodic,
gravity (g sin theta, -g cos theta) on an untilted mesh, sticking floor (wall friction 89
deg), void above (or a smooth lid), shear cap on, started from the exact steady state and
run 4 s. Exact solution of our model: Mohr-Coulomb in steady isochoric simple shear has
tau/sigma_n = sin(phi_m), so sin(phi_m) = tan(theta), one I = 0.0683 through the depth,
P = (2 + 2 nu)/3 |sigma_n|, surface speed 1.133 m/s.

## Decision

1. **Tilted gravity fix.** The hydrostatic face reconstruction (decision 0030) assumed a
   column along g: its jump -rho (g . d) gHat gHat is not periodic across a slope, and the
   layer slid uphill. With g not along an axis the jump is now that of a layer uniform
   across g's dominant axis e: -rho (d . e) (g e + e g - (g . e) e e). For g along an
   axis the old expression is kept, so every earlier result is unchanged bit for bit. A
   layer at 25 deg (below the stopping angle) now stays at rest below its top cell to
   2.4e-5 m/s (4.8e-4 with the old jump); `tests/test_inclined_layer.sh` checks it in CI.
2. **Measurement (done-when not met).** Speeds against the exact cell averages after 4 s,
   5 mm cells, free surface (lid in brackets where it differs):

   | E (MPa) | cell | body below the top cell (u / u_exact) | top cell | eta_num (Pa s) | eta_num / (rho c_s h) |
   | --- | --- | --- | --- | --- | --- |
   | 0.1 | 5 mm | 0.14 to 0.20 (0.12 to 0.18) | 0.14 | 2.8 (3.0) | 0.07 |
   | 1 | 5 mm | 0.05 to 0.09 | 0.85 (0.92) | 6.5 (8.0) | 0.05 to 0.06 |
   | 10 | 5 mm | 0.01 to 0.02 | 2.6 (2.7), still accelerating | 29 (22) | 0.05 to 0.07 |
   | 0.1 | 2.5 mm | 0.22 to 0.33 (0.20 to 0.29) | 0.22 | 1.7 | 0.08 |

   eta_num = sum(tau_exact - s_xy) / sum(gamma_dot) over the lower half (includes some
   inertia where the run is not yet steady). It scales with rho c_s h, as decision 0028
   reasoned: about 0.06 rho c_s h. The mu(I) viscosity at the base is 13 Pa s (secant)
   but only 0.8 Pa s as a tangent, and the tangent decides the speed: for 5 % on the
   surface speed eta_num must be below about 0.04 Pa s, 70 times lower than now at 0.1
   MPa and about 700 times at 10 MPa. **mu(I) alone cannot fix the E dependence** (the
   survey's estimate holds); WP 3.6c (plastic impedance at flowing faces) is needed.
3. **Top boundary cell slides.** The top cell of the layer, under void or a smooth lid,
   decouples and slides over the slowed body at 1 and 10 MPa (2.6 times the exact surface
   speed at 10 MPa, still accelerating); at rest at 25 deg it creeps at 2e-3 m/s while the
   body is at 2e-5. Likely cause (not checked): the shear cap holds its lower face at the
   static yield, below the slope's stress ratio (tan 31 = 0.60 > sin 35 = 0.57), and the
   top cell has little pressure. Recorded, not fixed: WP 3.6c revisits the cap.
4. The done-when stays as written; `tests/manual/test_inclined_layer.sh` is its test and
   fails now. WP 3.6c's done-when ("3.6b profile within 5 % at all three E") uses it.

## Consequences

- WP 3.6c should check its impedance against this table first (E = 0.1 MPa, 5 mm, 5 s
  run), then at 1 and 10 MPa.
- The shear cap of rule 8 uses the static phi even with mu(I); with mu(I) a flowing face
  needs phi_m. WP 3.6c should cap with the face's phi_m or drop the cap where the
  impedance replaces it.
- The fix in item 1 matters only for gravity off the mesh axes; silo, collapse and
  outlet runs are unchanged.
