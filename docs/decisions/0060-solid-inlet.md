# 0060 Solid inlet: a face that moves in, a swept volume filled with fresh material (WP 3.7a)

Date: 2026-11-18. Status: accepted (WP 3.7a done).

## Context

Sia (2026-11-18) named two things Abaqus CEL cannot do that celFoam should: an inlet
that feeds solid material into the domain, and axisymmetric runs on a wedge (WP 3.8).
Before this package a patch with a fixed velocity let material in only as a copy of the
cell beside it (first-order zero gradient, remap header), and with `freeSurface yes`
(every silo run) nothing at all came in: the remap's slab rule found no material
outside the mesh.

## Decision

1. **Dictionary** `inlets` in `constant/materialProperties` (needs `freeSurface yes`),
   built like `outlets` (decision 0033):

   ```
   inlets
   {
       feed
       {
           speed      0.5;      // Function1 of time [m/s], into the domain
           startTime  0.05;     // default 0
           endTime    0.35;     // default never
           density    1800;     // default rho
           pressure   0;        // inlet stress -p I [Pa], default 0
           friction   0;        // closed-wall friction [deg], default 0
       }
   }
   ```

   The patch must be a plain patch without a fixed velocity. Outside
   [startTime, endTime) the inlet is a closed wall with that friction (in
   `wallFriction_`, as a closed outlet). An inlet is open for the steps that start in
   that interval; its speed is taken at mid-step. Negative speeds stop the run (use an
   outlet to let material out).
2. **Lagrange step:** an open inlet face moves into the domain at the inlet speed,
   v* = -s n, so phiL = -s |S| (the Lagrangian cell is pushed in by a piston). The
   traction is the fed material's own, -p n, unless the cell's material is at the face
   and pressed against it: then the interior characteristic with v* (as a velocity
   patch), as long as that is more compressive than -p. No tension across the feed, so a
   falling stream detaches from the inlet.
3. **Remap:** the volume the face swept, dt s |S|, is fresh material: material volume
   in, mass dt s |S| rho_in, and the mass-weighted fields bring the inlet values
   (velocity -s n, stress -p I, no plastic strain) through the same increment formula as
   internal faces (`remapMassWeighted`, inlet values per patch face). The fed mass is
   counted in `charged()` (all processors), as `discharged()` counts the outflow;
   celFoam prints both with the totals.
4. **Slab side:** an open inlet counts as empty for the slab rule, as a traction-free
   patch. Counting it as full (fresh material beyond) was tried first: the cells of the
   stream then leaned up towards the inlet, their material never reached the face
   below, and the stream fell in place (2.9 m/s material speed, front moving at the
   inlet speed 0.5 m/s). As empty, the stream breaks into packets as it accelerates
   and lands.

## Consequences

- `tests/test_inlet.sh` (about 25 s, in CI): empty 0.2 x 0.3 m box, 10 mm, T0 rule
  set at E = 1 MPa, 0.04 m inlet in the lid at 0.5 m/s from 0.05 to 0.35 s. Nothing in
  while closed (0), fed mass = rho s w t to 1e-14, box mass = fed mass to 1e-14, no
  stress in empty cells, stream density = inlet density to 6e-4 (the nearly-full snap of
  decision 0045), deposit +0.82 % (limit 1 %), at rest (0.010 m/s) at 1.2 s, 2
  processors split through the inlet match serial to 1e-14.
- With plain Mohr-Coulomb (no flowImpedance) the landing packs the deposit to +4.3 %
  (overfilled cells, as decision 0053): use the T0 rule set when material lands.
- Without `inlets` every code path is unchanged (the remap's inlet arguments are
  null); `tests/test_void.sh` and `tests/test_outlet.sh` give the same numbers.
- The fed stream is a continuous column of whole cells at the inlet and packets
  below; its shape is set by the slab rules, not resolved. A real stream thins as it
  accelerates; a finer mesh or a narrower inlet does not change the fed rate, which is
  exact.
- Next (WP 3.7b): a silo fed from the top while it discharges, as a tutorial (first
  tutorial with `inlets`).
