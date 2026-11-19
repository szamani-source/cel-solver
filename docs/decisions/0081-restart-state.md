# 0081 Restart: state outside the fields, and a restart-proof time step (WP 6.2a)

Date: 2026-11-18. Status: accepted.

## Context

WP 6.2's done-when is "restarted run matches continuous run". celFoam already read
rho, U, sigma, epsP, epsPEq and alpha at the start time, but a restart lost:

- the bodies: free-body velocity and pose, the elements of deformable bodies (node
  positions and velocities, stresses, history, hourglass forces) and their surface;
  free bodies even refused to start after t = 0;
- the mass counters: discharged, fed, mass into each body, hand-over loss;
- the wall face velocities of the last step, which the flow impedance's wall ghost
  uses (decision 0044);
- the time steps: with `writeControl adjustableRunTime` OpenFOAM shortens each step to
  land on the next write time, counted from the start time; `0.1 - (t - 0.4)` and
  `0.5 - t` differ in the last digit, and the scheme's switches turn that into a
  different run (decision 0080).

## Decision

1. **State files.** At each write time, every processor writes
   `<time>/uniform/eulerianSolid` (counters and the wall face velocities per patch) and
   `<time>/uniform/rigidBodies` (each body's time, centre, rotation, velocity, omega;
   for deformable bodies also the surface points, surface velocity, element time and
   the element state), ASCII with 17 digits. Both are read when a run starts from that
   time. A free or deformable body without its state at t > 0 stops with an error;
   prescribed bodies start from their path as before.
2. **Time step.** celFoam computes the step itself: the stable step; with
   `writeControl adjustableRunTime` shortened (never lengthened) so that the steps land
   on the next write time, counted from time zero, and on endTime. Other write controls
   take the stable step unchanged, as before. OpenFOAM's own adjustment is off (`setDeltaT(dt, false)`); its write
   trigger still fires at those times. OpenFOAM rounds the number of steps to the
   nearest, which lengthened steps by up to half past the stable step; now it is
   rounded up.
3. **Output files.** postProcessing `bodyForces.dat` and `rigidBodies/<name>.dat` are
   appended to on a restart, not overwritten. bodyForces.dat now has 10 digits.
4. **Exactness needs exact fields**: `writeFormat binary` (or `writePrecision 17`). The
   test uses binary.

`tests/test_restart.sh` (about 40 s): the T0 tank restarted at 0.4 s while it
discharges, and a tank with a free steel block and a hex8 footing both in contact,
restarted at 0.05 s, serial and on 2 processors: fields, uniform/ state and
bodyForces.dat after the restart identical bit for bit to the continuous run.

## Consequences

- Restarting from a time that was not written by celFoam (or with a changed
  writeInterval) still works but lands on write times counted from zero.
- celFoam runs with adjustableRunTime (the tutorials) now take slightly different steps
  than before (rounded up, counted from zero); results move by round-off, which the switches can amplify to
  the scatter of decision 0046 (tutorial numbers are inside it).
- The restart state of a deformable body grows with its elements (about 1 kB per
  element in ASCII); fine for the body sizes in use.
- Mass scaling (the other half of the old WP 6.2) is split to WP 6.2b.
