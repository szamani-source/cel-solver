# WP 3.6f: the 1 MPa dip in the silo discharge rate

State: closed, finding recorded (2026-11-18)
Branch: `claude/project-thread-8hqs9j`

## Goal

Find why the T0 tank with mu(I) drains about 20 % slower at 1 MPa than at 0.1 and
10 MPa (decision 0045), using several runs per point to handle the scatter.

## Read first

Decisions 0029, 0044, 0045, 0046; `tests/discharge/scatter.py`.

## Done when

Each E measured as a band of three runs, and the cause of the dip named with runs that
switch it on and off.

## Log

- 2026-11-18: One session. Bands (three runs, maxCo 0.40 / 0.38 / 0.36): 65.8 / 52.0 /
  60.8 kg/s/m. Without the nearly-full snap: 65.5 / 60.5 / 63.1 (all bands within 5 %),
  but the inclined layer then fails; tolerances 1e-2 and 1e-7, and counting nearly-full
  cells as full in the Lagrange rules only, do not serve both. The rate is set by how
  loosened cells above the outlet are treated (full or partly filled), so the void
  closure is the lever, not the friction law. Proposed cure: compressible mu(I)/phi(I).
