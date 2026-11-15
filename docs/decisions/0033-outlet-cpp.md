# 0033 Outlet in C++: a wall that opens (WP 3.3)

Date: 2026-11-15. Status: accepted (WP 3.3). Copies the prototype outlet (rule 10,
decision 0026) with its fix (rule 10b, decision 0027) into the C++ cycle with void
(decision 0031), shear cap and rule 9 (decision 0032); neither rule 12 nor 13 (decision
0029).

## Context

T0 settles with the outlet closed, then opens it. The free-surface remap of decision 0031
already lets material out through any non-coupled face whose face velocity points out
(slab overlap, donor's Lagrangian density), and a traction patch already gives the
prototype's outlet face: zero traction, free one-sided velocity, empty ghost for the
slab sides. What was missing: a patch that is a wall first and a traction face later,
rule 10b, and a count of the mass that left.

## Decision

1. **`outlets` in `materialProperties`** (needs `freeSurface yes`):
   `outlets { outlet { friction 20; openTime 0.3; } }`. A closed outlet is a wall with
   that friction (it is put in the wall table, so the Lagrange step, the Hancock
   predictor and the slab sides treat it as one). It opens at the start of the first
   step at or after `openTime` (default never; `eulerianSolid::openOutlets()` opens all
   at once): it leaves the wall table and becomes a traction patch. Its velocity
   boundary condition must not fix a value.
2. **Opening ends settling: damping is set to 0**, as the prototype's discharge stage.
   One switch instead of a second time; settling without an outlet is unchanged.
3. **Rule 10b**: at an open outlet face whose cell holds material that does not touch the
   face (slab leaning away), the face takes the free one-sided velocity of the cell
   instead of 0, so the void under the slab passes out first. Only outlets: other
   traction faces keep velocity 0 there, so earlier cases run bit for bit.
4. **`discharged()`**: mass that left through non-coupled boundary faces, summed over
   processors (only open outlets let mass out; walls give no outward face velocity).
   `celVoidCheck` writes it as the last column of `history.dat`.
5. **Test** `tests/test_outlet.sh` (`tests/outlet/outlet_checks.py`, about 30 s): the
   prototype's outlet check case at 10 mm (W 0.4 m, H 0.6 m, D 0.06 m, E = 0.1 MPa,
   friction 20 deg, ramp 0.1 s, damping 250 1/s, opens at 0.3 s), outlet as its own
   patch (three blocks).

## Result (10 mm, E = 0.1 MPa, 22 s for 1.1 s of simulated time)

| Check | C++ | Prototype (no rules 12, 13) |
| --- | --- | --- |
| Mass through the closed outlet | 0 exactly; tank mass 4e-14 | 0 |
| First material out | 0.7 ms after opening | |
| Tank plus discharged | constant to 4e-14 | 6e-15 |
| Stress in empty cells | 0 | |
| Rate, 97 to 88 % of the mass | 74.1 kg/s/m, C 0.89 (t 0.24 to 0.77 s) | 81.35, C 0.98 (t 0.20 to 0.69 s) |
| Same without rule 10b | 27.4 kg/s/m, C 0.33 | about 0.3 (decision 0027) |

The C++ rate is 9 percent below the prototype's (test: within 20 percent, C at least
0.8). Its window starts 0.04 s later, so the start-up is slower too. The C++ column
collapse also runs a little shorter than the prototype (0.35 against 0.40 at 1 MPa,
decision 0032); the same differences (unsplit faces, 3D cap) are the likely cause. Not
chased here.

## Consequences

- WP 3.4 measures the rate at E = 0.1, 1, 10 MPa with this outlet and the collapse code;
  a rule 12 / 13 is added only if the rate depends on E (decision 0029).
- The outlet works per face, so it runs in parallel with the outlet split across
  processors; not tested here (WP 3.4 has the parallel check).
- Damping cannot stay on after an outlet opens. If a case ever needs that, make it a
  separate switch.
