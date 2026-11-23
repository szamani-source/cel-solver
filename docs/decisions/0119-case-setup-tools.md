# 0119: Case setup tools (celSetup, celCheck, GUI)

Date: 2026-11-23. Thread "Build case generator and mesh checker".

## Context

Sia: setting up simulations, above all a good mesh, is hard. The good Couto setup
(decisions 0105, 0107, 0115) lives in test scripts with the geometry hard-coded. Sia
approved a parametric case generator and a pre-run mesh checker, then added full
property setting with defaults and descriptions and a GUI with a live view, all before
the Linux installation package.

## Decision

- Python, numpy only, in `tools/celsetup/`; commands `tools/celSetup` and (S.3)
  `tools/celCheck`. The CI image and a bare OpenFOAM install have numpy but not PyYAML
  or matplotlib, so the input reader is our own small YAML subset (`miniyaml.py`,
  used always: PyYAML would read `1e-3` as a string) and pictures are optional.
- Input is YAML: geometry, walls, mesh, material, run. Templates turn geometry into
  blockMesh; the first, `conicalHopper`, is `aligned_buffered_mesh` with the Couto
  numbers turned into inputs. With the Couto inputs it reproduces the validated mesh
  and STL exactly (tested), so the generator inherits that validation.
- The material fill is computed by celSetup on the mesh (sub-sampled cells, alpha as
  the fraction of the open part outside the STL shell), not by setFields: the fill
  follows the rounded wall exactly and the mass placed is reported against the
  geometry.
- `startFrom latestTime` and an Allrun that skips meshing and decomposition when done:
  running Allrun again resumes from the last write (checkpoint rule for long runs).
- GUI (S.4): a local web page served by a small Python program, opened in the normal
  browser. Works offline, nothing to install beyond Python and numpy; Tk/Qt are often
  missing or differ between Linux distributions.

## Consequences

- Test scripts and celSetup share geometry by copy for now; the test checks they stay
  identical. Later templates (S.5) and the axis-column change (S.6) are options, so the
  validated mesh stays reachable.
