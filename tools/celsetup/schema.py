"""Every celSetup input with its default, unit and a one-line description.

One list serves the YAML reader (defaults), `celSetup --template` (a
commented starting file), the GUI (forms) and the manual. Each entry:
(key, default, unit, description). A default of REQUIRED must be given; None means left out.

"""

REQUIRED = "required"

MATERIAL = [
    ("model", "mohrCoulomb", "", "mohrCoulomb or druckerPrager"),
    ("rho", REQUIRED, "kg/m^3", "bulk density of the material as placed"),
    ("E", REQUIRED, "Pa", "Young's modulus"),
    ("nu", 0.3, "", "Poisson's ratio"),
    ("phi", REQUIRED, "deg", "internal friction angle"),
    ("psi", 1.0, "deg", "dilation angle (Sia: low, 0.5 to 1 deg, for discharge)"),
    ("c", 0.0, "Pa", "cohesion"),
    ("sigmaT", "auto", "Pa", "tensile strength (MC); auto: 0 without cohesion, else c cot(phi)"),
    ("fit", "plane_strain", "", "Drucker-Prager fit to Mohr-Coulomb: plane_strain (use this for silos), "
     "compression or extension"),
    ("d", REQUIRED, "m", "grain diameter: mu(I) law, flow impedance and the mesh checker's rules"),
    ("rhoS", "auto", "kg/m^3", "grain (solid) density; auto: rho / 0.62"),
    ("mu2", "auto", "", "mu(I) friction at high I; auto: tan(phi) + 0.26 (Jop et al. 2005)"),
    ("I0", 0.279, "", "mu(I) reference inertial number (Jop et al. 2005)"),
    ("nonlocalA", None, "", "nonlocal fluidity amplitude A (mesh-independent shear bands, "
     "decision 0111); leave empty for the local law"),
    ("rules", "recommended", "", "scheme options: recommended (shearCap, spreadLayers, flowImpedance, "
     "normalCap, fullBlend all on, as in the validated silo runs) or none"),
]

WALLS_RUN = [
    ("endTime", 2.0, "s", "simulated time"),
    ("openTime", 1.0, "s", "outlet opens at this time (settling before it)"),
    ("writeInterval", 0.05, "s", "write (and checkpoint) every this much simulated time"),
    ("maxCo", 0.4, "", "fraction of the stable (CFL) time step used"),
    ("gravityRamp", 0.1, "s", "gravity ramps in over this time at the start"),
    ("damping", 100.0, "1/s", "velocity damping while settling (switched off when the outlet opens)"),
    ("nProcs", 1, "", "processor cores for the run"),
]

SECTIONS = {"material": MATERIAL, "run": WALLS_RUN}


def apply_defaults(cfg, section, entries):
    d = dict(cfg.get(section) or {})
    for key, default, _, _ in entries:
        if d.get(key) is None:
            if default == REQUIRED:
                raise ValueError(f"{section}: '{key}' is required ({_desc(entries, key)})")
            d[key] = default
    cfg[section] = d
    return d


def _desc(entries, key):
    for k, _, unit, text in entries:
        if k == key:
            return f"{text}, {unit}" if unit else text
    return ""
