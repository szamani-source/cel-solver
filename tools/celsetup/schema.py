"""Every celSetup input with its default, unit and a one-line description
(WP S.2).

One list serves the YAML reader (defaults), `celSetup --template` (a
commented starting file), `celSetup --describe`, the GUI (forms) and the
manual. Each entry is E(key, default, unit, description, choices, level):
- default REQUIRED must be given; None means left out (off);
- choices, when given, are the allowed words;
- level "basic" is shown by default, "advanced" in the GUI's advanced
  part and commented out in a template file.

The material and run sections are shared by all templates; geometry,
walls and mesh come from the template module (its SCHEMA).
"""

from collections import namedtuple

REQUIRED = "required"
AUTO = "auto"

E = namedtuple("E", "key default unit text choices level")
E.__new__.__defaults__ = (None, "basic")

MATERIAL = [
    E("model", "mohrCoulomb", "", "material model", ["mohrCoulomb", "druckerPrager"]),
    E("rho", REQUIRED, "kg/m^3", "bulk density of the material as placed"),
    E("E", REQUIRED, "Pa", "Young's modulus (sets the time step: stiffer is slower)"),
    E("nu", 0.3, "", "Poisson's ratio"),
    E("phi", REQUIRED, "deg", "internal friction angle"),
    E("psi", 1.0, "deg", "dilation angle; keep it low (0.5 to 1 deg) for discharge"),
    E("c", 0.0, "Pa", "cohesion (0 for dry sand and grain)"),
    E("sigmaT", AUTO, "Pa", "tensile strength, Mohr-Coulomb only; auto: 0 without cohesion, "
      "else c cot(phi)"),
    E("fit", "plane_strain", "", "how the Drucker-Prager cone matches Mohr-Coulomb; plane_strain "
      "for silos (compression gave far too slow a discharge)",
      ["plane_strain", "compression", "extension"]),
    E("d", REQUIRED, "m", "grain diameter: mu(I) law, flow impedance and the mesh checker's rules"),
    E("rhoS", AUTO, "kg/m^3", "density of the grains themselves; auto: rho / 0.62"),
    E("mu2", AUTO, "", "mu(I) friction coefficient at high inertial number, above tan(phi); "
      "auto: tan(phi) + 0.26 (Jop et al. 2005)"),
    E("I0", 0.279, "", "mu(I) reference inertial number (0.279 for glass beads, Jop et al. 2005)"),
    E("nonlocalA", None, "", "nonlocal fluidity amplitude A (0.48 in the Couto runs): shear bands "
      "a few grains wide whatever the mesh (decision 0111); empty for the local law"),
    E("capLambda", None, "", "compaction cap: modified compression index lambda* from an "
      "oedometer test (decision 0113); empty for no cap", level="advanced"),
    E("capP0", 100.0, "Pa", "pressure at which the cap starts (with capLambda)", level="advanced"),
    E("muIAlpha", None, "", "Barker-Gray regularisation of mu(I), 0 to 2; empty for off",
      level="advanced"),
    E("rules", "recommended", "", "scheme options: recommended turns on the five switches below, "
      "as in the validated silo runs; none turns them off (each can still be set)",
      ["recommended", "none"]),
    E("shearCap", AUTO, "", "cap the face shear stress at yield; capped faces slip "
      "(auto: follows rules)", ["auto", "yes", "no"], "advanced"),
    E("normalCap", AUTO, "", "cap the face normal stress jump at yield (auto: follows rules)",
      ["auto", "yes", "no"], "advanced"),
    E("spreadLayers", AUTO, "", "partly filled top cells spread along the surface "
      "(auto: follows rules)", ["auto", "yes", "no"], "advanced"),
    E("flowImpedance", AUTO, "", "flowing faces use the mu(I) viscosity instead of the elastic "
      "impedance (auto: follows rules)", ["auto", "yes", "no"], "advanced"),
    E("fullBlend", AUTO, "", "nearly full cells blend smoothly into the partly filled rules; "
      "removes speed spikes (auto: follows rules)", ["auto", "yes", "no"], "advanced"),
    E("voidPlasticFlow", "no", "", "plastic expansion of a partly filled cell goes into its void "
      "(needed for strongly dilating materials, decision 0095)", ["yes", "no"], "advanced"),
]

SWITCHES = ["shearCap", "normalCap", "spreadLayers", "flowImpedance", "fullBlend"]

RUN = [
    E("endTime", 2.0, "s", "simulated time"),
    E("openTime", 1.0, "s", "the outlet opens at this time; the material settles before it"),
    E("writeInterval", 0.05, "s", "write results (and a restart point) every this much simulated time"),
    E("maxCo", 0.4, "", "fraction of the stable time step used (0.4 is tested)"),
    E("gravityRamp", 0.1, "s", "gravity grows from zero over this time, no shock at the start"),
    E("damping", 100.0, "1/s", "velocity damping while settling; switched off when the outlet opens"),
    E("dampingStart", 0.0, "s", "damping acts from this time", level="advanced"),
    E("nProcs", 1, "", "processor cores for the run (a few thousand cells per core or more)"),
    E("gravity", 9.81, "m/s^2", "gravity, along -y", level="advanced"),
]

SHARED = {"material": MATERIAL, "run": RUN}
SECTION_TEXT = {
    "geometry": "Shape and size",
    "walls": "Silo wall (an STL shell inside the mesh) and its friction",
    "mesh": "Mesh: cell sizes and grading (the rules of the manual's meshing guide)",
    "material": "Material",
    "run": "Run control",
    # template-specific texts: "template.section"
    "flatBottom.walls": "Silo walls (the mesh sides and floor) and their friction",
    "flatBottom.mesh": "Mesh: cell sizes and grading (fine at the outlet, growing away from it)",
    "userMesh.mesh": "Your mesh: the file and its format (2D one cell thick, or a wedge about y)",
    "userMesh.geometry": "Where the material starts: the cells inside this box",
}


def section_text(sec, name=None):
    return SECTION_TEXT.get(f"{name}.{sec}", SECTION_TEXT.get(sec, ""))


OMIT = object()


def override(entries, **changes):
    """A copy of a shared list with some defaults changed; OMIT drops the
    entry (for templates where it means nothing, e.g. openTime without an
    outlet)."""
    out = []
    for e in entries:
        if e.key in changes:
            v = changes[e.key]
            if v is OMIT:
                continue
            e = e._replace(default=v)
        out.append(e)
    return out


def sections(template_schema):
    """Ordered {section: entries} for a template: its own sections, then the
    shared material and run lists unless the template gives its own."""
    out = dict(template_schema)
    for k, v in SHARED.items():
        out.setdefault(k, v)
    return out


def _word(v):
    if v is True:
        return "yes"
    if v is False:
        return "no"
    return v


def apply_defaults(given, section, entries):
    d = dict(given.get(section) or {})
    for e in entries:
        if d.get(e.key) is None:
            if e.default == REQUIRED:
                raise ValueError(f"{section}: '{e.key}' is required ({describe(e)})")
            d[e.key] = e.default
        elif e.choices and str(_word(d[e.key])) not in e.choices:
            raise ValueError(f"{section}: '{e.key}' must be one of {', '.join(e.choices)} "
                             f"(got {d[e.key]!r})")
    return d


def describe(e):
    return f"{e.text}, {e.unit}" if e.unit else e.text


def _yaml_value(v):
    if v is None:
        return ""
    return str(_word(v))


def template_yaml(name, template_schema, values=None, advanced=True, boundary_types=None):
    """A commented YAML file with every entry: given values, else defaults.
    Required entries without a value are left for the user to fill in;
    advanced entries not given are commented out."""
    values = values or {}
    lines = [f"# celSetup input, template {name}. Every entry is listed with its",
             "# default; edit what you need and run: celSetup THIS_FILE CASE_DIR",
             "# REQUIRED entries must be filled in. Units are SI. Lines starting",
             "# with # are advanced options at their defaults.",
             f"case: {values.get('case', name + 'Case')}", f"template: {name}", ""]
    for sec, entries in sections(template_schema).items():
        lines.append(f"{sec + ':':<27}# {section_text(sec, name)}")
        given = values.get(sec) or {}
        for e in entries:
            v = given.get(e.key, e.default)
            unit = f" [{e.unit}]" if e.unit else ""
            ch = f" ({' | '.join(e.choices)})" if e.choices else ""
            text = f"{e.text}{unit}{ch}"
            if v == REQUIRED:
                lines.append(f"  {e.key + ':':<24} # REQUIRED: {text}")
            elif e.level == "advanced" and e.key not in given:
                if advanced:
                    lines.append(f"  # {(e.key + ': ' + _yaml_value(v)):<22} # {text}")
            else:
                lines.append(f"  {(e.key + ': ' + _yaml_value(v)):<24} # {text}")
        lines.append("")
    if boundary_types is not None:
        lines += boundaries_yaml(values.get("boundaries"), boundary_types)
    return "\n".join(lines)


def boundaries_yaml(given, btypes):
    """The boundaries section: one block per patch, or (none given yet) a
    commented list of the conditions and their entries."""
    out = [f"{'boundaries:':<27}# one entry per patch of the mesh: type and its entries"]
    if given:
        for p, d in given.items():
            d = d if isinstance(d, dict) else {"type": d}
            out.append(f"  {p}:")
            t = d.get("type")
            out.append(f"    {('type: ' + str(t)):<22} # {btypes[t][0] if t in btypes else ''}")
            for e in (btypes[t][1] if t in btypes else []):
                if d.get(e.key) not in (None, ""):
                    unit = f" [{e.unit}]" if e.unit else ""
                    out.append(f"    {(e.key + ': ' + _yaml_value(d[e.key])):<22} # {e.text}{unit}")
    else:
        out.append("  # patchName:")
        out.append("  #   type: wall            # then the entries of that type:")
    for t, (text, entries) in btypes.items():
        ent = ", ".join(e.key + ("" if e.default != REQUIRED else " (required)") for e in entries)
        out.append(f"  # {t}: {text}" + (f"; entries: {ent}" if ent else ""))
    out.append("")
    return out


def table(template_schema, name=None):
    """Plain-text listing for celSetup --describe."""
    out = []
    for sec, entries in sections(template_schema).items():
        out.append(f"\n{sec}: {section_text(sec, name)}")
        for e in entries:
            d = "REQUIRED" if e.default == REQUIRED else (_yaml_value(e.default) or "(off)")
            out.append(f"  {e.key:<18} {d:<14} {e.unit:<8} {e.text}")
    return "\n".join(out)


def to_json(template_schema):
    """The schema for the GUI."""
    return {sec: [dict(e._asdict(), section=sec) for e in entries]
            for sec, entries in sections(template_schema).items()}
