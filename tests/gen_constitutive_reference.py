"""WP 2.2: write reference strain paths from the Python constitutive models (WP 1.4).

The C++ kernel test (tests/constitutive/compareReference.C) replays the same strain
increments and compares stress and history step by step.

Usage: python3 gen_constitutive_reference.py <output file>

File format, one record per line, numbers written with 17 significant digits:

    case <name> <model> <n> <key> <value> ... (n parameter pairs)
    init <sigma 6>
    step <deps 6> <sigma 6> <epsP 6> <epsPEq> <apex 0|1>
    ...
    end

Voigt order xx, yy, zz, xy, yz, zx; strain with engineering shear (decision 0004).
<model> is the C++ model name (PY_MODELS maps it to the Python function).
Mixed stress/strain paths are recorded as the strain increments the Python driver found,
so the C++ side only needs a strain-driven replay.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "prototype", "constitutive"))

import driver  # noqa: E402
import models  # noqa: E402


# C++ model name -> Python reference model
PY_MODELS = {
    "linearElastic": models.linear_elastic,
    "vonMises": models.von_mises,
    "mohrCoulomb": models.mohr_coulomb,
}


def fmt(values):
    return " ".join("%.17g" % v for v in values)


def write_case(out, name, model_name, params, sigma0, deps_list):
    """Run the Python model along deps_list from sigma0 and write the record."""
    model = PY_MODELS[model_name]
    out.write("case %s %s %d %s\n" % (
        name, model_name, len(params),
        " ".join("%s %.17g" % (k, v) for k, v in params.items())))
    sigma = np.array(sigma0, dtype=float)
    history = models.init_history()
    out.write("init %s\n" % fmt(sigma))
    for deps in deps_list:
        deps = np.asarray(deps, dtype=float)
        sigma, history = model(deps, sigma, history, params)
        out.write("step %s %s %s %.17g %d\n" % (
            fmt(deps), fmt(sigma), fmt(history["eps_p"]), history["eps_p_eq"],
            1 if history["apex"] else 0))
    out.write("end\n")


def increments(run):
    """Strain increments applied by a driver path."""
    return np.diff(run["eps"], axis=0)


def main():
    E, nu = 20e6, 0.3
    elastic = dict(E=E, nu=nu)
    j2_perfect = dict(E=E, nu=nu, sigma_y=200e3, H=0.0)
    j2_hard = dict(E=E, nu=nu, sigma_y=200e3, H=2e6)
    pc = 100e3

    rng = np.random.default_rng(2022)

    def random_path(n):
        steps = []
        for _ in range(n):
            mag = 10 ** rng.uniform(-6, -2.5)
            steps.append(mag * rng.normal(size=6))
        return steps

    random_steps = random_path(1000)

    with open(sys.argv[1], "w") as out:
        for tag, name, p in (("elastic", "linearElastic", elastic),
                             ("j2perfect", "vonMises", j2_perfect),
                             ("j2hard", "vonMises", j2_hard)):
            m = PY_MODELS[name]

            # Uniaxial strain: load to 2 % compression, unload to 1 %, reload in tension
            d = np.zeros(6)
            load = [d + [0, 0, -2e-2 / 100, 0, 0, 0]] * 100
            unload = [d + [0, 0, 1e-2 / 50, 0, 0, 0]] * 50
            write_case(out, tag + "_uniaxialStrain", name, p, np.zeros(6),
                       load + unload + load[:30])

            # Triaxial compression and extension from 100 kPa (mixed control)
            r = driver.triaxial(m, p, pc, -0.02, 200)
            write_case(out, tag + "_triaxialCompression", name, p, -pc * models.I_V,
                       increments(r))
            r = driver.triaxial(m, p, pc, 0.02, 200, extension=True)
            write_case(out, tag + "_triaxialExtension", name, p, -pc * models.I_V,
                       increments(r))

            # Simple shear at constant normal stress (mixed control)
            r = driver.simple_shear(m, p, pc, 0.05, 200)
            write_case(out, tag + "_simpleShear", name, p, -pc * models.I_V,
                       increments(r))

            # 1000 random increments from an isotropic compressive state
            write_case(out, tag + "_random1000", name, p, -pc * models.I_V,
                       random_steps)

        # WP 2.5b: Mohr-Coulomb, frictional-cohesive and cohesionless sand
        mc = PY_MODELS["mohrCoulomb"]
        mc_sets = (("mcCohesive", dict(E=E, nu=nu, phi=30.0, c=10e3, psi=5.0)),
                   ("mcSand", dict(E=E, nu=nu, phi=35.0, c=0.0, psi=0.0)))
        for tag, p in mc_sets:
            # Uniaxial strain from 100 kPa: load in compression, unload, then pull
            # into tension until the apex is reached
            d = np.zeros(6)
            load = [d + [0, 0, -2e-2 / 100, 0, 0, 0]] * 100
            unload = [d + [0, 0, 1e-2 / 50, 0, 0, 0]] * 50
            pull = [d + [0, 0, 2e-2 / 100, 0, 0, 0]] * 100
            write_case(out, tag + "_uniaxialStrain", "mohrCoulomb", p,
                       -pc * models.I_V, load + unload + pull)

            # Triaxial compression and extension: edge returns (s1 = s2, s2 = s3)
            r = driver.triaxial(mc, p, pc, -0.02, 200)
            write_case(out, tag + "_triaxialCompression", "mohrCoulomb", p,
                       -pc * models.I_V, increments(r))
            r = driver.triaxial(mc, p, pc, 0.02, 200, extension=True)
            write_case(out, tag + "_triaxialExtension", "mohrCoulomb", p,
                       -pc * models.I_V, increments(r))

            # Simple shear at constant normal stress
            r = driver.simple_shear(mc, p, pc, 0.05, 200)
            write_case(out, tag + "_simpleShear", "mohrCoulomb", p,
                       -pc * models.I_V, increments(r))

            # 1000 random increments (same fixed-seed increments as above)
            write_case(out, tag + "_random1000", "mohrCoulomb", p,
                       -pc * models.I_V, random_steps)


if __name__ == "__main__":
    main()
