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


def fmt(values):
    return " ".join("%.17g" % v for v in values)


def write_case(out, name, model_name, params, sigma0, deps_list):
    """Run the Python model along deps_list from sigma0 and write the record."""
    model = models.MODELS[model_name]
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
            m = models.MODELS[name]

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


if __name__ == "__main__":
    main()
