"""WP 3.6a: reference paths of the mu(I) models (prototype/constitutive/mu_i.py) for the
C++ kernel (tests/constitutive/compareReference.C).

Same file format as gen_constitutive_reference.py, plus one line per case

    dt <time step>

after the case line: the C++ side then replays the increments with updateRate(dt).

Usage: python3 gen_mu_i_reference.py <output file>
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "prototype", "constitutive"))

import driver  # noqa: E402
import models  # noqa: E402
import mu_i  # noqa: E402

PY_MODELS = {"druckerPrager": mu_i.mu_i_drucker_prager,
             "mohrCoulomb": mu_i.mu_i_mohr_coulomb}
CPP_PARAM = {"sigma_t": "sigmaT", "rho_s": "rhoS", "muI_alpha": "muIAlpha"}


def fmt(values):
    return " ".join("%.17g" % v for v in values)


def write_case(out, name, model_name, params, sigma0, deps_list):
    model = PY_MODELS[model_name]
    numbers = {k: v for k, v in params.items() if not isinstance(v, str) and k != "dt"}
    options = {k: v for k, v in params.items() if isinstance(v, str)}
    out.write("case %s %s %d %s\n" % (
        name, model_name, len(numbers),
        " ".join("%s %.17g" % (CPP_PARAM.get(k, k), v) for k, v in numbers.items())))
    out.write("dt %.17g\n" % params["dt"])
    for k, v in options.items():
        out.write("option %s %s\n" % (k, v))
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


def main():
    E, nu, pc = 20e6, 0.3, 100e3
    sand = dict(mu_i.sand_placeholder(), E=E, nu=nu, dt=1e-4)
    glass = dict(phi=20.9, mu2=math.tan(math.radians(32.76)), I0=0.279, d=0.53e-3,
                 rho_s=2450.0, muI_alpha=mu_i.ALPHA_BARKER_GRAY, E=E, nu=nu, dt=1e-4)
    sets = (
        ("muIdpSand", "druckerPrager", dict(sand, c=0.0, psi=0.0, fit="compression")),
        ("muIdpCohesive", "druckerPrager",
         dict(sand, c=10e3, psi=5.0, fit="plane_strain")),
        ("muImcSand", "mohrCoulomb", dict(sand, c=0.0, psi=0.0)),
        ("muImcTension", "mohrCoulomb", dict(sand, c=10e3, psi=5.0, sigma_t=5e3)),
        ("muIdpGlassBG", "druckerPrager", dict(glass, c=0.0, psi=0.0, fit="compression")),
        ("muImcGlassBG", "mohrCoulomb", dict(glass, c=0.0, psi=0.0)),
    )
    rng = np.random.default_rng(2036)
    random_steps = [10 ** rng.uniform(-6, -2.5) * rng.normal(size=6) for _ in range(300)]
    d = np.zeros(6)

    with open(sys.argv[1], "w") as out:
        for tag, name, p in sets:
            m = PY_MODELS[name]
            s0 = -pc * models.I_V
            # Simple shear at fixed pressure (all normal stresses held), shear
            # rates 2 and 20 1/s
            for gam, n in ((2e-4, 100), (2e-3, 60)):
                r = driver.run_path(m, p, s0, n, d + [0, 0, 0, gam, 0, 0],
                                    [False, False, False, True, True, True], s0)
                write_case(out, "%s_shearFixedP_%g" % (tag, gam), name, p, s0,
                           np.diff(r["eps"], axis=0))
            # Triaxial compression (mixed control)
            r = driver.triaxial(m, p, pc, -0.02, 100)
            write_case(out, tag + "_triaxialCompression", name, p, s0,
                       np.diff(r["eps"], axis=0))
            # Uniaxial strain: compress, then pull to the apex (P -> 0, I -> infinity)
            write_case(out, tag + "_uniaxialPull", name, p, s0,
                       [d + [0, 0, -2e-4, 0, 0, 0]] * 50 + [d + [0, 0, 2e-4, 0, 0, 0]] * 100)
            # Random increments (shear rates 1e-2 to 30 1/s). Not named "random": the
            # WP 2.5a check that Drucker-Prager random paths reach the apex is for
            # its 1000-step paths (uniaxialPull reaches the apex here)
            write_case(out, tag + "_noise300", name, p, s0, random_steps)


if __name__ == "__main__":
    main()
