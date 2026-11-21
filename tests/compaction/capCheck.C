// WP K.1 (decision 0113): the compaction cap at one point. Mohr-Coulomb and
// Drucker-Prager with capLambda (lambda* = 0.01, E 50 MPa), from an
// isotropic 1 kPa: (a) isotropic compression to 100 kPa and back to 10 kPa,
// (b) oedometer (only eps_yy) to 100 kPa mean stress and back to 10 kPa.
// Prints per step: model test phase eps_v(compressive) P pc.
#include "celConstitutiveKernel.H"
#include "celCap.H"

#include <cmath>
#include <cstdio>

using namespace cel::constitutive;

int main()
{
    Parameters p
    {
        {"E", 5e7}, {"nu", 0.3}, {"c", 0}, {"phi", 30}, {"psi", 0},
        {"capLambda", 0.01}, {"capP0", 100}
    };
    for (const char* name : {"mohrCoulomb", "druckerPrager"})
    for (const char* test : {"isotropic", "oedometer"})
    {
        auto m = makeModel(name, p);
        const bool iso = test[0] == 'i';
        Voigt s = {-1000, -1000, -1000, 0, 0, 0};
        History h;
        double ev = 0;
        const double de = 2e-5;
        for (int phase = 0; phase < 2; ++phase)
        {
            const double sign = phase == 0 ? 1 : -1;
            for (int i = 0; i < 100000; ++i)
            {
                const double P = -(s[0] + s[1] + s[2])/3;
                if (phase == 0 && P >= 1e5) break;
                if (phase == 1 && P <= 1e4) break;
                Voigt d = {0, 0, 0, 0, 0, 0}, sn;
                if (iso)
                {
                    d[0] = d[1] = d[2] = -sign*de/3;
                }
                else
                {
                    d[1] = -sign*de;
                }
                History hn;
                m->updateRate(d, 1e-4, s, h, sn, hn);
                for (int k = 0; k < 6; ++k) s[k] = sn[k];
                h = hn;
                ev += sign*de;
                std::printf
                (
                    "%s %s %d %.15g %.15g %.15g\n", name, test, phase, ev,
                    -(s[0] + s[1] + s[2])/3, h.pc
                );
            }
        }
    }
    return 0;
}
