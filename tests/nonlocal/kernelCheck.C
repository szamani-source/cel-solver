// WP N.1 (decision 0111): homogeneous simple shear at fixed pressure with the
// nonlocal Mohr-Coulomb and Drucker-Prager kernels at g = g_loc (uniform, so
// the nonlocal term vanishes), dt 1e-6 and 1e-4 s (the local part is
// implicit). Prints model, dt, rate, mu at the end, mu_s + b I.
#include "celConstitutiveKernel.H"
#include "celNonlocal.H"

#include <cmath>
#include <cstdio>

using namespace cel::constitutive;

int main()
{
    const double pi = std::acos(-1.0);
    Parameters p
    {
        {"E", 1e7}, {"nu", 0.3}, {"c", 0}, {"phi", std::atan(0.3819)*180/pi},
        {"psi", 0}, {"mu2", 0.6435}, {"I0", 0.278}, {"nonlocalA", 0.48},
        {"d", 1e-3}, {"rhoS", 2450}
    };
    const double P = 1000;
    for (const char* name : {"mohrCoulomb", "druckerPrager"})
    {
        auto m = makeModel(name, p);
        const Nonlocal& nl = *m->nonlocal();
        for (double dt : {1e-6, 1e-4})
        for (double rate : {1.0, 10.0, 50.0})
        {
            Voigt s = {-P, -P, -P, 0, 0, 0};
            History h;
            for (int i = 0; i < int(0.2/dt); ++i)
            {
                double mu, Pc;
                nl.muP(s, mu, Pc);
                const double g = nl.gLoc(mu, Pc);
                Voigt d = {0, 0, 0, rate*dt, 0, 0}, sn;
                History hn;
                m->updateFluidity(d, dt, g, s, h, sn, hn);
                for (int k = 0; k < 6; ++k) s[k] = sn[k];
                h = hn;
            }
            const double I = rate*nl.d()/std::sqrt(P/2450.0);
            std::printf("%s %g %g %.15g %.15g\n", name, dt, rate, s[3]/P, nl.muS() + nl.b()*I);
        }
    }
    return 0;
}
