#include "celVonMises.H"

#include <cmath>
#include <stdexcept>

namespace cel
{
namespace constitutive
{

VonMises::VonMises(const Parameters& p)
:
    el_(parameter(p, "E"), parameter(p, "nu")),
    sigmaY_(parameter(p, "sigma_y")),
    H_(parameter(p, "H", 0.0))
{
    if (sigmaY_ <= 0)
    {
        throw std::invalid_argument("vonMises: sigma_y must be positive");
    }
}


void VonMises::update
(
    const Voigt deps,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    const double G = el_.G;

    Voigt trial;
    el_.addStress(sigOld, deps, trial);

    hNew = hOld;
    hNew.apex = false;

    const double q = std::sqrt(3.0*J2(trial));
    const double sy = sigmaY_ + H_*hOld.epsPEq;
    const double f = q - sy;

    if (f <= yieldTol*sy)
    {
        for (int i = 0; i < 6; ++i)
        {
            sigNew[i] = trial[i];
        }
        return;
    }

    // Radial return
    const double dgam = f/(3.0*G + H_);
    const double p = meanStress(trial);
    const double factor = 1.0 - 3.0*G*dgam/q;

    Voigt s;
    deviator(trial, s);
    for (int i = 0; i < 6; ++i)
    {
        sigNew[i] = s[i]*factor + (i < 3 ? p : 0.0);
    }

    // Plastic strain from the corrector: deps_p = C (trial - sigNew)
    Voigt dsig, depsP;
    for (int i = 0; i < 6; ++i)
    {
        dsig[i] = trial[i] - sigNew[i];
    }
    el_.strain(dsig, depsP);
    for (int i = 0; i < 6; ++i)
    {
        hNew.epsP[i] = hOld.epsP[i] + depsP[i];
    }

    // Exact equivalent plastic strain increment for J2
    hNew.epsPEq = hOld.epsPEq + dgam;
}


double VonMises::yield(const Voigt sigma, const History& h) const
{
    return std::sqrt(3.0*J2(sigma)) - sigmaY_ - H_*h.epsPEq;
}

} // End namespace constitutive
} // End namespace cel
