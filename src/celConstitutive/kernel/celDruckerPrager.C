#include "celDruckerPrager.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cel
{
namespace constitutive
{

namespace
{

//- (eta, xi) for one angle and fit, as coeffs() in dp_coefficients (models.py)
void coefficients
(
    DruckerPrager::Fit fit,
    double angleDeg,
    double& eta,
    double& xi
)
{
    const double a = angleDeg*M_PI/180.0;
    const double sa = std::sin(a), ca = std::cos(a), ta = std::tan(a);

    switch (fit)
    {
        case DruckerPrager::Fit::compression:
        {
            const double d = std::sqrt(3.0)*(3.0 - sa);
            eta = 6.0*sa/d;
            xi = 6.0*ca/d;
            break;
        }
        case DruckerPrager::Fit::extension:
        {
            const double d = std::sqrt(3.0)*(3.0 + sa);
            eta = 6.0*sa/d;
            xi = 6.0*ca/d;
            break;
        }
        case DruckerPrager::Fit::planeStrain:
        {
            const double d = std::sqrt(9.0 + 12.0*ta*ta);
            eta = 3.0*ta/d;
            xi = 3.0/d;
            break;
        }
    }
}

} // End anonymous namespace


DruckerPrager::Fit DruckerPrager::fitFromName(const std::string& name)
{
    if (name == "compression") return Fit::compression;
    if (name == "extension") return Fit::extension;
    if (name == "plane_strain") return Fit::planeStrain;
    throw std::invalid_argument
    (
        "druckerPrager: unknown fit '" + name
      + "' (valid: compression extension plane_strain)"
    );
}


const char* DruckerPrager::fitName(Fit fit)
{
    switch (fit)
    {
        case Fit::compression: return "compression";
        case Fit::extension: return "extension";
        case Fit::planeStrain: return "plane_strain";
    }
    return "?";
}


DruckerPrager::DruckerPrager(const Parameters& p, const Options& opts)
:
    el_(parameter(p, "E"), parameter(p, "nu")),
    c_(parameter(p, "c")),
    fit_(fitFromName(option(opts, "fit", "compression"))),
    eta_(0),
    xi_(0),
    etaBar_(0)
{
    const double phi = parameter(p, "phi");
    const double psi = parameter(p, "psi");

    if (c_ < 0)
    {
        throw std::invalid_argument("druckerPrager: c must not be negative");
    }
    if (phi < 0 || phi >= 90 || psi < 0 || psi > phi)
    {
        throw std::invalid_argument
        (
            "druckerPrager: need 0 <= psi <= phi < 90 degrees"
        );
    }
    if (phi == 0 && c_ == 0)
    {
        throw std::invalid_argument
        (
            "druckerPrager: phi = 0 needs c > 0 (no strength otherwise)"
        );
    }

    double unused;
    coefficients(fit_, phi, eta_, xi_);
    coefficients(fit_, psi, etaBar_, unused);
}


void DruckerPrager::update
(
    const Voigt deps,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    const double K = el_.K, G = el_.G;

    Voigt trial;
    el_.addStress(sigOld, deps, trial);

    const double p = meanStress(trial);
    Voigt s;
    deviator(trial, s);
    const double sqj = std::sqrt(J2(trial));
    const double f = sqj + eta_*p - xi_*c_;

    if (f <= yieldTol*std::max({std::abs(p), sqj, c_, 1.0}))
    {
        hNew = hOld;
        hNew.apex = false;
        std::copy(trial, trial + 6, sigNew);
        return;
    }

    const double dgam = f/(G + K*eta_*etaBar_);

    if (sqj - G*dgam >= 0.0 || eta_ == 0.0)
    {
        // Return to the smooth cone
        const double factor = sqj > 0.0 ? 1.0 - G*dgam/sqj : 1.0;
        const double pNew = p - K*etaBar_*dgam;
        for (int i = 0; i < 6; ++i)
        {
            sigNew[i] = s[i]*factor + (i < 3 ? pNew : 0.0);
        }
        plasticCorrector(el_, trial, sigNew, hOld, false, hNew);
        return;
    }

    // Apex: zero deviator, mean stress at the cone tip (zero if cohesionless)
    const double pApex = xi_*c_/eta_;
    for (int i = 0; i < 6; ++i)
    {
        sigNew[i] = i < 3 ? pApex : 0.0;
    }
    plasticCorrector(el_, trial, sigNew, hOld, true, hNew);
}


double DruckerPrager::yield(const Voigt sigma) const
{
    return std::sqrt(J2(sigma)) + eta_*meanStress(sigma) - xi_*c_;
}

} // End namespace constitutive
} // End namespace cel
