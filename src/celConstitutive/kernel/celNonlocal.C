#include "celNonlocal.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cel
{
namespace constitutive
{

Nonlocal::Nonlocal
(
    double muS, double b, double A, double d, double rhoS, double pC
)
:
    muS_(muS),
    b_(b),
    A_(A),
    d_(d),
    rhoS_(rhoS),
    pC_(pC)
{}


std::unique_ptr<Nonlocal> Nonlocal::New(const Parameters& p, double phiDeg)
{
    const double A = parameter(p, "nonlocalA", 0);
    if (A <= 0)
    {
        return nullptr;
    }
    const double pi = std::acos(-1.0);
    const double muS = std::tan(phiDeg*pi/180.0);
    const double mu2 = parameter(p, "mu2");
    const double I0 = parameter(p, "I0");
    const double d = parameter(p, "d");
    const double rhoS = parameter(p, "rhoS");
    if (!(mu2 > muS) || !(I0 > 0) || !(d > 0) || !(rhoS > 0))
    {
        throw std::invalid_argument
        (
            "nonlocal: needs mu2 > tan(phi), I0 > 0, d > 0, rhoS > 0"
        );
    }
    // Cohesion c (Mohr-Coulomb's c, Drucker-Prager's c as given): the
    // friction law acts on P* = P + c cot(phi), so yield at mu_s is the
    // cohesive line tau = c' + mu_s P (decision 0112)
    const double c = parameter(p, "c", 0);
    const double pC = (c > 0 && muS > 0) ? c/muS : 0;
    return std::unique_ptr<Nonlocal>
    (
        new Nonlocal(muS, (mu2 - muS)/I0, A, d, rhoS, pC)
    );
}


void Nonlocal::muP(const Voigt s, double& mu, double& P) const
{
    P = -meanStress(s) + pC_;
    // P so small that tau/P or rhoS/P would overflow (a free-surface cell,
    // as FrictionLaw::inertialNumber): no friction, no fluidity
    // (mu capped at 1e3, far above any admissible friction, so k(mu) stays
    // finite)
    mu = P > 1e-300*rhoS_ ? std::min(std::sqrt(J2(s))/P, 1e3) : 0;
}


double Nonlocal::gLoc(double mu, double P) const
{
    if (!(mu > muS_) || !(P > 1e-300*rhoS_))
    {
        return 0;
    }
    return std::sqrt(P/(rhoS_*d_*d_))*(mu - muS_)/(b_*mu);
}


double Nonlocal::k(double mu) const
{
    return std::abs(mu - muS_)/((A_*d_)*(A_*d_));
}


bool Nonlocal::update
(
    const Elastic& el,
    const Voigt sigOld,
    const Voigt trial,
    double dt,
    double g,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    const double P = -meanStress(trial) + pC_;
    if (!(P > 1e-300*rhoS_))
    {
        return false;
    }
    Voigt s;
    deviator(trial, s);
    const double tauT = std::sqrt(J2(s));
    if (tauT <= 0)
    {
        for (int i = 0; i < 6; ++i) sigNew[i] = trial[i];
        hNew = hOld;
        return true;
    }

    // Nonlocal part of the fluidity, lagged: delta = g - g_loc(old stress).
    // The local part follows the new stress (implicit), so a homogeneous
    // flow is the local law at any dt; g_loc mu = c (mu - mu_s)/b is linear
    // in mu above mu_s, and F(tau) below increases, so the root in
    // [0, tau_trial] is unique (Illinois method).
    double muO = 0, PO = 0;
    muP(sigOld, muO, PO);
    const double delta = g - gLoc(muO, PO);
    const double h = el.G*dt;
    auto F = [&](double tau)
    {
        const double mu = tau/P;
        return tau - tauT + h*std::max(gLoc(mu, P) + delta, 0.0)*mu;
    };
    double a = 0, b = tauT, fa = F(a), fb = F(b);
    double tau = b;
    if (fb > 0 && fa < 0)
    {
        int side = 0;
        for (int it = 0; it < 200; ++it)
        {
            tau = (a*fb - b*fa)/(fb - fa);
            const double ft = F(tau);
            if (std::abs(ft) <= 1e-14*tauT || b - a <= 1e-15*tauT)
            {
                break;
            }
            if (ft*fb > 0)
            {
                b = tau; fb = ft;
                if (side == -1) fa *= 0.5;
                side = -1;
            }
            else
            {
                a = tau; fa = ft;
                if (side == 1) fb *= 0.5;
                side = 1;
            }
        }
    }
    else if (fa >= 0)
    {
        tau = 0;
    }
    const double f = tau/tauT;
    for (int i = 0; i < 3; ++i)
    {
        sigNew[i] = pC_ - P + f*s[i];
        sigNew[i + 3] = f*s[i + 3];
    }
    plasticCorrector(el, trial, sigNew, hOld, false, hNew);
    return true;
}


} // End namespace constitutive
} // End namespace cel
