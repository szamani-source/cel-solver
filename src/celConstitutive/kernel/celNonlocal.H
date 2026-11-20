/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::Nonlocal

Description
    Nonlocal granular fluidity (WP N.1, decision 0111; Kamrin and Koval
    2012, Henann and Kamrin 2013, 2014). Reference:
    prototype/constitutive/ngf.py.

        mu = tau/P, tau = sqrt(J2), P = -tr(sigma)/3
        gamma_dot_p = g mu                      (flow rule, g the fluidity)
        g_loc = sqrt(P/(rhoS d^2)) (mu - mu_s)/(b mu)   mu > mu_s, else 0
        k g - lap(g) = k g_loc,  k = |mu - mu_s|/(A d)^2

    b = (mu2 - mu_s)/I0 (the small-I slope of the mu(I) entries), A the
    entry nonlocalA. The field equation is solved by the Eulerian solver;
    this class gives g_loc, k and the point update at a given g:

        elastic trial, then the deviator scaled down at fixed P (plastic
        flow along the deviator, no volume change) so that
        tau = tau_trial - G dt g mu, with g = g_loc(mu_new) + (g - g_loc(old)):
        the local part implicit (an explicit g is unstable once
        G dt mu g'(mu)/P > 2), the nonlocal part lagged.

    P <= 0 (tension) falls back to the base model's rate-independent return.

SourceFiles
    celNonlocal.C
\*---------------------------------------------------------------------------*/

#ifndef celNonlocal_H
#define celNonlocal_H

#include "celConstitutiveKernel.H"

#include <memory>

namespace cel
{
namespace constitutive
{

class Nonlocal
{
    double muS_;
    double b_;
    double A_;
    double d_;
    double rhoS_;
    //- Cohesion as an added pressure, c cot(phi) (WP C.1, decision 0112):
    //  the law reads P* = P + pC everywhere it reads P
    double pC_;

public:

    Nonlocal
    (
        double muS, double b, double A, double d, double rhoS,
        double pC = 0
    );

    //- The law for a base model with friction angle phiDeg, or null if p
    //  has no nonlocalA (or 0). Needs mu2, I0, d, rhoS; throws otherwise.
    static std::unique_ptr<Nonlocal> New(const Parameters& p, double phiDeg);

    double muS() const { return muS_; }
    double b() const { return b_; }
    double A() const { return A_; }
    double d() const { return d_; }
    double pC() const { return pC_; }
    double rhoS() const { return rhoS_; }

    //- mu = sqrt(J2)/P* and P* = P + c cot(phi) of a stress (P* <= 0:
    //  mu = 0)
    void muP(const Voigt s, double& mu, double& P) const;

    //- Local fluidity [1/s] at mu and P* (as muP returns)
    double gLoc(double mu, double P) const;

    //- Coefficient k = |mu - mu_s|/(A d)^2 [1/m^2]
    double k(double mu) const;

    //- Point update at fluidity g over dt from the elastic trial: the
    //  nonlocal part g - g_loc(sigOld) is held, the local part follows the
    //  new stress. Returns false (nothing done) if the trial pressure is not
    //  positive.
    bool update
    (
        const Elastic& el,
        const Voigt sigOld,
        const Voigt trial,
        double dt,
        double g,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const;
};


} // End namespace constitutive
} // End namespace cel

#endif
