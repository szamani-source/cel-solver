/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::FrictionLaw

Description
    Rate-dependent friction mu(I) for Drucker-Prager and Mohr-Coulomb (WP 3.6a,
    decision 0037). Reference: prototype/constitutive/mu_i.py.

    Jop, Forterre and Pouliquen (2005) law

        mu_J(I) = mu_s + (mu_2 - mu_s)/(I0/I + 1),  I = gamma_dot d/sqrt(P/rho_s)

    with the optional Barker-Gray (2017) low-I branch (alpha > 0)

        mu(I) = sqrt(alpha/ln(A/I)) for I <= I1N,  A = I1N exp(alpha/mu_J(I1N)^2)

    where I1N is the first zero of the well-posedness function C(I) (Barker
    and Gray 2017, eq. 3.9). For mu_s = tan 35 deg the Jop law has no
    well-posed range, so the branch cannot be built and T0 runs alpha = 0.

    muIReturn wraps a rate-independent return. mu sets the mobilised angle
    phi_m = atan(mu); the base return at phi_m gives the plastic shear
    increment dgamma = sqrt(2 e:e) of the deviatoric plastic strain and the
    pressure P after the return, hence I = dgamma/dt d sqrt(rho_s/P). The root
    of R(mu) = mu - mu(I(mu)), which increases with mu, is found by the
    Illinois method to 1e-14 mu_2.

    Dictionary parameters (with the base model's, which give mu_s = tan(phi)):
        mu2 [-], I0 [-], d [m], rhoS [kg/m^3], optional muIAlpha [-] in (0, 2)
    mu2 absent or equal to tan(phi): the rate-independent model, unchanged.

SourceFiles
    celMuI.C
\*---------------------------------------------------------------------------*/

#ifndef celMuI_H
#define celMuI_H

#include "celConstitutiveKernel.H"

#include <functional>
#include <memory>

namespace cel
{
namespace constitutive
{

class FrictionLaw
{
    double muS_;
    double mu2_;
    double I0_;
    double alpha_;
    double d_;
    double rhoS_;
    double I1N_;
    double A_;

public:

    FrictionLaw
    (
        double muS,
        double mu2,
        double I0,
        double alpha,
        double d,
        double rhoS
    );

    //- The law for a base model with friction angle phiDeg, or null if p
    //  has no mu2 or mu2 equals tan(phi). Throws if mu2 < tan(phi) or a
    //  needed parameter is missing.
    static std::unique_ptr<FrictionLaw> New(const Parameters& p, double phiDeg);

    double muS() const { return muS_; }
    double mu2() const { return mu2_; }
    double alpha() const { return alpha_; }
    double I1N() const { return I1N_; }

    double muJop(double I) const;

    //- C(I) of Barker and Gray (2017) eq. 3.9 for the Jop law (> 0: ill posed)
    double wellPosedness(double I) const;

    //- mu(I); I = +infinity gives mu_2
    double mu(double I) const;

    //- I from a plastic shear increment, the stress after the return and dt
    double inertialNumber(double dgam, const Voigt sigNew, double dt) const;
};


//- sqrt(2 e:e) of the deviatoric plastic strain C (trial - sigNew), tensor
//  components (plastic_shear_increment in mu_i.py)
double plasticShearIncrement
(
    const Elastic& el,
    const Voigt trial,
    const Voigt sigNew
);


//- Rate-independent return at friction mu = tan(phi_m). Must give the
//  model's own update when mu == law.muS().
typedef std::function<void(double mu, Voigt sigNew, History& hNew)>
    FrictionReturn;


//- Rate-dependent return (rate_dependent in mu_i.py)
void muIReturn
(
    const FrictionLaw& law,
    const Elastic& el,
    const Voigt trial,
    double dt,
    const FrictionReturn& ret,
    Voigt sigNew,
    History& hNew
);


} // End namespace constitutive
} // End namespace cel

#endif
