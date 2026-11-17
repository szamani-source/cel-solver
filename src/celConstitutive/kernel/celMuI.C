#include "celMuI.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cel
{
namespace constitutive
{

namespace
{
constexpr double pi = 3.14159265358979323846;
}


FrictionLaw::FrictionLaw
(
    double muS,
    double mu2,
    double I0,
    double alpha,
    double d,
    double rhoS
)
:
    muS_(muS),
    mu2_(mu2),
    I0_(I0),
    alpha_(alpha),
    d_(d),
    rhoS_(rhoS),
    I1N_(0),
    A_(0)
{
    if (!(mu2_ > muS_ && muS_ > 0 && I0_ > 0))
    {
        throw std::invalid_argument("mu(I): need mu2 > tan(phi) > 0 and I0 > 0");
    }
    if (!(d_ > 0 && rhoS_ > 0))
    {
        throw std::invalid_argument("mu(I): need d > 0 and rhoS > 0");
    }
    if (alpha_ < 0 || alpha_ >= 2)
    {
        throw std::invalid_argument
        (
            "mu(I): muIAlpha must be in (0, 2) (0 or absent: no regularisation)"
        );
    }
    if (alpha_ == 0)
    {
        return;
    }

    // First neutral point: scan log10 I in [-12, 3] by 0.01, then bisect
    const int n = 1500;
    double prev = std::pow(10.0, -12.0);
    if (wellPosedness(prev) <= 0)
    {
        throw std::invalid_argument
        (
            "mu(I): already well posed at I = 1e-12, no regularisation needed"
        );
    }
    for (int k = 1; k <= n; ++k)
    {
        const double I = std::pow(10.0, -12.0 + 15.0*k/n);
        if (wellPosedness(I) <= 0)
        {
            double lo = std::log(prev), hi = std::log(I);
            for (int it = 0; it < 200; ++it)
            {
                const double mid = 0.5*(lo + hi);
                if (wellPosedness(std::exp(mid)) > 0)
                {
                    lo = mid;
                }
                else
                {
                    hi = mid;
                }
                if (hi - lo <= 1e-15)
                {
                    break;
                }
            }
            I1N_ = std::exp(hi);
            const double m = muJop(I1N_);
            A_ = I1N_*std::exp(alpha_/(m*m));
            return;
        }
        prev = I;
    }
    throw std::invalid_argument
    (
        "mu(I): these parameters give no well-posed range of I, so the "
        "Barker-Gray regularisation cannot be built; leave muIAlpha out"
    );
}


std::unique_ptr<FrictionLaw> FrictionLaw::New
(
    const Parameters& p,
    double phiDeg
)
{
    if (!p.count("mu2"))
    {
        return nullptr;
    }
    const double muS = std::tan(phiDeg*pi/180.0);
    const double mu2 = parameter(p, "mu2");
    if (mu2 == muS)
    {
        return nullptr;
    }
    if (mu2 < muS)
    {
        throw std::invalid_argument("mu(I): mu2 must not be below tan(phi)");
    }
    return std::unique_ptr<FrictionLaw>
    (
        new FrictionLaw
        (
            muS,
            mu2,
            parameter(p, "I0"),
            parameter(p, "muIAlpha", 0),
            parameter(p, "d"),
            parameter(p, "rhoS")
        )
    );
}


double FrictionLaw::muJop(double I) const
{
    if (std::isinf(I))
    {
        return mu2_;
    }
    return muS_ + (mu2_ - muS_)*I/(I0_ + I);
}


double FrictionLaw::wellPosedness(double I) const
{
    const double mu = muJop(I);
    const double dmu = (mu2_ - muS_)*I0_/((I0_ + I)*(I0_ + I));
    const double nu = I*dmu/mu;
    const double t = 1.0 - 0.5*nu;
    return 4.0*nu*nu - 4.0*nu + mu*mu*t*t;
}


double FrictionLaw::mu(double I) const
{
    if (I <= 0)
    {
        return alpha_ > 0 ? 0.0 : muS_;
    }
    if (alpha_ > 0 && I <= I1N_)
    {
        return std::sqrt(alpha_/std::log(A_/I));
    }
    return muJop(I);
}


double FrictionLaw::inertialNumber
(
    double dgam,
    const Voigt sigNew,
    double dt
) const
{
    if (dgam == 0)
    {
        return 0;
    }
    const double P = -meanStress(sigNew);
    // P <= 0, or so small that rhoS/P would overflow (a floating-point
    // trap in OpenFOAM; seen at a free surface, WP 3.6s): I is infinite
    if (!(P > 1e-300*rhoS_))
    {
        return std::numeric_limits<double>::infinity();
    }
    return dgam/dt*d_*std::sqrt(rhoS_/P);
}


double plasticShearIncrement
(
    const Elastic& el,
    const Voigt trial,
    const Voigt sigNew
)
{
    Voigt dsig, depsP;
    for (int i = 0; i < 6; ++i)
    {
        dsig[i] = trial[i] - sigNew[i];
    }
    el.strain(dsig, depsP);
    double e[6];
    for (int i = 0; i < 6; ++i)
    {
        e[i] = i < 3 ? depsP[i] : 0.5*depsP[i];
    }
    const double ev = (e[0] + e[1] + e[2])/3.0;
    double n2 = 0;
    for (int i = 0; i < 3; ++i)
    {
        n2 += (e[i] - ev)*(e[i] - ev) + 2.0*e[i + 3]*e[i + 3];
    }
    return std::sqrt(2.0*n2);
}


namespace
{

struct Eval
{
    double r;
    Voigt sig;
    History h;
};

} // End anonymous namespace


void muIReturn
(
    const FrictionLaw& law,
    const Elastic& el,
    const Voigt trial,
    double dt,
    const FrictionReturn& ret,
    Voigt sigNew,
    History& hNew
)
{
    if (!(dt > 0))
    {
        throw std::invalid_argument("mu(I): the time step must be positive");
    }

    double lastI = 0;
    auto at = [&](double mu, Eval& e)
    {
        ret(mu, e.sig, e.h);
        lastI = law.inertialNumber
        (
            plasticShearIncrement(el, trial, e.sig), e.sig, dt
        );
        e.r = mu - law.mu(lastI);
    };
    auto take = [&](const Eval& e)
    {
        std::copy(e.sig, e.sig + 6, sigNew);
        hNew = e.h;
    };

    Eval ea, eb;
    at(law.muS(), ea);
    if (law.alpha() == 0 && lastI == 0)
    {
        take(ea);       // elastic below mu_s (Jop law)
        return;
    }
    if (ea.r == 0)
    {
        take(ea);
        return;
    }

    double a, b;
    if (ea.r < 0)
    {
        a = law.muS();
        b = law.mu2();
        at(b, eb);
    }
    else
    {
        a = 1e-3*law.muS();
        b = law.muS();
        eb = ea;
        at(a, ea);
        if (ea.r >= 0)
        {
            take(ea);   // creep below double precision
            return;
        }
    }
    if (eb.r <= 0)
    {
        take(eb);       // root at mu_2 (I infinite, e.g. P = 0)
        return;
    }

    // Illinois (modified regula falsi), as _illinois in mu_i.py
    const double tol = 1e-14*law.mu2();
    double fa = ea.r, fb = eb.r;
    int side = 0;
    Eval ec;
    for (int it = 0; it < 200; ++it)
    {
        double c = (a*fb - b*fa)/(fb - fa);
        if (!(a < c && c < b))
        {
            c = 0.5*(a + b);
        }
        at(c, ec);
        if (ec.r == 0)
        {
            take(ec);
            return;
        }
        if (ec.r < 0)
        {
            a = c;
            fa = ec.r;
            ea = ec;
            if (side == -1)
            {
                fb *= 0.5;
            }
            side = -1;
        }
        else
        {
            b = c;
            fb = ec.r;
            eb = ec;
            if (side == 1)
            {
                fa *= 0.5;
            }
            side = 1;
        }
        if (b - a <= tol)
        {
            take(std::abs(fa) < std::abs(fb) ? ea : eb);
            return;
        }
    }
    throw std::runtime_error("mu(I): return did not converge");
}

} // End namespace constitutive
} // End namespace cel
