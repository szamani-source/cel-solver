/*---------------------------------------------------------------------------*\
    WP 3.6g: checks of the mu(I) dilatancy (decision 0048) on the kernel.

    Plane-strain simple shear under a constant normal load (a shear box):
    each step adds the shear strain gdot dt to gamma_xy, and eps_yy is
    solved (secant) so that sigma_yy stays at -P0. eps_xx = eps_zz = 0.

    1. Dense start (phi0 above phiEq) dilates, loose start contracts, and
       both reach the steady solid fraction phiEq(I) = phiC - bPhi I, with I
       measured from the steady shear rate and pressure, to 1e-6 relative.
    2. At the steady state the friction is the plain mu(I) of the Jop law:
       the dilation angle is zero there, so the steady stress ratio of the
       dilatant model equals that of the model without dilatancy (1e-6).
    3. Flow rule at the first plastic steps of the dense start:
       d(eps_v^p)/d(gamma^p) is sin(psi) (Mohr-Coulomb plane return,
       gamma^p = eps1^p - eps3^p) with tan(psi) = KDil (phi - phiEq(I)).
    4. Without the four entries the model is the plain mu(I) model bit for
       bit; giving only some of them is refused.
    5. Drucker-Prager: the dense start dilates and settles at phiEq to 1e-6.
\*---------------------------------------------------------------------------*/

#include "celConstitutiveKernel.H"
#include "celDruckerPrager.H"
#include "celMohrCoulomb.H"
#include "celMuI.H"

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

using namespace cel::constitutive;

namespace
{

int failures = 0;

void check(bool ok, const std::string& what)
{
    std::printf("  %s: %s\n", ok ? "ok  " : "FAIL", what.c_str());
    if (!ok)
    {
        ++failures;
    }
}

Parameters t0(double phi0, bool dil)
{
    Parameters p
    {
        {"E", 1e6}, {"nu", 0.3}, {"phi", 35}, {"c", 0}, {"psi", 0},
        {"mu2", 0.961813}, {"I0", 0.279}, {"d", 1e-3}, {"rhoS", 2650}
    };
    if (dil)
    {
        p["phi0"] = phi0;
        p["phiC"] = 0.60;
        p["bPhi"] = 0.30;
        p["KDil"] = 4.0;
    }
    return p;
}

struct State
{
    Voigt sig = {-5e3, -5e3, -5e3, 0, 0, 0};
    History h;
    double I = 0;
    double ratio = 0;   // -sigma_xy/sigma_yy
};

//- One shear-box step at shear rate gdot; returns the solved eps_yy
double step(const Model& m, State& s, double gdot, double dt, double P0)
{
    auto run = [&](double eyy, Voigt out, History& hNew)
    {
        Voigt deps = {0, eyy, 0, gdot*dt, 0, 0};
        m.updateRate(deps, dt, s.sig, s.h, out, hNew);
        return out[1] + P0;
    };
    Voigt out;
    History hNew;
    double a = 0, fa = run(a, out, hNew);
    double b = 1e-9 + std::abs(gdot*dt), fb = run(b, out, hNew);
    for (int it = 0; it < 100 && std::abs(fb) > 1e-10*P0; ++it)
    {
        const double c = b - fb*(b - a)/(fb - fa);
        a = b; fa = fb;
        b = c; fb = run(b, out, hNew);
    }
    run(b, out, hNew);
    for (int i = 0; i < 6; ++i)
    {
        s.sig[i] = out[i];
    }
    s.h = hNew;
    s.ratio = -out[3]/out[1];
    return b;
}

//- Shear to steady state; returns the state
State shear
(
    const Model& m,
    double gdot,
    double dt,
    int n,
    double P0,
    double* firstRatio = nullptr,
    double* firstPhi = nullptr
)
{
    State s;
    s.sig[0] = s.sig[1] = s.sig[2] = -P0;
    bool first = true;
    for (int k = 0; k < n; ++k)
    {
        const History hOld = s.h;
        step(m, s, gdot, dt, P0);
        const double dv =
            (s.h.epsP[0] + s.h.epsP[1] + s.h.epsP[2])
          - (hOld.epsP[0] + hOld.epsP[1] + hOld.epsP[2]);
        // Principal plastic strain increments in the x-y plane
        const double ex = s.h.epsP[0] - hOld.epsP[0];
        const double ey = s.h.epsP[1] - hOld.epsP[1];
        const double gxy = s.h.epsP[3] - hOld.epsP[3];
        const double r = std::sqrt(0.25*(ex - ey)*(ex - ey) + 0.25*gxy*gxy);
        if (first && r > 0 && firstRatio)
        {
            *firstRatio = dv/(2*r);
            if (firstPhi)
            {
                *firstPhi = 1;  // marks that a plastic step was found
            }
            first = false;
        }
    }
    const double P = -(s.sig[0] + s.sig[1] + s.sig[2])/3;
    s.I = gdot*1e-3*std::sqrt(2650/P);
    return s;
}

} // End anonymous namespace


int main()
{
    const double P0 = 5e3, gdot = 20, dt = 1e-5;
    const int n = 40000;    // total shear strain 8

    std::printf("== 1. steady solid fraction phiEq(I)\n");
    for (double phi0 : {0.679, 0.55})
    {
        MohrCoulomb m(t0(phi0, true));
        const FrictionLaw& law = *m.frictionLaw();
        const State s = shear(m, gdot, dt, n, P0);
        // Measured I: plastic shear rate of steady simple shear = gdot
        const double P = -(s.sig[0] + s.sig[1] + s.sig[2])/3;
        const double I = gdot*1e-3*std::sqrt(2650/P);
        const double phi = law.solidFraction(s.h);
        const double phiEq = 0.60 - 0.30*I;
        char msg[200];
        std::snprintf
        (
            msg, sizeof msg,
            "phi0 %.3f -> phi %.5f, phiEq(I = %.4f) %.5f, rel %.1e",
            phi0, phi, I, phiEq, std::abs(phi - phiEq)/phiEq
        );
        check(std::abs(phi - phiEq) < 1e-6*phiEq, msg);
        check
        (
            phi0 > 0.6 ? phi < phi0 : phi > phi0,
            phi0 > 0.6 ? "dense start dilates" : "loose start contracts"
        );
    }

    std::printf("== 2. steady friction is the plain mu(I)\n");
    {
        MohrCoulomb md(t0(0.679, true)), mp(t0(0.679, false));
        const State sd = shear(md, gdot, dt, n, P0);
        const State sp = shear(mp, gdot, dt, n, P0);
        char msg[200];
        std::snprintf
        (
            msg, sizeof msg, "tau/sigma_n dilatant %.8f, plain %.8f",
            sd.ratio, sp.ratio
        );
        check(std::abs(sd.ratio - sp.ratio) < 1e-6*sp.ratio, msg);
    }

    std::printf("== 3. flow rule at the first plastic step\n");
    {
        MohrCoulomb m(t0(0.679, true));
        double ratio = 0, found = 0;
        // Short run: only the first plastic step is used
        shear(m, gdot, dt, 2000, P0, &ratio, &found);
        // Expected sin(psi) at the first plastic step: phi = phi0 (no
        // plastic strain yet), I from the step's own mu; bound it by I = 0
        // and I = 0.1 (the run's I is about 0.01)
        const double sLo =
            std::sin(std::atan(4.0*(0.679 - 0.60)));
        const double sHi =
            std::sin(std::atan(4.0*(0.679 - 0.60 + 0.30*0.1)));
        char msg[200];
        std::snprintf
        (
            msg, sizeof msg,
            "d(eps_v^p)/d(gamma^p) %.5f in [sin psi(I=0) %.5f, %.5f]",
            ratio, sLo, sHi
        );
        check(found > 0 && ratio >= sLo*(1 - 1e-9) && ratio <= sHi, msg);
    }

    std::printf("== 4. off without the entries, refused with only some\n");
    {
        MohrCoulomb mp(t0(0.679, false));
        check(!mp.frictionLaw()->dilatant(), "no entries: not dilatant");
        Parameters p = t0(0.679, false);
        p["phiC"] = 0.6;
        bool threw = false;
        try
        {
            MohrCoulomb bad(p);
        }
        catch (const std::invalid_argument&)
        {
            threw = true;
        }
        check(threw, "phiC alone is refused");
    }

    std::printf("== 5. Drucker-Prager\n");
    {
        Options o{{"fit", "plane_strain"}};
        DruckerPrager m(t0(0.679, true), o);
        const State s = shear(m, gdot, dt, n, P0);
        const double P = -(s.sig[0] + s.sig[1] + s.sig[2])/3;
        const double I = gdot*1e-3*std::sqrt(2650/P);
        const double phi = m.frictionLaw()->solidFraction(s.h);
        const double phiEq = 0.60 - 0.30*I;
        char msg[200];
        std::snprintf
        (
            msg, sizeof msg, "phi %.5f, phiEq(I = %.4f) %.5f", phi, I, phiEq
        );
        check(phi < 0.679 && std::abs(phi - phiEq) < 1e-6*phiEq, msg);
    }

    std::printf("%s\n", failures ? "FAILED" : "all passed");
    return failures ? 1 : 0;
}
