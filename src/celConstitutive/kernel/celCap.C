#include "celCap.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cel
{
namespace constitutive
{

Cap::Cap(double lambda, double p0)
:
    lambda_(lambda),
    p0_(p0)
{}


std::unique_ptr<Cap> Cap::New(const Parameters& p)
{
    const double lambda = parameter(p, "capLambda", 0);
    if (lambda <= 0)
    {
        return nullptr;
    }
    const double p0 = parameter(p, "capP0", 100);
    if (!(p0 > 0))
    {
        throw std::invalid_argument("cap: capP0 must be positive");
    }
    return std::unique_ptr<Cap>(new Cap(lambda, p0));
}


void Cap::apply
(
    const Elastic& el,
    const Voigt sigOld,
    Voigt sig,
    History& h
) const
{
    // p_c is at least capP0 and at least the old stress's P (the old stress
    // lies inside the cap). A p_c of 0 means "not set": fed material, or a
    // cell not loaded yet. The remap mixes such cells into loaded ones by
    // mass (and with a limiter), so p_c can arrive anywhere between 0 and
    // the loaded value, and which one can hinge on round-off. Taking the
    // largest of the three is continuous in p_c and gives a cell that was
    // never loaded p_c = max(P, capP0), as before (WP F.8, decision 0133).
    const double pc =
        std::max(h.pc, std::max(-meanStress(sigOld), p0_));
    h.pc = pc;
    const double P = -meanStress(sig);
    if (!(P > pc))
    {
        return;
    }
    const double pcNew = pc*std::exp(std::min((P - pc)/(el.K*lambda_), 50.0));
    if (!(pcNew < P))
    {
        // Elastic line softer than the compression line: no compaction
        return;
    }
    Voigt before;
    std::copy(sig, sig + 6, before);
    const double f = pcNew/P;
    for (int i = 0; i < 6; ++i)
    {
        sig[i] *= f;
    }
    const History hBefore = h;
    plasticCorrector(el, before, sig, hBefore, h.apex, h);
    h.pc = pcNew;
}


} // End namespace constitutive
} // End namespace cel
