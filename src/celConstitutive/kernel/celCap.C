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
    const double pc = h.pc > 0 ? h.pc : std::max(-meanStress(sigOld), p0_);
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
