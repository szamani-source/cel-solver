#include "celLinearElastic.H"

namespace cel
{
namespace constitutive
{

LinearElastic::LinearElastic(const Parameters& p)
:
    el_(parameter(p, "E"), parameter(p, "nu"))
{}


void LinearElastic::update
(
    const Voigt deps,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    el_.addStress(sigOld, deps, sigNew);
    hNew = hOld;
    hNew.apex = false;
}

} // End namespace constitutive
} // End namespace cel
