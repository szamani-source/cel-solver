#include "celConstitutiveKernel.H"
#include "celLinearElastic.H"
#include "celVonMises.H"

#include <stdexcept>

namespace cel
{
namespace constitutive
{

Elastic::Elastic(double E, double nu)
:
    K(E/(3.0*(1.0 - 2.0*nu))),
    G(E/(2.0*(1.0 + nu)))
{}


void Elastic::addStress(const Voigt sigma, const Voigt deps, Voigt out) const
{
    const double lambda = K - 2.0*G/3.0;
    const double ev = deps[0] + deps[1] + deps[2];
    for (int i = 0; i < 3; ++i)
    {
        out[i] = sigma[i] + (lambda*ev + 2.0*G*deps[i]);
        out[i + 3] = sigma[i + 3] + G*deps[i + 3];
    }
}


void Elastic::strain(const Voigt dsigma, Voigt out) const
{
    // Normal part: tr/(9K) + dev/(2G); shear (engineering): tau/G
    const double p = (dsigma[0] + dsigma[1] + dsigma[2])/3.0;
    for (int i = 0; i < 3; ++i)
    {
        out[i] = p/(3.0*K) + (dsigma[i] - p)/(2.0*G);
        out[i + 3] = dsigma[i + 3]/G;
    }
}


double meanStress(const Voigt s)
{
    return (s[0] + s[1] + s[2])/3.0;
}


void deviator(const Voigt s, Voigt out)
{
    const double p = meanStress(s);
    for (int i = 0; i < 3; ++i)
    {
        out[i] = s[i] - p;
        out[i + 3] = s[i + 3];
    }
}


double J2(const Voigt s)
{
    Voigt d;
    deviator(s, d);
    return
        0.5*(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
      + (d[3]*d[3] + d[4]*d[4] + d[5]*d[5]);
}


double parameter(const Parameters& p, const std::string& key)
{
    const auto it = p.find(key);
    if (it == p.end())
    {
        throw std::invalid_argument("missing material parameter '" + key + "'");
    }
    return it->second;
}


double parameter(const Parameters& p, const std::string& key, double def)
{
    const auto it = p.find(key);
    return it == p.end() ? def : it->second;
}


std::unique_ptr<Model> makeModel(const std::string& name, const Parameters& p)
{
    if (name == "linearElastic")
    {
        return std::unique_ptr<Model>(new LinearElastic(p));
    }
    if (name == "vonMises")
    {
        return std::unique_ptr<Model>(new VonMises(p));
    }
    throw std::invalid_argument
    (
        "unknown constitutive model '" + name
      + "' (valid: linearElastic vonMises)"
    );
}

} // End namespace constitutive
} // End namespace cel
