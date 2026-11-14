/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::LinearElastic

Description
    Isotropic linear elasticity, sigma_new = sigma_old + D deps.
    Parameters: E [Pa], nu [-]. Reference: linear_elastic in models.py.
\*---------------------------------------------------------------------------*/

#ifndef celLinearElastic_H
#define celLinearElastic_H

#include "celConstitutiveKernel.H"

namespace cel
{
namespace constitutive
{

class LinearElastic
:
    public Model
{
    Elastic el_;

public:

    explicit LinearElastic(const Parameters& p);

    const char* name() const override
    {
        return "linearElastic";
    }

    void update
    (
        const Voigt deps,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const override;

    double waveModulus() const override
    {
        return el_.waveModulus();
    }

    double shearModulus() const override
    {
        return el_.G;
    }
};

} // End namespace constitutive
} // End namespace cel

#endif
