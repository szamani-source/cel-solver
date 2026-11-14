/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::VonMises

Description
    von Mises (J2) plasticity with linear isotropic hardening, radial return
    (exact for linear hardening). Parameters: E [Pa], nu [-], sigma_y [Pa],
    H [Pa] (hardening modulus, default 0 = perfectly plastic).
    Reference: von_mises in models.py.
\*---------------------------------------------------------------------------*/

#ifndef celVonMises_H
#define celVonMises_H

#include "celConstitutiveKernel.H"

namespace cel
{
namespace constitutive
{

class VonMises
:
    public Model
{
    Elastic el_;
    double sigmaY_;
    double H_;

public:

    explicit VonMises(const Parameters& p);

    const char* name() const override
    {
        return "vonMises";
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

    //- Yield function q - (sigma_y + H epsPEq), for tests
    double yield(const Voigt sigma, const History& h) const;
};

} // End namespace constitutive
} // End namespace cel

#endif
