/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::DruckerPrager

Description
    Perfectly plastic Drucker-Prager with non-associated flow (WP 2.5a).

        f = sqrt(J2) + eta p - xi c,    g = sqrt(J2) + etaBar p

    p is the mean stress, positive in tension. eta and xi come from phi,
    etaBar from psi, by one of three matches to Mohr-Coulomb (option "fit"):
    compression (outer cone, default), extension (inner cone), plane_strain.
    Return to the smooth cone, or to the apex (zero deviator, p = xi c/eta)
    when the cone return would overshoot the axis.

    Parameters: E [Pa], nu [-], phi [deg], c [Pa], psi [deg].
    Options:    fit = compression | extension | plane_strain.
    Reference: dp_coefficients and drucker_prager in models.py.

SourceFiles
    celDruckerPrager.C
\*---------------------------------------------------------------------------*/

#ifndef celDruckerPrager_H
#define celDruckerPrager_H

#include "celConstitutiveKernel.H"

namespace cel
{
namespace constitutive
{

class DruckerPrager
:
    public Model
{
public:

    //- Match of the cone to Mohr-Coulomb
    enum class Fit
    {
        compression,
        extension,
        planeStrain
    };

    //- Fit from its dictionary word; throws std::invalid_argument if unknown
    static Fit fitFromName(const std::string& name);

    //- Dictionary word of a fit
    static const char* fitName(Fit fit);


private:

    Elastic el_;
    double c_;
    Fit fit_;
    double eta_;
    double xi_;
    double etaBar_;


public:

    DruckerPrager(const Parameters& p, const Options& opts);

    const char* name() const override
    {
        return "druckerPrager";
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

    Fit fit() const
    {
        return fit_;
    }

    double eta() const
    {
        return eta_;
    }

    double xi() const
    {
        return xi_;
    }

    double etaBar() const
    {
        return etaBar_;
    }

    //- Yield function sqrt(J2) + eta p - xi c, for tests
    double yield(const Voigt sigma) const;
};

} // End namespace constitutive
} // End namespace cel

#endif
