/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::MohrCoulomb

Description
    Perfectly plastic Mohr-Coulomb with non-associated flow (WP 2.5b).
    Parameters: E [Pa], nu [-], phi [deg] friction angle, c [Pa] cohesion,
    psi [deg] dilation angle.

    Yield function in sorted principal stresses s1 >= s2 >= s3 (tension
    positive):
        f = (s1 - s3) + (s1 + s3) sin(phi) - 2 c cos(phi)
    and the plastic potential with psi in place of phi.

    Return mapping in principal space, as mohr_coulomb in
    prototype/constitutive/models.py (decision 0004): main plane, then the
    edge the main return crosses (triaxial compression corner s1 = s2 or
    extension corner s2 = s3), then the apex at c cot(phi). Each return is a
    1x1 or 2x2 linear solve with the active planes, exact for perfect
    plasticity. The principal directions come from symmetricEigen (celEigen.H).

    No tension cut-off: the Python reference has none (decision 0004).
    phi = 0 gives Tresca (no apex; a failed corner return throws).
\*---------------------------------------------------------------------------*/

#ifndef celMohrCoulomb_H
#define celMohrCoulomb_H

#include "celConstitutiveKernel.H"

namespace cel
{
namespace constitutive
{

class MohrCoulomb
:
    public Model
{
public:

    //- Which return the last plastic step took (for tests)
    enum class Return
    {
        elastic,
        plane,
        compressionEdge,    // s1 = s2
        extensionEdge,      // s2 = s3
        apex
    };

private:

    Elastic el_;
    double c_;
    double sinPhi_;
    double cosPhi_;
    double sinPsi_;

public:

    explicit MohrCoulomb(const Parameters& p);

    const char* name() const override
    {
        return "mohrCoulomb";
    }

    void update
    (
        const Voigt deps,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const override;

    //- As update, and report which return was taken
    Return updateReturn
    (
        const Voigt deps,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const;

    double waveModulus() const override
    {
        return el_.waveModulus();
    }

    //- Yield function at sigma (mc_yield in models.py), for tests
    double yield(const Voigt sigma) const;
};

} // End namespace constitutive
} // End namespace cel

#endif
