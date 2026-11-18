/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::MohrCoulomb

Description
    Perfectly plastic Mohr-Coulomb with non-associated flow (WP 2.5b).
    Parameters: E [Pa], nu [-], phi [deg] friction angle, c [Pa] cohesion,
    psi [deg] dilation angle, and optional sigmaT [Pa] tensile strength
    (Rankine tension cut-off, WP 2.5c; 0 <= sigmaT <= c cot(phi)).

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

    With sigmaT (WP 2.5c, decision 0015) three Rankine planes s_i <= sigmaT
    with associated flow are added. The Mohr-Coulomb return above stands if
    it lands inside them, so paths that never reach the cut-off are
    unchanged. Otherwise every active set of one to three planes out of
    {main, s1 = s2 edge, s2 = s3 edge, T1, T2, T3} is tried in a fixed order
    and the first admissible one is taken (multipliers >= 0, order kept,
    all yield functions <= 0), as _mc_tension_return in models.py.
    phi = 0 gives Tresca (no apex; a failed corner return throws).

    Optional mu(I) friction (WP 3.6a, celMuI.H): parameters mu2, I0, d, rhoS,
    muIAlpha. The return then uses phi_m = atan(mu(I)) for phi; sigmaT is
    lowered to c cot(phi_m) where it would lie above that apex. With the
    mu(I) dilatancy entries (phi0, phiC, bPhi, KDil; WP 3.6g) psi follows
    the state (FrictionLaw::tanPsi) in every return, also at mu_s.
\*---------------------------------------------------------------------------*/

#ifndef celMohrCoulomb_H
#define celMohrCoulomb_H

#include "celConstitutiveKernel.H"
#include "celMuI.H"

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
        apex,
        tensionPlane,       // s1 = sigmaT
        tensionEdge,        // s1 = s2 = sigmaT
        tensionApex,        // s1 = s2 = s3 = sigmaT
        tensionCorner       // Mohr-Coulomb and tension planes both active
    };

    //- Number of Return values (for counting)
    static constexpr int nReturns = 9;

private:

    Elastic el_;
    double c_;
    double sinPhi_;
    double cosPhi_;
    double sinPsi_;
    bool tension_;
    double sigmaT_;

    //- mu(I) friction, or null (rate-independent)
    std::unique_ptr<FrictionLaw> law_;

    //- Friction-dependent constants of one return
    struct Friction
    {
        double sinPhi;
        double cosPhi;
        double sigmaT;
        double sinPsi;
    };

    //- Return mapping with the given friction constants
    Return returnAt
    (
        const Friction& fr,
        const Voigt deps,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const;

    //- Active-set return onto Mohr-Coulomb plus Rankine planes, sorted
    //  principal trial w -> s. Returns the kind of return taken.
    Return tensionReturn
    (
        const Friction& fr,
        const double w[3],
        double scale,
        double s[3]
    ) const;

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

    void updateRate
    (
        const Voigt deps,
        double dt,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const override;

    bool rateDependent() const override
    {
        return bool(law_);
    }

    //- The mu(I) law, or null
    const FrictionLaw* frictionLaw() const override
    {
        return law_.get();
    }

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

    double shearModulus() const override
    {
        return el_.G;
    }

    //- Yield function at sigma (mc_yield in models.py), for tests
    double yield(const Voigt sigma) const;

    //- True if a tension cut-off sigmaT was given
    bool hasTensionCutOff() const
    {
        return tension_;
    }

    //- Tension cut-off yield function max(s_i) - sigmaT (rankine_yield in
    //  models.py); -infinity without a cut-off
    double tensionYield(const Voigt sigma) const;
};

} // End namespace constitutive
} // End namespace cel

#endif
