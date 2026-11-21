/*---------------------------------------------------------------------------*\
Class
    cel::constitutive::Cap

Description
    Compaction cap (WP K.1, decision 0113): a volumetric yield p <= p_c on
    top of a frictional model, with p_c hardening so that loading along
    p = p_c follows the normal compression line

        eps_v = lambda* ln(p/p0)        (eps_v compressive, natural strain)

    lambda* the modified compression index (lambda/(1 + e0), entry
    capLambda), measured in an oedometer or isotropic test. Unloading and
    reloading below p_c are elastic with the model's K (E, nu).

    Applied after the frictional update to its stress sigma' (P' = -tr/3):
    if P' > p_c, the total volumetric strain of the step is shared so that
    the new stress sits on the line: p_c,new = p_c exp((P' - p_c)/(K lambda*))
    (closed form of K lambda* ln(p_c,new/p_c) = P' - p_c), and the stress is
    scaled by p_c,new/P' (towards the apex, so it stays inside the cone of
    a cohesionless model). Needs K lambda* > P (the elastic line stiffer
    than the compression line); otherwise the step stays elastic.

    p_c starts at max(P of the old stress, capP0) (entry capP0, default
    100 Pa): material is normally consolidated where it is first loaded.

SourceFiles
    celCap.C
\*---------------------------------------------------------------------------*/

#ifndef celCap_H
#define celCap_H

#include "celConstitutiveKernel.H"

#include <memory>

namespace cel
{
namespace constitutive
{

class Cap
{
    double lambda_;
    double p0_;

public:

    Cap(double lambda, double p0);

    //- The cap of parameters p, or null if p has no capLambda (or 0)
    static std::unique_ptr<Cap> New(const Parameters& p);

    double lambda() const { return lambda_; }
    double p0() const { return p0_; }

    //- Cap the stress sig (in place) and update h (pc, plastic strain);
    //  sigOld gives the starting p_c where h.pc is not set yet
    void apply(const Elastic& el, const Voigt sigOld, Voigt sig, History& h)
        const;
};


} // End namespace constitutive
} // End namespace cel

#endif
