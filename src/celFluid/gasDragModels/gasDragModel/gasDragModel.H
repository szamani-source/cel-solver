/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModel

Description
    Run-time selected drag law between the pore gas and the grains (version
    2, WP F.2, decision 0125 item 7 and decision 0127). Each law gives the
    momentum exchange coefficient K [kg/(m^3 s)] in

        phi_g grad p = -K (u_g - u_s)

    for a packing with gas fraction e, slip speed u = |u_g - u_s|, grain
    diameter d, gas viscosity mu and gas density rhoG. In Stage A (pore
    pressure, F.2 to F.5) the superficial flux q = e (u_g - u_s) obeys
    q = -lambda grad p with the mobility lambda = e^2/K(e, u), where u
    follows from |grad p| = K(e, u) u / e (slip()).

    Selected in constant/gasProperties by `drag <name>;` (default
    ergunWenYu, the law every test uses, Sia 2026-11-24); optional
    coefficients in `<name>Coeffs`. F.2 has ergunWenYu; F.7 (decision 0132)
    adds ergun, wenYu, diFelice, syamlalOBrien and beetstra, each checked
    against tools/gas_reference.py by tests/test_drag_laws.sh.

SourceFiles
    gasDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModel_H
#define Foam_gasDragModel_H

#include "dictionary.H"
#include "autoPtr.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class gasDragModel
{
public:

    //- Runtime type information
    TypeName("gasDragModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            gasDragModel,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        explicit gasDragModel(const dictionary&) {}

        gasDragModel(const gasDragModel&) = delete;
        void operator=(const gasDragModel&) = delete;


    // Selectors

        //- The law called name, with coefficients from dict (may be empty)
        static autoPtr<gasDragModel> New
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~gasDragModel() = default;


    // Member Functions

        //- Momentum exchange coefficient K [kg/(m^3 s)], 0 < e < 1,
        //  slip u >= 0 (u = 0 gives the Darcy limit)
        virtual scalar K
        (
            const scalar e,
            const scalar u,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const = 0;

        //- dK/du [kg/(m^3 s) per m/s]; the default is a central
        //  difference (forward at u = 0)
        virtual scalar dKdu
        (
            const scalar e,
            const scalar u,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;

        //- Slip speed u >= 0 at which K(e, u) u / e = G (pressure gradient
        //  magnitude G [Pa/m]). The default brackets the root between 0
        //  and the Darcy speed G e / K(e, 0), doubled until it brackets
        //  (K may fall with u, WP F.7), and solves it by the Illinois
        //  method to 1e-14 relative.
        virtual scalar slip
        (
            const scalar e,
            const scalar G,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;

        //- True where K = K0 + K1 u exactly (an Ergun form), with K0 and
        //  K1 set; the default is false (WP F.7: closed-form slip of the
        //  pore gas momentum)
        virtual bool affine
        (
            const scalar e,
            const scalar d,
            const scalar mu,
            const scalar rhoG,
            scalar& K0,
            scalar& K1
        ) const
        {
            return false;
        }

        //- Zero-slip permeability mu e^2 / K(e, 0) [m^2]; 0 for e <= 0,
        //  GREAT for e >= 1 (no grains)
        scalar permeability
        (
            const scalar e,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;

        //- Mobility lambda = e^2 / K(e, slip(G)) [m^2/(Pa s)] at pressure
        //  gradient G: q = -lambda grad p. 0 for e <= 0, GREAT/mu for
        //  e >= 1.
        scalar mobility
        (
            const scalar e,
            const scalar G,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;

        //- Tangent mobility d|q|/d|grad p| = e^2/(K + u dK/du) at the slip
        //  of G [m^2/(Pa s)], the Newton matrix of the pore pressure
        //  (decision 0127); equals mobility() where K does not depend on u
        scalar tangentMobility
        (
            const scalar e,
            const scalar G,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;
};


/*---------------------------------------------------------------------------*\
                      Class inertialGasDrag Declaration
\*---------------------------------------------------------------------------*/

//- A drag law with the inertia of the pore gas over a step dt added (WP
//  F.7, decision 0132): K + e rhoG/dt. The gas momentum in the pores,
//  e rhoG (w - w_o)/dt = -e grad p - K w (w the slip), is then
//  K~ w = e rhoG w_o/dt - e grad p, a drag law in K~ against the gradient
//  less the old slip's part. Not run-time selectable; poreGas makes it.
class inertialGasDrag
:
    public gasDragModel
{
    const gasDragModel& base_;

    //- 1/dt [1/s]
    scalar rDt_;

public:

    TypeName("inertialGasDrag");

    inertialGasDrag(const gasDragModel& base, const scalar dt)
    :
        gasDragModel(dictionary::null),
        base_(base),
        rDt_(1/dt)
    {}

    virtual ~inertialGasDrag() = default;

    virtual scalar K
    (
        const scalar e,
        const scalar u,
        const scalar d,
        const scalar mu,
        const scalar rhoG
    ) const
    {
        return base_.K(e, u, d, mu, rhoG) + e*rhoG*rDt_;
    }

    virtual scalar dKdu
    (
        const scalar e,
        const scalar u,
        const scalar d,
        const scalar mu,
        const scalar rhoG
    ) const
    {
        return base_.dKdu(e, u, d, mu, rhoG);
    }

    virtual bool affine
    (
        const scalar e,
        const scalar d,
        const scalar mu,
        const scalar rhoG,
        scalar& K0,
        scalar& K1
    ) const
    {
        if (base_.affine(e, d, mu, rhoG, K0, K1))
        {
            K0 += e*rhoG*rDt_;
            return true;
        }
        return false;
    }

    //- Closed form where the law is affine in u (K u/e = G is then a
    //  quadratic), else the default root finder
    virtual scalar slip
    (
        const scalar e,
        const scalar G,
        const scalar d,
        const scalar mu,
        const scalar rhoG
    ) const
    {
        scalar K0, K1;
        if (G > 0 && e > 0 && e < 1 && affine(e, d, mu, rhoG, K0, K1))
        {
            const scalar a = K0/e, b = K1/e;
            return 2*G/(a + Foam::sqrt(a*a + 4*b*G));
        }
        return gasDragModel::slip(e, G, d, mu, rhoG);
    }
};


} // End namespace Foam

#endif

// ************************************************************************* //
