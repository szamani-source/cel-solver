/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModels::ergunWenYu

Description
    Ergun-Wen-Yu drag (Gidaspow 1994), the default law and the one every
    test uses (Sia, 2026-11-24; decision 0125 item 7). With e the gas
    fraction and u the slip speed:

        e <  0.8:  K = 150 (1-e)^2 mu/(e d^2) + 1.75 (1-e) rhoG u/d  (Ergun)
        e >= 0.8:  K = 0.75 Cd e (1-e) rhoG u/d e^-2.65              (Wen-Yu)
                   Cd = 24/Re (1 + 0.15 Re^0.687) below Re = 1000, 0.44
                   above, Re = e rhoG u d/mu

    as tools/gas_reference.py drag_ergun_wen_yu. The Wen-Yu branch is
    written with Cd Re, so u = 0 gives its Stokes limit
    K = 18 mu (1-e) e^-2.65/d^2. On the Ergun branch the slip at a given
    pressure gradient is the root of a quadratic (closed form).

    \verbatim
        drag    ergunWenYu;     // constant/gasProperties (the default)
    \endverbatim

SourceFiles
    ergunWenYu.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModels_ergunWenYu_H
#define Foam_gasDragModels_ergunWenYu_H

#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{

class ergunWenYu
:
    public gasDragModel
{
public:

    //- Runtime type information
    TypeName("ergunWenYu");

    //- Ergun below this gas fraction, Wen and Yu at and above it
    static constexpr scalar eSwitch = 0.8;


    // Constructors

        explicit ergunWenYu(const dictionary& dict);


    //- Destructor
    virtual ~ergunWenYu() = default;


    // Member Functions

        virtual scalar K
        (
            const scalar e,
            const scalar u,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;

        //- Ergun branch: K = K0 + K1 u (WP F.7)
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
            if (!(e < eSwitch))
            {
                return false;
            }
            K0 = 150*(1 - e)*(1 - e)*mu/(e*d*d);
            K1 = 1.75*(1 - e)*rhoG/d;
            return true;
        }

        //- Exact on the Ergun branch (K linear in u), the default above
        virtual scalar dKdu
        (
            const scalar e,
            const scalar u,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;

        //- Closed form on the Ergun branch, the default root finder on
        //  the Wen-Yu branch
        virtual scalar slip
        (
            const scalar e,
            const scalar G,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;
};


} // End namespace gasDragModels
} // End namespace Foam

#endif

// ************************************************************************* //
