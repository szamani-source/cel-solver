/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModels::wenYu

Description
    Wen and Yu (1966) at every gas fraction:

        K = 0.75 Cd e (1-e) rhoG u/d e^-2.65,
        Cd = 24/Re (1 + 0.15 Re^0.687) below Re = 1000, 0.44 above,
        Re = e rhoG u d/mu

    (the dilute branch of ergunWenYu), written with Cd Re so u = 0 gives
    the Stokes limit 18 mu (1-e) e^-2.65/d^2.

    A user option (WP F.7, decision 0125 item 7); every test uses the
    default ergunWenYu (Sia, 2026-11-24). Checked against
    tools/gas_reference.py by tests/test_drag_laws.sh.

    \verbatim
        drag    wenYu;        // constant/gasProperties
    \endverbatim

SourceFiles
    wenYu.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModels_wenYu_H
#define Foam_gasDragModels_wenYu_H

#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{

class wenYu
:
    public gasDragModel
{
public:

    //- Runtime type information
    TypeName("wenYu");


    // Constructors

        explicit wenYu(const dictionary& dict);


    //- Destructor
    virtual ~wenYu() = default;


    // Member Functions

        virtual scalar K
        (
            const scalar e,
            const scalar u,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;
};


} // End namespace gasDragModels
} // End namespace Foam

#endif

// ************************************************************************* //
