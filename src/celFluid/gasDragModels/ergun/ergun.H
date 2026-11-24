/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModels::ergun

Description
    Ergun (1952) at every gas fraction:

        K = 150 (1-e)^2 mu/(e d^2) + 1.75 (1-e) rhoG u/d

    (the dense branch of ergunWenYu; it tends to 0, not to one
    sphere's drag, as e -> 1, so it is meant for packed beds).

    A user option (WP F.7, decision 0125 item 7); every test uses the
    default ergunWenYu (Sia, 2026-11-24). Checked against
    tools/gas_reference.py by tests/test_drag_laws.sh.

    \verbatim
        drag    ergun;        // constant/gasProperties
    \endverbatim

SourceFiles
    ergun.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModels_ergun_H
#define Foam_gasDragModels_ergun_H

#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{

class ergun
:
    public gasDragModel
{
public:

    //- Runtime type information
    TypeName("ergun");


    // Constructors

        explicit ergun(const dictionary& dict);


    //- Destructor
    virtual ~ergun() = default;


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
