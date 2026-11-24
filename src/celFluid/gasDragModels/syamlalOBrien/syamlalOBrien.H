/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModels::syamlalOBrien

Description
    Syamlal and O'Brien (1988), as in the MFiX theory guide (Syamlal,
    Rogers and O'Brien 1993):

        K = 0.75 e (1-e) rhoG u/(Vr^2 d) Cd(Re/Vr),
        Cd(Re/Vr) = (0.63 + 4.8 sqrt(Vr/Re))^2,  Re = rhoG u d/mu,
        Vr = (A - 0.06 Re + sqrt((0.06 Re)^2 + 0.12 Re (2B - A) + A^2))/2,
        A = e^4.14,  B = 0.8 e^1.28 (e <= 0.85), e^2.65 (e > 0.85)

    written as K = 0.75 e (1-e) mu/(Vr^2 d^2) (0.63 sqrt(Re) + 4.8 sqrt(Vr))^2,
    finite at u = 0.

    A user option (WP F.7, decision 0125 item 7); every test uses the
    default ergunWenYu (Sia, 2026-11-24). Checked against
    tools/gas_reference.py by tests/test_drag_laws.sh.

    \verbatim
        drag    syamlalOBrien;        // constant/gasProperties
    \endverbatim

SourceFiles
    syamlalOBrien.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModels_syamlalOBrien_H
#define Foam_gasDragModels_syamlalOBrien_H

#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{

class syamlalOBrien
:
    public gasDragModel
{
public:

    //- Runtime type information
    TypeName("syamlalOBrien");


    // Constructors

        explicit syamlalOBrien(const dictionary& dict);


    //- Destructor
    virtual ~syamlalOBrien() = default;


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
