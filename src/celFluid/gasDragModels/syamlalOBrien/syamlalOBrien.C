/*---------------------------------------------------------------------------*\
    syamlalOBrien: gas-grain drag law (WP F.7, decision 0132)
\*---------------------------------------------------------------------------*/

#include "syamlalOBrien.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{
    defineTypeNameAndDebug(syamlalOBrien, 0);
    addToRunTimeSelectionTable(gasDragModel, syamlalOBrien, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasDragModels::syamlalOBrien::syamlalOBrien(const dictionary& dict)
:
    gasDragModel(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModels::syamlalOBrien::K
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    const scalar Re = rhoG*u*d/mu;
    const scalar A = Foam::pow(e, 4.14);
    const scalar B = e <= 0.85 ? 0.8*Foam::pow(e, 1.28) : Foam::pow(e, 2.65);
    const scalar Vr =
        0.5
       *(
            A - 0.06*Re
          + Foam::sqrt(Foam::sqr(0.06*Re) + 0.12*Re*(2*B - A) + A*A)
        );
    const scalar s = 0.63*Foam::sqrt(Re) + 4.8*Foam::sqrt(Vr);
    return 0.75*e*(1 - e)*mu/(Vr*Vr*d*d)*s*s;
}


// ************************************************************************* //
