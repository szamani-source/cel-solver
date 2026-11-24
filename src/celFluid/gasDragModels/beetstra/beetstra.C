/*---------------------------------------------------------------------------*\
    beetstra: gas-grain drag law (WP F.7, decision 0132)
\*---------------------------------------------------------------------------*/

#include "beetstra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{
    defineTypeNameAndDebug(beetstra, 0);
    addToRunTimeSelectionTable(gasDragModel, beetstra, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasDragModels::beetstra::beetstra(const dictionary& dict)
:
    gasDragModel(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModels::beetstra::K
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    const scalar phi = 1 - e;
    const scalar Re = e*rhoG*u*d/mu;
    scalar F = 10*phi/(e*e) + e*e*(1 + 1.5*Foam::sqrt(phi));
    if (Re > 0)
    {
        F +=
            0.413*Re/(24*e*e)
           *(1/e + 3*e*phi + 8.4*Foam::pow(Re, -0.343))
           /(1 + Foam::pow(10.0, 3*phi)*Foam::pow(Re, -0.5*(1 + 4*phi)));
    }
    return 18*mu*e*phi*F/(d*d);
}


// ************************************************************************* //
