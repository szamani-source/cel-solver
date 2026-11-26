/*---------------------------------------------------------------------------*\
    ergun: gas-grain drag law (WP F.7, decision 0132)
\*---------------------------------------------------------------------------*/

#include "ergun.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{
    defineTypeNameAndDebug(ergun, 0);
    addToRunTimeSelectionTable(gasDragModel, ergun, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasDragModels::ergun::ergun(const dictionary& dict)
:
    gasDragModel(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModels::ergun::K
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    return
        150*(1 - e)*(1 - e)*mu/(e*d*d)
      + 1.75*(1 - e)*rhoG*u/d;
}


// ************************************************************************* //
