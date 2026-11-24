/*---------------------------------------------------------------------------*\
    diFelice: gas-grain drag law (WP F.7, decision 0132)
\*---------------------------------------------------------------------------*/

#include "diFelice.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{
    defineTypeNameAndDebug(diFelice, 0);
    addToRunTimeSelectionTable(gasDragModel, diFelice, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasDragModels::diFelice::diFelice(const dictionary& dict)
:
    gasDragModel(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModels::diFelice::K
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    const scalar Re = e*rhoG*u*d/mu;
    const scalar s = 0.63*Foam::sqrt(Re) + 4.8;
    const scalar beta =
        Re > 0
      ? 3.7 - 0.65*Foam::exp(-0.5*Foam::sqr(1.5 - Foam::log10(Re)))
      : 3.7;
    return 0.75*s*s*mu*(1 - e)/(d*d)*Foam::pow(e, 1 - beta);
}


// ************************************************************************* //
