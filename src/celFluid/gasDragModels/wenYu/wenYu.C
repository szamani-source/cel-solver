/*---------------------------------------------------------------------------*\
    wenYu: gas-grain drag law (WP F.7, decision 0132)
\*---------------------------------------------------------------------------*/

#include "wenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{
    defineTypeNameAndDebug(wenYu, 0);
    addToRunTimeSelectionTable(gasDragModel, wenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasDragModels::wenYu::wenYu(const dictionary& dict)
:
    gasDragModel(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModels::wenYu::K
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    const scalar Re = e*rhoG*u*d/mu;
    const scalar CdRe =
        Re < 1000
      ? 24*(1 + 0.15*Foam::pow(Re, 0.687))
      : 0.44*Re;
    return 0.75*CdRe*mu*(1 - e)/(d*d)*Foam::pow(e, -2.65);
}


// ************************************************************************* //
