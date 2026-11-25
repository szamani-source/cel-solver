/*---------------------------------------------------------------------------*\
    ergunWenYu: Gidaspow's Ergun / Wen-Yu drag (WP F.2, decision 0127)
\*---------------------------------------------------------------------------*/

#include "ergunWenYu.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{
    defineTypeNameAndDebug(ergunWenYu, 0);
    addToRunTimeSelectionTable(gasDragModel, ergunWenYu, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gasDragModels::ergunWenYu::ergunWenYu(const dictionary& dict)
:
    gasDragModel(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModels::ergunWenYu::K
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    if (e < eSwitch)
    {
        return
            150*(1 - e)*(1 - e)*mu/(e*d*d)
          + 1.75*(1 - e)*rhoG*u/d;
    }

    // Cd Re of Schiller and Naumann (0.44 Re above Re = 1000)
    const scalar Re = e*rhoG*u*d/mu;
    const scalar CdRe =
        Re < 1000
      ? 24*(1 + 0.15*Foam::pow(Re, 0.687))
      : 0.44*Re;

    // 0.75 Cd e (1-e) rhoG u/d e^-2.65 = 0.75 Cd Re mu (1-e)/d^2 e^-2.65
    return 0.75*CdRe*mu*(1 - e)/(d*d)*Foam::pow(e, -2.65);
}


Foam::scalar Foam::gasDragModels::ergunWenYu::dKdu
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    if (e < eSwitch)
    {
        return 1.75*(1 - e)*rhoG/d;
    }
    return gasDragModel::dKdu(e, u, d, mu, rhoG);
}


Foam::scalar Foam::gasDragModels::ergunWenYu::slip
(
    const scalar e,
    const scalar G,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    if (G <= 0 || e <= 0 || e >= 1)
    {
        return 0;
    }
    if (e >= eSwitch)
    {
        return gasDragModel::slip(e, G, d, mu, rhoG);
    }

    // K u / e = a u + b u^2 = G (Ergun), stable root
    const scalar a = 150*(1 - e)*(1 - e)*mu/(e*e*d*d);
    const scalar b = 1.75*(1 - e)*rhoG/(e*d);
    return 2*G/(a + Foam::sqrt(a*a + 4*b*G));
}


// ************************************************************************* //
