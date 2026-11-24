/*---------------------------------------------------------------------------*\
    gasDragModel: run-time selected gas-grain drag (WP F.2, decision 0127)
\*---------------------------------------------------------------------------*/

#include "gasDragModel.H"
#include "error.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gasDragModel, 0);
    defineRunTimeSelectionTable(gasDragModel, dictionary);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::gasDragModel> Foam::gasDragModel::New
(
    const word& name,
    const dictionary& dict
)
{
    auto* ctorPtr = dictionaryConstructorTable(name);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "gasDragModel",
            name,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<gasDragModel>(ctorPtr(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::gasDragModel::dKdu
(
    const scalar e,
    const scalar u,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    // Step from the scale of the Darcy-to-inertial transition, Re ~ 1
    const scalar h = 1e-6*max(u, mu/(rhoG*d));
    if (u < h)
    {
        return (K(e, u + h, d, mu, rhoG) - K(e, u, d, mu, rhoG))/h;
    }
    return (K(e, u + h, d, mu, rhoG) - K(e, u - h, d, mu, rhoG))/(2*h);
}


Foam::scalar Foam::gasDragModel::slip
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

    // f(u) = K(e, u) u / e - G grows with u; f(0) = -G < 0 and, since
    // for most laws K(e, u) >= K(e, 0), so f(uD) >= 0 at the Darcy speed uD
    auto f = [&](const scalar u) { return K(e, u, d, mu, rhoG)*u/e - G; };

    scalar a = 0, fa = -G;
    scalar b = G*e/K(e, 0, d, mu, rhoG), fb = f(b);
    // A law whose K falls with the slip (syamlalOBrien, WP F.7) can leave
    // the root above the Darcy speed: widen the bracket
    for (int it = 0; fb < 0 && it < 200; ++it)
    {
        a = b;
        fa = fb;
        b *= 2;
        fb = f(b);
    }
    if (fb == 0)
    {
        return b;
    }

    // Illinois (modified regula falsi)
    int side = 0;
    for (int it = 0; it < 200; ++it)
    {
        const scalar c = (a*fb - b*fa)/(fb - fa);
        const scalar fc = f(c);
        if (fc == 0 || mag(b - a) <= 1e-14*mag(c))
        {
            return c;
        }
        if ((fc > 0) == (fb > 0))
        {
            b = c;
            fb = fc;
            if (side == -1)
            {
                fa *= 0.5;
            }
            side = -1;
        }
        else
        {
            a = c;
            fa = fc;
            if (side == 1)
            {
                fb *= 0.5;
            }
            side = 1;
        }
    }
    return 0.5*(a + b);
}


Foam::scalar Foam::gasDragModel::permeability
(
    const scalar e,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    if (e <= 0)
    {
        return 0;
    }
    if (e >= 1)
    {
        return GREAT;
    }
    return mu*e*e/K(e, 0, d, mu, rhoG);
}


Foam::scalar Foam::gasDragModel::mobility
(
    const scalar e,
    const scalar G,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    if (e <= 0)
    {
        return 0;
    }
    if (e >= 1)
    {
        return GREAT/mu;
    }
    return e*e/K(e, slip(e, G, d, mu, rhoG), d, mu, rhoG);
}


Foam::scalar Foam::gasDragModel::tangentMobility
(
    const scalar e,
    const scalar G,
    const scalar d,
    const scalar mu,
    const scalar rhoG
) const
{
    if (e <= 0)
    {
        return 0;
    }
    if (e >= 1)
    {
        return GREAT/mu;
    }
    // G = K(u) u / e, so du/dG = e/(K + u K'(u)) and d|q|/dG = e du/dG
    const scalar u = slip(e, G, d, mu, rhoG);
    return e*e/(K(e, u, d, mu, rhoG) + u*dKdu(e, u, d, mu, rhoG));
}


// ************************************************************************* //
