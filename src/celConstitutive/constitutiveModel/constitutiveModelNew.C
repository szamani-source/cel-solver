#include "constitutiveModel.H"
#include "error.H"
#include "celMuI.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::constitutiveModel> Foam::constitutiveModel::New
(
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Constitutive model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "constitutiveModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    autoPtr<constitutiveModel> model(ctorPtr(dict));

    // WP 3.6a: report the mu(I) friction if it is on
    if (const auto* law = model->kernel().frictionLaw())
    {
        Info<< "Friction mu(I): mu_s " << law->muS() << ", mu_2 " << law->mu2();
        if (law->alpha() > 0)
        {
            Info<< ", Barker-Gray alpha " << law->alpha()
                << ", I_1^N " << law->I1N() << endl;
        }
        else
        {
            Info<< ", no regularisation" << endl;
        }
    }

    return model;
}


// ************************************************************************* //
