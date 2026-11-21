#include "druckerPragerModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveModels
{
    defineTypeNameAndDebug(druckerPrager, 0);
    addToRunTimeSelectionTable(constitutiveModel, druckerPrager, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveModels::druckerPrager::druckerPrager(const dictionary& dict)
:
    constitutiveModel
    (
        dict,
        typeName,
        wordList({"E", "nu", "phi", "c", "psi"}),
        wordList({"mu2", "I0", "d", "rhoS", "muIAlpha", "phi0", "phiC", "bPhi", "KDil", "nonlocalA", "capLambda", "capP0"}),
        wordList({"fit"})
    )
{}


// ************************************************************************* //
