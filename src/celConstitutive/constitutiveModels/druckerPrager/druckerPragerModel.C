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
        wordList(),
        wordList({"fit"})
    )
{}


// ************************************************************************* //
