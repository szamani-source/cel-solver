#include "vonMisesModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveModels
{
    defineTypeNameAndDebug(vonMises, 0);
    addToRunTimeSelectionTable(constitutiveModel, vonMises, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveModels::vonMises::vonMises(const dictionary& dict)
:
    constitutiveModel
    (
        dict,
        typeName,
        wordList({"E", "nu", "sigma_y"}),
        wordList({"H"})
    )
{}


// ************************************************************************* //
