#include "linearElasticModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveModels
{
    defineTypeNameAndDebug(linearElastic, 0);
    addToRunTimeSelectionTable(constitutiveModel, linearElastic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveModels::linearElastic::linearElastic(const dictionary& dict)
:
    constitutiveModel
    (
        dict,
        typeName,
        wordList({"E", "nu"}),
        wordList()
    )
{}


// ************************************************************************* //
