#include "mohrCoulombModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace constitutiveModels
{
    defineTypeNameAndDebug(mohrCoulomb, 0);
    addToRunTimeSelectionTable(constitutiveModel, mohrCoulomb, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveModels::mohrCoulomb::mohrCoulomb(const dictionary& dict)
:
    constitutiveModel
    (
        dict,
        typeName,
        wordList({"E", "nu", "phi", "c", "psi"}),
        wordList()
    )
{}


// ************************************************************************* //
