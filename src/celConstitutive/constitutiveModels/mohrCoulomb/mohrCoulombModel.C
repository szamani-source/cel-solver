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
        wordList({"sigmaT", "mu2", "I0", "d", "rhoS", "muIAlpha", "phi0", "phiC", "bPhi", "KDil", "nonlocalA", "capLambda", "capP0"})
    )
{}


// ************************************************************************* //
