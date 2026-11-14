/*---------------------------------------------------------------------------*\
Class
    Foam::constitutiveModels::linearElastic

Description
    Linear elastic: E, nu.
    Wraps cel::constitutive kernel model "linearElastic".

SourceFiles
    linearElasticModel.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_constitutiveModels_linearElastic_H
#define Foam_constitutiveModels_linearElastic_H

#include "constitutiveModel.H"

namespace Foam
{
namespace constitutiveModels
{

class linearElastic
:
    public constitutiveModel
{
public:

    //- Runtime type information
    TypeName("linearElastic");

    //- Construct from dictionary
    explicit linearElastic(const dictionary& dict);

    //- Destructor
    virtual ~linearElastic() = default;
};

} // End namespace constitutiveModels
} // End namespace Foam

#endif

// ************************************************************************* //
