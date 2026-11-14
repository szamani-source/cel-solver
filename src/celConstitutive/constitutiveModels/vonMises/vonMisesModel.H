/*---------------------------------------------------------------------------*\
Class
    Foam::constitutiveModels::vonMises

Description
    von Mises (J2) with linear isotropic hardening: E, nu, sigma_y, H (optional).
    Wraps cel::constitutive kernel model "vonMises".

SourceFiles
    vonMisesModel.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_constitutiveModels_vonMises_H
#define Foam_constitutiveModels_vonMises_H

#include "constitutiveModel.H"

namespace Foam
{
namespace constitutiveModels
{

class vonMises
:
    public constitutiveModel
{
public:

    //- Runtime type information
    TypeName("vonMises");

    //- Construct from dictionary
    explicit vonMises(const dictionary& dict);

    //- Destructor
    virtual ~vonMises() = default;
};

} // End namespace constitutiveModels
} // End namespace Foam

#endif

// ************************************************************************* //
