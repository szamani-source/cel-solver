/*---------------------------------------------------------------------------*\
Class
    Foam::constitutiveModels::druckerPrager

Description
    Perfectly plastic Drucker-Prager, non-associated: E [Pa], nu, phi [deg],
    c [Pa], psi [deg], fit (optional word: compression | extension |
    plane_strain, default compression).
    Wraps cel::constitutive kernel model "druckerPrager".

SourceFiles
    druckerPragerModel.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_constitutiveModels_druckerPrager_H
#define Foam_constitutiveModels_druckerPrager_H

#include "constitutiveModel.H"

namespace Foam
{
namespace constitutiveModels
{

class druckerPrager
:
    public constitutiveModel
{
public:

    //- Runtime type information
    TypeName("druckerPrager");

    //- Construct from dictionary
    explicit druckerPrager(const dictionary& dict);

    //- Destructor
    virtual ~druckerPrager() = default;
};

} // End namespace constitutiveModels
} // End namespace Foam

#endif

// ************************************************************************* //
