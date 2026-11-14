/*---------------------------------------------------------------------------*\
Class
    Foam::constitutiveModels::mohrCoulomb

Description
    Perfectly plastic Mohr-Coulomb, non-associated: E, nu, phi [deg], c [Pa],
    psi [deg], and the optional tensile strength sigmaT [Pa] (Rankine
    tension cut-off, 0 <= sigmaT <= c cot(phi); absent = no cut-off, WP 2.5c).
    Wraps cel::constitutive kernel model "mohrCoulomb".

SourceFiles
    mohrCoulombModel.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_constitutiveModels_mohrCoulomb_H
#define Foam_constitutiveModels_mohrCoulomb_H

#include "constitutiveModel.H"

namespace Foam
{
namespace constitutiveModels
{

class mohrCoulomb
:
    public constitutiveModel
{
public:

    //- Runtime type information
    TypeName("mohrCoulomb");

    //- Construct from dictionary
    explicit mohrCoulomb(const dictionary& dict);

    //- Destructor
    virtual ~mohrCoulomb() = default;
};

} // End namespace constitutiveModels
} // End namespace Foam

#endif

// ************************************************************************* //
