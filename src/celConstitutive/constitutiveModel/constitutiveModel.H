/*---------------------------------------------------------------------------*\
Class
    Foam::constitutiveModel

Description
    OpenFOAM wrapper around the plain C++ constitutive kernel
    (cel::constitutive::Model, kernel/celConstitutiveKernel.H), with run-time
    selection from a dictionary:

    \verbatim
        type    vonMises;   // linearElastic | vonMises
        E       2e7;        // Pa
        nu      0.3;
        sigma_y 2e5;        // Pa (vonMises)
        H       0;          // Pa, optional hardening modulus (vonMises)
    \endverbatim

    update() takes OpenFOAM symmTensors. The strain increment is the tensor
    strain increment (e.g. symm(gradU)*dt), so its shear components are
    eps_ij, not gamma_ij; the wrapper doubles them for the kernel. Component
    order is converted between symmTensor (xx, xy, xz, yy, yz, zz) and the
    kernel's Voigt order (xx, yy, zz, xy, yz, zx).

SourceFiles
    constitutiveModel.C
    constitutiveModelNew.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_constitutiveModel_H
#define Foam_constitutiveModel_H

#include "dictionary.H"
#include "symmTensor.H"
#include "wordList.H"
#include "autoPtr.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

#include "celConstitutiveKernel.H"

#include <memory>

namespace Foam
{

class constitutiveModel
{
    // Private Data

        //- The kernel model doing the work
        std::unique_ptr<cel::constitutive::Model> kernel_;


public:

    //- History variables stored per point
    typedef cel::constitutive::History history;

    //- Runtime type information
    TypeName("constitutiveModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            constitutiveModel,
            dictionary,
            (const dictionary& dict),
            (dict)
        );


    // Constructors

        //- Read the required and optional scalar parameters from dict and
        //  build the kernel model of the given type
        constitutiveModel
        (
            const dictionary& dict,
            const word& modelType,
            const wordList& required,
            const wordList& optional
        );

        //- No copy construct
        constitutiveModel(const constitutiveModel&) = delete;

        //- No copy assignment
        void operator=(const constitutiveModel&) = delete;


    //- Select from dictionary keyword "type"
    static autoPtr<constitutiveModel> New(const dictionary& dict);


    //- Destructor
    virtual ~constitutiveModel() = default;


    // Member Functions

        //- The kernel model
        const cel::constitutive::Model& kernel() const
        {
            return *kernel_;
        }

        //- Elastic P-wave modulus K + 4G/3 [Pa]
        scalar waveModulus() const
        {
            return kernel_->waveModulus();
        }

        //- One increment at one point. deps is the tensor strain increment.
        void update
        (
            const symmTensor& deps,
            const symmTensor& sigOld,
            const history& hOld,
            symmTensor& sigNew,
            history& hNew
        ) const;


    // Conversions

        //- symmTensor to Voigt (xx, yy, zz, xy, yz, zx), shear times shearFactor
        static void toVoigt
        (
            const symmTensor& t,
            double v[6],
            const scalar shearFactor = 1
        );

        //- Voigt (tensor shear components) to symmTensor
        static symmTensor fromVoigt(const double v[6]);
};


} // End namespace Foam

#endif

// ************************************************************************* //
