/*---------------------------------------------------------------------------*\
Class
    Foam::eulerianSolid

Description
    Eulerian state of one solid material on the fixed cell-centred mesh
    (WP 2.1). Owns the fields and the constitutive model, and provides one
    method per step of the explicit cycle, so later work packages fill in a
    method rather than reshape the solver:

    - stableDeltaT()   implemented (WP 2.1)
    - lagrangeStep(dt) implemented (WP 2.3, eulerianSolidLagrange.C)
    - remap(dt)        empty, WP 2.4 (advection back onto the fixed mesh)

    Fields (cell-centred, read from the time directory if present, otherwise
    created with the default shown, written at output times):

    \verbatim
        rho     volScalarField      [kg/m^3]  default: rho in materialProperties
        U       volVectorField      [m/s]     default: 0
        sigma   volSymmTensorField  [Pa]      default: 0 (positive in tension)
        epsP    volSymmTensorField  [-]       default: 0 (plastic strain, tensor
                                              shear components)
        epsPEq  volScalarField      [-]       default: 0 (equivalent plastic strain)
    \endverbatim

    Fields created from defaults get zeroGradient patches.

    Dictionaries read:
    - constant/materialProperties: rho [kg/m^3] plus the constitutive model
      entries (type, E, nu, ...; see constitutiveModel.H).
    - system/controlDict:
        maxCo (default 0.4, decision 0016), the fraction of the CFL
        limit (stableDeltaT);
        solveMomentum (default yes): no holds U fixed (prescribed motion);
        secondOrder (default yes): no uses first-order (piecewise constant)
        face states in the Lagrange step.

SourceFiles
    eulerianSolid.C
    eulerianSolidLagrange.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_eulerianSolid_H
#define Foam_eulerianSolid_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "constitutiveModel.H"

namespace Foam
{

class eulerianSolid
{
    // Private Data

        //- The fixed Eulerian mesh
        const fvMesh& mesh_;

        //- Material dictionary (constant/materialProperties)
        IOdictionary materialProperties_;

        //- Constitutive model, selected by keyword "type"
        autoPtr<constitutiveModel> material_;

        //- Fraction of the CFL limit used for the time step
        scalar maxCo_;

        //- Update U in the Lagrange step (no: U held fixed, prescribed
        //  motion for kinematic tests)
        Switch solveMomentum_;

        //- Limited linear face states plus the half-step predictor
        //  (no: first order)
        Switch secondOrder_;

        //- Density [kg/m^3]
        volScalarField rho_;

        //- Velocity [m/s]
        volVectorField U_;

        //- Cauchy stress [Pa], positive in tension
        volSymmTensorField sigma_;

        //- Plastic strain (tensor shear components) [-]
        volSymmTensorField epsP_;

        //- Equivalent plastic strain [-]
        volScalarField epsPEq_;

        //- Face volume flux of the last Lagrange step [m^3/s]: face velocity
        //  dotted with the face area vector, positive out of the owner cell.
        //  This is the contract between lagrangeStep() and remap():
        //  after lagrangeStep(dt) the cell fields hold the values of the
        //  moved (Lagrangian) cells, whose volume is
        //      V_L = V + dt * (sum of phiL out of the cell),
        //  with rho the Lagrangian density (cell mass unchanged). remap(dt)
        //  moves the volume dt*phiL back across each face and returns all
        //  fields to the fixed mesh. Not read or written. Zero until
        //  WP 2.3 fills lagrangeStep().
        surfaceScalarField phiL_;


public:

    // Constructors

        //- Read materialProperties and the fields of the current time
        explicit eulerianSolid(const fvMesh& mesh);

        //- No copy construct
        eulerianSolid(const eulerianSolid&) = delete;

        //- No copy assignment
        void operator=(const eulerianSolid&) = delete;


    //- Destructor
    ~eulerianSolid() = default;


    // Access

        const fvMesh& mesh() const { return mesh_; }
        const constitutiveModel& material() const { return *material_; }
        scalar maxCo() const { return maxCo_; }
        bool solveMomentum() const { return solveMomentum_; }
        bool secondOrder() const { return secondOrder_; }

        const volScalarField& rho() const { return rho_; }
        const volVectorField& U() const { return U_; }
        const volSymmTensorField& sigma() const { return sigma_; }
        const volSymmTensorField& epsP() const { return epsP_; }
        const volScalarField& epsPEq() const { return epsPEq_; }
        const surfaceScalarField& phiL() const { return phiL_; }

    // Non-const access for test programs that set a state by hand
    // (for example a prescribed phiL to test remap() alone)

        volScalarField& rhoRef() { return rho_; }
        volVectorField& URef() { return U_; }
        volSymmTensorField& sigmaRef() { return sigma_; }
        volSymmTensorField& epsPRef() { return epsP_; }
        volScalarField& epsPEqRef() { return epsPEq_; }
        surfaceScalarField& phiLRef() { return phiL_; }


    // Explicit cycle

        //- Stable explicit time step [s]:
        //  maxCo * min over cells of cbrt(V)/c, with c = sqrt(M/rho) the
        //  elastic P-wave speed in the cell and M the model's wave modulus.
        //  On a uniform mesh this is maxCo * min(cbrt(V)) / max(c).
        //  Parallel-safe.
        scalar stableDeltaT() const;

        //- Lagrange step over dt (WP 2.3, decision 0013, boundary
        //  conditions decision 0016): MUSCL-Hancock face states, acoustic
        //  Riemann solver per face, momentum, Hughes-Winget rotation and one
        //  constitutive call per cell. Fills phiL and leaves the Lagrangian
        //  state (decision 0014). Parallel-safe; coupled faces are solved on
        //  the owner side only, so both sides hold identical face values.
        void lagrangeStep(const scalar dt);

        //- Remap the Lagrange-step state back onto the fixed mesh (WP 2.4).
        //  Empty in WP 2.1.
        void remap(const scalar dt);


    // Diagnostics (parallel-safe)

        //- Total mass sum(rho V) [kg]
        scalar totalMass() const;

        //- Total momentum sum(rho U V) [kg m/s]
        vector totalMomentum() const;
};


} // End namespace Foam

#endif

// ************************************************************************* //
