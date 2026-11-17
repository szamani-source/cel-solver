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
    - remap(dt)        implemented (WP 2.4, eulerianSolidRemap.C)

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
      Optional subdictionary remapLimiters (WP 2.4), for example
      remapLimiters { rho MC; U MC; sigma MC; epsP superbee; epsPEq superbee; }
      (these are the defaults; also minmod and upwind).
    - system/controlDict:
        maxCo (default 0.4, decision 0018), the fraction of the CFL
        limit (stableDeltaT);
        solveMomentum (default yes): no holds U fixed (prescribed motion;
        face velocities are then the mean of the face states, decision 0022);
        secondOrder (default yes): no uses first-order (piecewise constant)
        face states in the Lagrange step;
        gravityRamp (default 0) [s]: gravity grows from 0 as (1 - cos) over
        this time (WP 3.1, decision 0030);
        damping (default 0) [1/s]: U <- U/(1 + damping dt) after each
        momentum update, for settling (setDamping() changes it).
    - constant/g (optional, standard OpenFOAM uniformDimensionedVectorField):
      gravity [m/s^2], default zero. Balanced in the face states by the
      hydrostatic stress jump of each cell (decision 0030).
    - constant/materialProperties, optional subdictionary walls (WP 3.1):
      patch name -> { friction <deg>; } makes that patch a rigid wall: no
      normal velocity, no tension (material separates), Coulomb shear
      traction at most tan(friction) times the wall pressure (0: smooth).
    - constant/materialProperties, optional freeSurface (default no; WP 3.2a,
      decision 0031): one material plus void. Field alpha (material volume
      fraction, read if present, default 1); rho is then the bulk density
      (cell mass over cell volume, 0 in empty cells) and U, sigma, epsP,
      epsPEq are per unit material mass. Cartesian hex meshes only.
      Material with c 0 or sigmaT 0 separates instead of carrying tension.
    - constant/materialProperties, optional subdictionary outlets (WP 3.3,
      decision 0033; needs freeSurface yes): patch name ->
      { friction <deg>; openTime <s>; } makes that patch a closed wall with
      that friction until the start of the first step at or after openTime
      (default: never; openOutlets() opens them at once), then a
      traction-free outflow face (prototype rules 10 and 10b): the material
      swept across it leaves the domain and is counted in discharged().
      Opening an outlet ends the settling stage: damping is set to 0.

SourceFiles
    eulerianSolid.C
    eulerianSolidLagrange.C
    eulerianSolidRemap.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_eulerianSolid_H
#define Foam_eulerianSolid_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "constitutiveModel.H"
#include "HashTable.H"

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

        //- Remap limiter per field (rho, U, sigma, epsP, epsPEq):
        //  upwind, minmod, MC or superbee. Optional subdictionary
        //  remapLimiters in materialProperties; defaults MC for rho, U and
        //  sigma, superbee for epsP and epsPEq (decision 0013 item 5).
        HashTable<word> remapLimiter_;

        //- Gravity [m/s^2] (constant/g, zero if absent)
        vector g_;

        //- Gravity ramp time [s] (controlDict gravityRamp)
        scalar gravityRamp_;

        //- Velocity damping [1/s] (controlDict damping)
        scalar damping_;

        //- Wall patches: name -> tan(friction angle) (materialProperties walls)
        HashTable<scalar> wallFriction_;

        //- Face velocity of the last Lagrange step [m/s] (diagnostics)
        surfaceVectorField vStar_;

        //- Face force of the last Lagrange step [N]: traction times area on
        //  the owner side; on boundary faces the force of the boundary on
        //  the cell (diagnostics: wall and bottom forces)
        surfaceVectorField Fstar_;

        //- Free surface (WP 3.2a, decision 0031): materialProperties
        //  freeSurface yes; one material plus void on a Cartesian hex mesh
        Switch freeSurface_;

        //- Material volume fraction [-] (1 full, 0 empty), freeSurface only
        volScalarField alpha_;

        //- Reference (material) density [kg/m^3], materialProperties rho
        scalar rho0_;

        //- Material without tensile strength (c 0 or sigmaT 0): separates
        bool noTension_;

        //- Rule 8 (decision 0023; C++ WP 3.2b, decision 0032):
        //  materialProperties shearCap yes; at contact faces the shear jump
        //  term of the Riemann solution is scaled back so the face stress
        //  stays inside the Mohr-Coulomb surface (phi, c of the material);
        //  a capped face slips (each side sees its own tangential velocity)
        Switch shearCap_;

        //- Rule 12 (decision 0028; C++ WP 3.4b, decision 0035):
        //  materialProperties normalCap yes; at contact faces the normal
        //  jump term is capped at Mohr-Coulomb yield as rule 8 caps the
        //  shear jump; a full cell strains with its own characteristic
        //  velocity at capped faces (deviatoric part, volume kept)
        Switch normalCap_;

        //- sin(phi) and c cos(phi) of the shear cap
        scalar capSinPhi_;
        scalar capCCos_;

        //- WP 3.6c (decision 0040): materialProperties flowImpedance yes
        //  (Mohr-Coulomb with mu(I) only). In flowing cells (stress at the
        //  Mohr-Coulomb surface of phi_m >= phi) the face shear jump term
        //  uses the mu(I) tangent viscosity instead of the elastic S
        //  impedance, the Hancock predictor drops the elastic stress rate
        //  along the plastic flow direction, and the shear cap is skipped
        //  between two flowing cells and uses phi_m elsewhere
        Switch flowImpedance_;

        //- mu_s, mu_2, I0, d, rho_s of the flow impedance
        FixedList<scalar, 5> flowLaw_;

        //- Rule 9 (decision 0024; decision 0032): materialProperties
        //  spreadLayers yes (freeSurface only); a partly filled cell
        //  (alpha >= 0.1) whose alpha gradient points more along gravity
        //  than along a horizontal axis is spread along that axis
        Switch spreadLayers_;

        //- Number of faces capped by rule 8 in the last Lagrange step
        //  (this processor, owner side of coupled faces)
        label nCapped_;

        //- Number of faces capped by rule 12 in the last Lagrange step
        label nNormalCapped_;

        //- Outlet patches (WP 3.3): opening time [s]; a closed outlet is
        //  in wallFriction_, an open one is not
        HashTable<scalar> outletOpenTime_;

        //- Mass that left through non-coupled boundary faces (this
        //  processor) [kg]
        scalar discharged_;

        //- Cell lengths along x, y, z (Cartesian hex cells, freeSurface)
        vectorField hCell_;

        //- Passed from lagrangeStep to remap (freeSurface): slab side per
        //  axis (+1, -1 or 0 spread), material fraction of the Lagrangian
        //  cell, and 1 where the cell and all its neighbours are full
        //  (second order allowed)
        vectorField lean_;
        scalarField alphaL_;
        scalarField smooth_;


    // Private Member Functions

        //- Cell lengths, checking that every face is normal to an axis
        void setCartesianGeometry();

        //- Remap with void (freeSurface): slab volume fluxes, decision 0031
        void remapFreeSurface(const scalar dt);

        //- Limiters from materialProperties, with the defaults
        static HashTable<word> readRemapLimiters(const dictionary& dict);


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
        const surfaceVectorField& vStar() const { return vStar_; }
        const surfaceVectorField& Fstar() const { return Fstar_; }
        const vector& g() const { return g_; }
        scalar damping() const { return damping_; }
        const HashTable<scalar>& wallFriction() const { return wallFriction_; }
        bool freeSurface() const { return freeSurface_; }
        const volScalarField& alpha() const { return alpha_; }
        scalar rho0() const { return rho0_; }
        bool noTension() const { return noTension_; }
        bool shearCap() const { return shearCap_; }

        //- Flow impedance on (WP 3.6c)
        bool flowImpedance() const { return flowImpedance_; }
        bool spreadLayers() const { return spreadLayers_; }
        label nCapped() const { return nCapped_; }
        label nNormalCapped() const { return nNormalCapped_; }

        //- Is the patch an outlet, and is it open (WP 3.3)?
        bool isOutlet(const word& patch) const
        {
            return outletOpenTime_.found(patch);
        }
        bool outletOpen(const word& patch) const
        {
            return isOutlet(patch) && !wallFriction_.found(patch);
        }

        //- Open every outlet now (damping off)
        void openOutlets();

        //- Mass that left through non-coupled boundaries, all processors
        //  [kg]
        scalar discharged() const;

        //- Gravity at time t, with the ramp
        vector gravityAt(const scalar t) const;

        //- Change the velocity damping (for example off once settled)
        void setDamping(const scalar alpha) { damping_ = alpha; }

    // Non-const access for test programs that set a state by hand
    // (for example a prescribed phiL to test remap() alone)

        volScalarField& rhoRef() { return rho_; }
        volVectorField& URef() { return U_; }
        volSymmTensorField& sigmaRef() { return sigma_; }
        volSymmTensorField& epsPRef() { return epsP_; }
        volScalarField& epsPEqRef() { return epsPEq_; }
        surfaceScalarField& phiLRef() { return phiL_; }
        volScalarField& alphaRef() { return alpha_; }


    // Explicit cycle

        //- Stable explicit time step [s]:
        //  maxCo * min over cells of cbrt(V)/c, with c = sqrt(M/rho) the
        //  elastic P-wave speed in the cell and M the model's wave modulus.
        //  On a uniform mesh this is maxCo * min(cbrt(V)) / max(c).
        //  Parallel-safe.
        scalar stableDeltaT() const;

        //- Lagrange step over dt (WP 2.3, decision 0013, boundary
        //  conditions decision 0018): MUSCL-Hancock face states, acoustic
        //  Riemann solver per face, momentum, Hughes-Winget rotation and one
        //  constitutive call per cell. Fills phiL and leaves the Lagrangian
        //  state (decision 0014). Parallel-safe; coupled faces are solved on
        //  the owner side only, so both sides hold identical face values.
        void lagrangeStep(const scalar dt);

        //- Remap the Lagrange-step state back onto the fixed mesh (WP 2.4,
        //  eulerianSolidRemap.C): moves the volume dt*phiL back across each
        //  face (internal and coupled; across other boundary faces material
        //  leaves or enters with the boundary cell's own values).
        //  Density by volume flux, U, sigma, epsP and epsPEq mass-weighted,
        //  limited linear reconstruction at the swept-region centroid.
        //  Conserves mass and sum(m q) for every field; a cell with no flux
        //  through any face is left bit-identical (decision 0017).
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
