/*---------------------------------------------------------------------------*\
Class
    Foam::hex8Body

Description
    A deformable Lagrangian body of 8-node hexahedra, explicit in time
    (WP 5.1, decision 0070).

    Element: one-point (uniform strain) hexahedron of Flanagan and
    Belytschko (1981). The gradient operator is the mean of the shape
    function gradients over the element,

        B_I = (1/V) integral of grad N_I dV = (1/V) sum over faces of
              integral of N_I n dA,

    with each face integral taken exactly (2x2 Gauss points on the bilinear
    face), so sum_I x_I B_I = I for any element shape and the element
    passes the patch test. Internal force f_I = V sigma . B_I.

    Hourglass control: Flanagan-Belytschko base vectors gamma, orthogonal to
    every linear velocity field; stiffness form (as LS-DYNA type 4): the
    hourglass forces grow with the integrated hourglass velocity,
        Q_a += dt k q_a,  q_a = sum_I v_I gamma_aI,  f_I = -sum_a Q_a gamma_aI,
    k = hourglass M V sum_I |B_I|^2 / 8 (M the P-wave modulus,
    hourglass default 0.05).

    Time integration: central difference with staggered velocities.
    step(dt): forces at x_n from sigma_n; v += dt a (constrained components
    set to their prescribed velocity); x += dt v; the velocity gradient on
    the mid-step configuration gives the strain increment and spin; stress
    rotated (Hughes-Winget) and updated by the constitutive model
    (celConstitutive, the same models as the Eulerian side). Lumped mass.

    Dictionary:

        density     7800;               // kg/m^3
        material    { type vonMises; E 2e11; nu 0.3; sigma_y 4e8; }
        hourglass   0.05;               // optional
        hourglassType stiffness;        // optional: stiffness (default,
                                        // conserves energy; elastic bodies)
                                        // or viscous (LS-DYNA type 1: large
                                        // plastic flow, decision 0071)
        planeStrain yes;                // optional: no z velocity (2D runs,
                                        // one element thick)
        damping     0.1;                // optional, default 0: viscous
                                        // stress 2 damping rho c L D (L the
                                        // element size); damps ringing of
                                        // coupled bodies (decision 0072)

SourceFiles
    hex8Body.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_hex8Body_H
#define Foam_hex8Body_H

#include "constitutiveModel.H"
#include "pointField.H"
#include "vectorField.H"
#include "symmTensorField.H"
#include "FixedList.H"
#include "List.H"
#include "dictionary.H"

namespace Foam
{

class hex8Body
{
public:

    typedef FixedList<label, 8> hexElement;


private:

    // Private Data

        word name_;

        //- Material and density
        autoPtr<constitutiveModel> model_;
        scalar rho_;
        scalar hourglass_;
        bool planeStrain_;

        //- Viscous hourglass control (else stiffness) and the wave speed
        bool viscousHg_;
        scalar cWave_;

        //- Damping (decision 0072): viscous stress damping rho c L D, as a
        //  fraction of critical for the shortest mode; and that stress
        scalar damping_;
        symmTensorField sigmaV_;

        //- Elements (OpenFOAM hex node order)
        List<hexElement> elems_;

        //- Nodal positions, velocities (half step), lumped masses
        pointField x_;
        vectorField v_;
        scalarField m_;

        //- Constrained velocity components per node (1 fixed, 0 free) and
        //  their prescribed values
        vectorField fixMask_;
        vectorField fixVel_;

        //- External nodal forces (applied each step) and gravity
        vectorField fExt_;
        vector g_;

        //- Element state: stress, history, hourglass forces (4 modes)
        symmTensorField sigma_;
        List<constitutiveModel::history> hist_;
        List<FixedList<vector, 4>> Qhg_;

        //- Element volumes at the start
        scalarField V0_;

        //- Internal work (stress and hourglass) summed over the run [J]
        scalar workInt_;
        scalar workExt_;


    // Private Member Functions

        //- Volume and gradient operator of element e at positions x
        scalar gradients
        (
            const hexElement& e,
            const pointField& x,
            FixedList<vector, 8>& B
        ) const;

        //- Internal (stress + hourglass) forces at the current state
        void internalForces(vectorField& f) const;


public:

    //- Runtime type information
    ClassName("hex8Body");


    // Constructors

        //- From nodes and elements (OpenFOAM hex node order; negative
        //  volume elements are mirrored) and the body's dictionary
        hex8Body
        (
            const word& name,
            const pointField& points,
            const List<hexElement>& elems,
            const dictionary& dict
        );


    // Member Functions

        const word& name() const noexcept { return name_; }
        const pointField& points() const noexcept { return x_; }
        const vectorField& velocity() const noexcept { return v_; }
        const scalarField& mass() const noexcept { return m_; }
        const List<hexElement>& elements() const noexcept { return elems_; }
        const symmTensorField& sigma() const noexcept { return sigma_; }
        const List<constitutiveModel::history>& history() const noexcept
        {
            return hist_;
        }

        //- Set the nodal velocities (initial condition)
        void setVelocity(const vectorField& v);

        //- Prescribe velocity components of node i: mask components 1 are
        //  held at the value of vel
        void fix(const label i, const vector& mask, const vector& vel);

        //- External nodal forces, applied every step until changed
        vectorField& externalForce() noexcept { return fExt_; }

        void setGravity(const vector& g) { g_ = g; }

        //- Set element stresses (initial condition)
        void setStress(const symmTensor& s) { sigma_ = s; }

        //- Current element volumes
        tmp<scalarField> volumes() const;

        //- Internal forces at the current state (stress and hourglass;
        //  force of the elements on the nodes, so minus the residual)
        tmp<vectorField> forces() const;

        //- Stable time step: smallest V / largest face area over c
        scalar stableDeltaT() const;

        //- One explicit step
        void step(const scalar dt);

        //- Kinetic energy [J]
        scalar kineticEnergy() const;

        //- Energy held by the hourglass forces, sum |Q|^2 / 2k [J]
        scalar hourglassEnergy() const;

        //- Internal work done on the elements since the start [J]
        scalar internalWork() const noexcept { return workInt_; }

        //- Work done by external forces, gravity and prescribed velocities
        //  since the start [J]
        scalar externalWork() const noexcept { return workExt_; }

        //- Write the state a restart needs (WP 6.2, decision 0081): nodal
        //  positions and velocities, element stresses, history, hourglass
        //  forces and damping stresses, work sums
        void writeState(dictionary& dict) const;

        //- Read the state written by writeState
        void readState(const dictionary& dict);
};


} // End namespace Foam

#endif
