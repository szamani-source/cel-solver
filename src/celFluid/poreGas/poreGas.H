/*---------------------------------------------------------------------------*\
Class
    Foam::poreGas

Description
    Gas in the pores of the granular material (version 2, WP F.1, decision
    0125 and 0126). F.1 is the skeleton: it reads constant/gasProperties and
    computes, after every remap, per cell

        phiS = alpha rho_m / rhoP = rho / rhoP     solids volume fraction
        phiG = 1 - phiS                            gas volume fraction
        k                                          permeability of the
                                                   filled part [m^2]

    where rho is the cell density of eulerianSolid (mass over the open cell
    volume, so alpha rho_m with rho_m = rho/alpha the material density; Sia's
    phi_s = alpha rho_bulk / rho_particle). k is the zero-slip permeability
    of the filled part of the cell, e_m = 1 - rho/(alpha rhoP):
    k = mu e_m^2 / K(e_m, slip 0) with K the Ergun-Wen-Yu (Gidaspow)
    momentum exchange coefficient, which below e_m = 0.8 is Ergun's
    k = e^3 d^2 / (150 (1 - e)^2). Void cells (alpha below alphaMin) and
    cells packed past rhoP have k = 0 (the void holds the boundary pressure,
    F.4).

    WP F.2 (decision 0127) solves the pore pressure each step, implicitly,
    with the solid frozen (no grain velocity, no force on the solid; F.3):

        d(phi_g rho_g)/dt + div(rho_g q) = 0,   q = -lambda grad p,
        rho_g = (p0 + p)/(R T),  lambda = e^2/K(e, slip)

    The field p is the gauge pressure (absolute minus p0), which keeps the
    round-off of the solve far below the pressure differences.
    with K from the run-time selected gasDragModel (drag, default
    ergunWenYu). The slip follows from |grad p| at each face, so Ergun's
    inertial term enters through lambda. Each step is solved by Newton
    iterations on the correction: the residual of the mass balance with the
    flux law (secant mobility lambda) at the current p, the matrix with the
    tangent mobility d|q|/d|grad p| (rho_g lagged), Anderson-accelerated,
    until the largest correction is below pTolerance. The step ends with
    p = p_old - dt div(F)/S from the face fluxes F of that state, so the
    gas mass changes exactly by the boundary fluxes. Face coefficients: harmonic mean of the cell mobilities
    (0 if either side has none) times the linear face gas density, times the
    open face fraction with bodies. The equation is conservative: the gas
    mass changes only by the boundary fluxes (printed each step).

    Boundaries, from 0/p (optional; without it p = 0 and every patch is
    sealed): fixedValue = fixed pressure, zeroGradient = sealed, and
    fixedGradient patches listed in `inflow` get the gradient that blows
    gas in at the given superficial velocity (at p0 and T).

    constant/gasProperties:
    \verbatim
        model   poreGas;        // none | poreGas (none: as if absent)
        rhoP    2650;           // grain (particle) density [kg/m^3]
        d       0.001;          // grain diameter [m]
        T       293.15;         // gas temperature [K], default 293.15
        p0      101325;         // reference (ambient) pressure [Pa], default
        drag    ergunWenYu;     // run-time selected drag law (gasDragModel);
                                // F.7 adds the others (decision 0125 item 7)

        // WP F.2, all optional:
        inflow                  // fixedGradient patches of 0/p blowing gas
        {                       // in: superficial velocity [m/s] at p0, T
            bottom  0.05;
        }
        nOuter          50;     // Newton iterations per step, at most
        nAnderson       5;      // Anderson acceleration depth (0 off)
        pTolerance      1e-9;   // [Pa] largest correction when converged
                                // (default 1e-14 p0)
        laplacian       "Gauss linear corrected";
        solver                  // linear solver for each correction
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-30;
            relTol          0.1;    // inexact Newton
            maxIter         1000;
        }
    \endverbatim

    Air: ideal gas R = 287.05 J/(kg K), Sutherland viscosity, as
    tools/gas_reference.py.

SourceFiles
    poreGas.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_poreGas_H
#define Foam_poreGas_H

#include "eulerianSolid.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "autoPtr.H"
#include "surfaceFields.H"
#include "HashTable.H"
#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class poreGas
{
    // Private data

        const eulerianSolid& solid_;

        //- constant/gasProperties
        IOdictionary dict_;

        //- Grain density [kg/m^3], grain diameter [m]
        scalar rhoP_;
        scalar d_;

        //- Gas temperature [K] and reference pressure [Pa]
        scalar T_;
        scalar p0_;

        //- Drag law name and the law (WP F.2)
        word drag_;
        autoPtr<gasDragModel> dragModel_;

        //- Gas density at p0 and T [kg/m^3], viscosity at T [Pa s]
        scalar rhoG_;
        scalar mu_;

        //- Solids and gas volume fractions, permeability of the filled part
        volScalarField phiS_;
        volScalarField phiG_;
        volScalarField k_;

        //- Gas fraction of the filled part e_m = 1 - rho/(alpha rhoP)
        volScalarField eM_;

        // WP F.2: pore pressure

            //- Gauge gas pressure p_abs - p0 [Pa], from 0/p if present
            autoPtr<volScalarField> pPtr_;

            //- Superficial gas flux relative to the grains [m/s]
            volVectorField qGas_;

            //- Gas mass flux through the faces of the last solve [kg/s]
            surfaceScalarField massFlux_;

            //- Inflow patches: superficial velocity at p0 [m/s]
            HashTable<scalar> inflow_;

            //- Newton iterations at most; converged when the largest
            //  correction is below pTol [Pa]
            label nOuter_;
            scalar pTol_;

            //- Depth of the Anderson acceleration of the iterations
            label nAnderson_;

            //- Laplacian scheme and linear solver controls
            string laplacian_;
            dictionary solverDict_;

            //- Gas mass that came in through the boundaries so far [kg]
            scalar netIn_;


public:

    //- Gas constant of air [J/(kg K)]
    static constexpr scalar Rair = 287.05;

    //- Cells with alpha below this are void (no permeability)
    static constexpr scalar alphaMin = 1e-6;

    //- Ergun-Wen-Yu switch: Ergun below this gas fraction (Gidaspow)
    static constexpr scalar eSwitch = 0.8;


    // Constructors

        //- Read gasProperties (which must exist) and compute the fractions
        poreGas(const eulerianSolid& solid, const IOobject& io);

        //- The pore gas of the case: null if constant/gasProperties is
        //  absent or its model is none (v1 behaviour, nothing changes)
        static autoPtr<poreGas> New(const eulerianSolid& solid);

        poreGas(const poreGas&) = delete;
        void operator=(const poreGas&) = delete;


    // Gas properties (pure functions, shared with the tests)

        //- Air viscosity, Sutherland [Pa s]
        static scalar sutherland(const scalar T);

        //- Zero-slip permeability [m^2] of a packing with gas fraction e,
        //  grain diameter d: mu e^2 / K(e, 0) of Ergun-Wen-Yu (independent
        //  of mu). Below eSwitch e^3 d^2 / (150 (1-e)^2) (Ergun), above
        //  d^2 e^4.65 / (18 (1-e)) (Wen and Yu, Stokes limit). 0 for e <= 0.
        //  (F.1; the solver now asks the selected gasDragModel.)
        static scalar permeability(const scalar e, const scalar d);


    // Access

        scalar rhoP() const noexcept { return rhoP_; }
        scalar d() const noexcept { return d_; }
        scalar T() const noexcept { return T_; }
        scalar p0() const noexcept { return p0_; }
        scalar rhoG() const noexcept { return rhoG_; }
        scalar mu() const noexcept { return mu_; }
        const word& drag() const noexcept { return drag_; }

        const volScalarField& phiS() const noexcept { return phiS_; }
        const volScalarField& phiG() const noexcept { return phiG_; }
        const volScalarField& k() const noexcept { return k_; }
        const volScalarField& p() const noexcept { return *pPtr_; }
        const gasDragModel& dragModel() const noexcept { return *dragModel_; }

        //- Gas mass in the domain, all processors [kg]
        scalar gasMass() const;

        //- Gas mass that came in through the boundaries so far [kg]
        scalar netIn() const noexcept { return netIn_; }


    // Cycle

        //- Fractions and permeability from the solid's current rho and
        //  alpha; call after remap() (and after moving walls). The fields
        //  are written with the case's other fields.
        void correctFractions();

        //- Face coefficients rho_g lambda_f of the flux law (secant
        //  mobility, DfS) and of its derivative (tangent mobility, DfT),
        //  each side's mobility taken at the face pressure gradient
        void faceCoefficients
        (
            const volVectorField& gradP,
            surfaceScalarField& DfS,
            surfaceScalarField& DfT
        ) const;

        //- Pore pressure over the step dt, solid frozen (WP F.2); call
        //  after correctFractions()
        void solve(const scalar dt);
};


} // End namespace Foam

#endif

// ************************************************************************* //
