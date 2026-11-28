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
    until the error left, the largest correction (or, once the corrections
    shrink tenfold per iteration, the next one estimated from their ratio;
    WP F.3), is below pTolerance. The first guess is p extrapolated from the
    last two steps (WP F.3). The step ends with
    p = p_old - dt div(F)/S from the face fluxes F of that state, so the
    gas mass changes exactly by the boundary fluxes. Face coefficients: harmonic mean of the cell mobilities
    (0 if either side has none) times the linear face gas density, times the
    open face fraction with bodies. The equation is conservative: the gas
    mass changes only by the boundary fluxes (printed each step).

    WP F.3 (decision 0128) couples the two both ways (coupled yes, the
    default): sigma of eulerianSolid is the effective stress, and each
    cell's momentum also gets -sum(p_f S_f) over its faces, so the total
    stress sigma - p I balances the loads; the gas moves with the grains,

        d(phi_g rho_g)/dt + div(phi_g rho_g u_s) + div(rho_g q) = 0,

    discretised to match the solid's remap: the gas volume crossing a face
    in a step is the Lagrangian volume dt phiL less the grain volume
    dm/rhoP of the remap, at the donor's gas density of the start of the
    step, and the change of phi_g in each cell comes from the solid's new
    density. A compressed cell therefore squeezes its gas by exactly its
    Lagrangian volume change (undrained limit p0 dV/(phi_g V)).

    WP F.4 (decision 0129): void, free surface and outlets. With a
    freeSurface solid, cells with alpha below alphaVoid (1e-4) are void:
    they hold the void pressure (p = 0 gauge, the ambient p0) and exchange
    gas with the atmosphere, which the mass bookkeeping counts as inflow
    (voidIn). A face between a gas cell and a void cell carries the
    gas cell's mobility over the distance from its centre to the face (a
    pressure face, like a fixedValue boundary), and the solid sees the
    void pressure on it (pFace). Partly filled cells keep k of their filled
    part and store the gas of the whole open volume. Gas crosses outlets
    and inlets with the material (dt phiL - dm/rhoP, at the boundary
    pressure when it comes in). Outlets (materialProperties outlets) must be
    fixedValue in 0/p (without 0/p they are fixedValue 0) and are sealed
    for the gas until they open.

    WP F.6 (decision 0131): with model gasFlow the void cells carry a gas
    velocity Ug and their own pressure (low Mach, implicit, the solid's
    time step): momentum predictor, then nCorr pressure corrections that
    are the same Newton solve with rho_f rAU_f on the faces between void
    cells and the H/A mass flux (Rhie-Chow with ddtCorr) as a source. Faces
    between bed and void keep the F.4 rule against the void cell's own
    pressure; cells with grains hold the grain velocity plus q/phi_g. Open
    patches (fixedValue p) are the atmosphere. Needs a freeSurface solid.
    Implementation in poreGasFlow.C.

    WP F.7 (decision 0132): with model gasFlow the gas in the pores has
    momentum too (poreInertia yes, the default): e rho_g (w - w_o)/dt =
    -e grad p - K w for the slip w, the drag implicit, solved in the same
    Newton iterations as the drag law K + e rho_g/dt against the gradient
    less rho_g q_o/(e dt) plus the old flux's face mass flux (faces between
    gas cells, gas and void, and fixed-pressure boundaries; inflow patches
    keep Darcy). The grains feel -sum(p_f S_f) as before; the reaction of the
    pore gas inertia (inertiaForce, default no) is left out.
    Implementation in poreGasFlow.C.

    Boundaries, from 0/p (optional; without it p = 0 and every patch is
    sealed): fixedValue = fixed pressure, zeroGradient = sealed, and
    fixedGradient patches listed in `inflow` get the gradient that blows
    gas in at the given superficial velocity (at p0 and T).

    constant/gasProperties:
    \verbatim
        model   poreGas;        // none | poreGas | gasFlow (WP F.6:
                                // gas momentum in the void)
        rhoP    2650;           // grain (particle) density [kg/m^3]
        d       0.001;          // grain diameter [m]
        T       293.15;         // gas temperature [K], default 293.15
        p0      101325;         // reference (ambient) pressure [Pa], default
        drag    ergunWenYu;     // run-time selected drag law (gasDragModel):
                                // ergunWenYu (default, every test), ergun,
                                // wenYu, diFelice, syamlalOBrien, beetstra
                                // (WP F.7, decision 0132)

        coupled yes;            // WP F.3: two-way with the grains
                                // (no: F.2, the solid frozen for the gas)
        alphaVoid 1e-4;         // WP F.4: void below this alpha

        // WP F.2, all optional:
        inflow                  // fixedGradient patches of 0/p blowing gas
        {                       // in: superficial velocity [m/s] at p0, T,
            bottom  0.05;       // or a Function1 of time (WP F.7), e.g.
                                // bottom table ((0 0) (1 0.05));
        }
        nOuter          50;     // Newton iterations per step, at most
        nAnderson       5;      // Anderson acceleration depth (0 off)
        pTolerance      1e-4;   // [Pa] error left when converged (decision 0141)
                                // (default 1e-14 p0)
        laplacian       "Gauss linear corrected";
        parallelCheck   no;     // WP F.9a (decision 0134), debug: compare
                                // every gas and solid field across each
                                // processor cut (stale halos, one-sided
                                // face values); prints "Parallel check"
        parallelCheckEvery 1;   // ... every this many steps
        liftCheck       yes;    // decision 0135: report the share of the
                                // grain mass whose weight the gas carries
                                // (a lifted or expanding bed is outside the
                                // model's range); prints "Gas lift" at
                                // write times, warns once past liftWarn
        liftRatio       1;      // ... a cell is carried when the upward gas
                                // force is at least this times its weight
        liftWarn        0.01;   // ... warn when this share is carried
        liftTime        0.01;   // ... for this long [s]

        // WP F.6, model gasFlow only, all optional:
        nCorr           1;      // pressure corrections per step (PISO)
        poreInertia     yes;    // WP F.7: gas momentum in the pores (its
                                // inertia, drag implicit; no = F.6 Darcy)
        inertiaForce    no;     // WP F.7: its reaction on the grains
                                // (unstable as a bed lifts, decision 0132)
        divUg           "Gauss upwind";     // convection of Ug
        USolver                 // linear solver of the gas momentum
        {
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-30;
            relTol          1e-10;
        }
        solver                  // linear solver for each correction
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-30;
            relTol          1e-3;   // inexact Newton (F.2: 0.1)
            maxIter         1000;
        }
    \endverbatim

    Field 0/Ug (gasFlow, optional): gas velocity; without it no slip on
    walls, zeroGradient where p is fixedValue.

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
#include "HashPtrTable.H"
#include "Function1.H"
#include "gasDragModel.H"
#include "Switch.H"
#include "fvMatrices.H"

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

            //- Inflow patches: superficial velocity at p0 [m/s], a
            //  Function1 of time (WP F.7: ramps)
            HashPtrTable<Function1<scalar>> inflow_;

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

        // WP F.3: coupling with the grains (decision 0128)

            //- Two-way coupling (default yes): the gas pressure acts on the
            //  solid (sigma is the effective stress) and the grains carry
            //  and squeeze the gas. no: the F.2 frozen-solid gas.
            Switch coupled_;

            //- Gas fraction at the last solve (before this step's remap)
            scalarField phiG0_;

            //- Pressure before the last step and that step's length, for
            //  the first guess of the next solve (linear extrapolation)
            scalarField pOldOld_;
            scalar dtOld_;

        // WP F.4: void, free surface and outlets (decision 0129)

            //- Free-surface solid: void cells hold the void pressure
            bool freeSurface_;

            //- Cells with alpha below this are void for the gas (default
            //  1e-4, the v1 void threshold, decision 0129)
            scalar alphaVoid_;

            //- 1 in void cells (alpha below alphaVoid_), else 0; its
            //  processor-neighbour values mark void across the cuts
            volScalarField void_;

            //- Gas pressure on the faces as the solid sees it: linear,
            //  but the void pressure on faces between material and void
            surfaceScalarField pf_;

            //- Gas mass the void (atmosphere) gave this run [kg]
            scalar voidIn_;

        // WP F.6: gas momentum in the void (model gasFlow, decision 0131)

            //- Gas flow in the void: void cells carry a gas velocity and
            //  their own pressure instead of holding the void pressure
            bool flow_;

            //- Pressure corrections per step (PISO), at least 1
            label nCorr_;

            //- Convection scheme of the gas momentum
            string divUg_;

            //- Linear solver of the gas momentum
            dictionary USolverDict_;

            //- Gas velocity Ug [m/s]: solved in void cells; in the others
            //  the grain velocity plus the pore flux over phi_g
            autoPtr<volVectorField> UgPtr_;

            //- Momentum matrix of this step and its 1/A, H/A
            autoPtr<fvVectorMatrix> UEqnPtr_;
            autoPtr<volScalarField> rAUPtr_;
            autoPtr<volVectorField> HbyAPtr_;

            //- 1/A at the faces [m^3 s/kg], and the gas mass flux of H/A
            //  through faces of void cells [kg/s]
            autoPtr<surfaceScalarField> rAUfPtr_;
            autoPtr<surfaceScalarField> mHbyAPtr_;

            //- Velocity and volume flux at the start of the step (ddtCorr)
            autoPtr<volVectorField> Uold_;
            autoPtr<surfaceScalarField> phiOld_;

        // WP F.7: gas momentum in the pores, drag on the grains (decision
        // 0132)

            //- Inertia of the pore gas (gasFlow only; default yes there)
            Switch poreInertia_;

            //- The drag law with the pore gas inertia of this step
            //  (K + e rhoG/dt); null when off
            autoPtr<gasDragModel> inertialDrag_;

        // WP F.8: gas carried across the free surface (decision 0133)

            //- The surface's own motion carries gas across faces between a
            //  gas cell and a void cell (the vacated or newly filled part
            //  of the cell exchanges gas with the void), as on every other
            //  face; no: only grains crossing displace gas there (F.4)
            Switch surfaceCarry_;

            //- Smallest gas fraction of the filled part (WP F.8, F.8c,
            //  decision 0133): keyword eMin, default 0.2
            scalar eMin_;

        // WP F.8e: pore pressure as an advected-diffused field (decision
        // 0133, route B; Iverson & Denlinger 2001, Gueugneau et al. 2017,
        // Aravena et al. 2021)

            //- gasProperties diffusivity [m^2/s]; 0 (default): off. When
            //  set, the pore pressure moves with the grains and diffuses
            //  with this diffusivity: the storage is rho_g0 k/(mu kappa)
            //  per unit of filled volume, and the grains' volume change
            //  does not act on it (no compression or suction)
            scalar kappa_;

            //- gasProperties frictionOnly (WP F.8g, needs diffusivity):
            //  the pore pressure only lowers the friction of the grains
            //  (Gueugneau 2017, Aravena 2021) instead of acting as a force
            Switch frictionOnly_;

            //- WP F.8h: total mean pressure of each cell at the last solve
            //  (frictionOnly carries the support ratio p/P); empty: not set
            scalarField Pprev_;

            //- Superficial flux relative to the grains through each face
            //  at the end of the last step, q.n [m/s], and 1 on the faces
            //  that were void flow faces then (their flux is absolute)
            autoPtr<surfaceScalarField> qnOld_;
            autoPtr<surfaceScalarField> voidFaceOld_;

            //- This step's old momentum at the faces: rho_g q_o.n/(e dt)
            //  |Sf| [kg/(m^2 s^2)], which times the face coefficient is
            //  its mass flux, and the old flux as a pressure gradient,
            //  rho_g q_o/(e dt) [Pa/m] (normal part from the face flux)
            autoPtr<surfaceScalarField> sInertia_;
            autoPtr<surfaceVectorField> gInertia_;

            //- Old relative flux per cell (qGas at the start of the step)
            //  and the void flag then
            vectorField qCellOld_;
            scalarField voidOld_;

            //- Mass flux the old pore-gas momentum carries through the
            //  faces at the current face coefficients [kg/s]
            autoPtr<surfaceScalarField> mInertia_;

            //- Reaction of the pore gas inertia on the grains, per cell
            //  [N]: -V (rho_g/dt) (q - q_o); the drag K w is the pressure
            //  part -phi_g grad p less this (decision 0132)
            autoPtr<volVectorField> gasForce_;

        // WP F.9: parallel consistency check (debug)

            //- gasProperties parallelCheck (default no): after each solve,
            //  compare every field's processor-boundary values with its
            //  neighbour's (stale halo values, one-sided face values)
            Switch parallelCheck_;

            //- Every this many steps (gasProperties parallelCheckEvery,
            //  default 1)
            label parallelCheckEvery_;

            //- Largest mismatch seen so far, and where
            scalar parallelWorst_;

        // Range of validity (decision 0135): the gas carrying the grains

            //- gasProperties liftCheck (default yes): each step, the share
            //  of the grain mass whose weight the gas carries
            //  (poreGasLift.C)
            Switch liftCheck_;

            //- A cell counts as carried when the upward gas force on it is
            //  at least liftRatio times its weight (default 1)
            scalar liftRatio_;

            //- Warn once the carried share of the grain mass passes this
            //  (gasProperties liftWarn, default 0.01)
            scalar liftWarn_;

            //- ... and stays past it for this long [s] (gasProperties
            //  liftTime, default 0.01): a short spike (an outlet opening)
            //  is not a lift
            scalar liftTime_;

            //- Time the share last went past liftWarn (-1: not past it)
            scalar liftSince_;

            //- The warning was given
            bool liftWarned_;

            //- Carried share of the grain mass, last step
            scalar liftShare_;


    // Private Member Functions (WP F.6, poreGasFlow.C)

        //- Create Ug and the flow controls (constructor)
        void initFlow();

        //- 1 on the faces of the void flow (both sides void; open
        //  boundary faces of void cells), else 0
        tmp<surfaceScalarField> voidFaces() const;

        //- Momentum matrix of the step (predictor solved), from the old
        //  velocity, the last mass flux and the pressure p
        void predictMomentum(const scalar dt, const scalarField& pOld);

        //- rAU, H/A and its mass flux through faces of void cells
        void flowFluxes();

        //- Gradient of the face pressures pFace (the void cell's own
        //  pressure on faces between bed and void), per cell
        tmp<vectorField> flowGradP();

        //- Void-cell velocities from H/A and the pressure gradient; the
        //  others from the grains and the pore flux
        void correctVelocity();

    // Private Member Functions (WP F.7, poreGasFlow.C)

        //- Old relative face fluxes and the inertial drag law of the step
        void startInertia(const scalar dt);

        //- mInertia_ at the face coefficients DfS
        void inertiaFlux(const surfaceScalarField& DfS);

        //- Store this step's face fluxes and the force on the grains
        void endInertia(const scalar dt);


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

        //- Gas mass flux through the faces of the last solve, relative to
        //  the grains (Darcy-Ergun, and in the void Ug) [kg/s], out of the
        //  owner (WP T.4: the gas temperature is carried by it)
        const surfaceScalarField& massFlux() const noexcept
        {
            return massFlux_;
        }

        //- Superficial gas flux relative to the grains per cell [m/s]
        const volVectorField& qGas() const noexcept { return qGas_; }
        const gasDragModel& dragModel() const noexcept { return *dragModel_; }

        //- Gas mass in the domain, all processors [kg]
        scalar gasMass() const;

        //- Gas mass that came in through the boundaries so far [kg]
        scalar netIn() const noexcept { return netIn_; }

        //- Two-way coupled with the solid (WP F.3)
        bool coupled() const noexcept { return coupled_; }

        //- WP F.8g: the pore pressure only lowers the friction
        bool frictionOnly() const noexcept { return frictionOnly_; }

        //- WP F.8i: the zero pressure on the material's surface (route B
        //  with frictionOnly and a free surface)
        bool drainAtSurface() const noexcept
        {
            return frictionOnly_ && kappa_ > 0 && coupled_ && freeSurface_;
        }

        //- Face pressures for the solid (WP F.4): linear between cells,
        //  the void pressure (0 gauge) on faces between material and void
        const surfaceScalarField& pFace() const noexcept { return pf_; }

        //- Void cell flag field (WP F.4)
        const volScalarField& voidCells() const noexcept { return void_; }

        //- Gas mass the void (atmosphere) gave so far [kg] (part of netIn)
        scalar voidIn() const noexcept { return voidIn_; }

        //- Gas flow in the void (WP F.6, model gasFlow)
        bool flow() const noexcept { return flow_; }

        //- Gas velocity (WP F.6; gasFlow only)
        const volVectorField& Ug() const noexcept { return *UgPtr_; }

        //- Pore gas inertia on (WP F.7)
        bool poreInertia() const noexcept { return bool(inertialDrag_); }

        //- Force of the pore gas inertia on the grains per cell [N] (WP
        //  F.7; null without pore inertia)
        const volVectorField* gasForcePtr() const noexcept
        {
            return gasForce_.get();
        }


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

        //- Face pressures pFace() from p (WP F.4)
        void correctFacePressure();

        //- Restart state (WP F.9, decision 0141) to <time>/uniform/poreGas
        //  (each processor): the gas totals, the last step's length and
        //  pressure (first guess), the support pressures of frictionOnly
        //  and the lift check's clock; call at write times
        void writeState() const;

        //- Read it back from the start time if present (constructor)
        void readState();

        //- Parallel consistency check (WP F.9, poreGasParallelCheck.C):
        //  the largest mismatch of the processor-boundary values of the
        //  gas and solid fields against their neighbours', relative to
        //  each field's largest value; printed, and kept in
        //  parallelWorst_
        void checkParallel(const word& where);

        //- Range-of-validity check (decision 0135, poreGasLift.C): the
        //  share of the grain mass whose weight the gas carries; warns
        //  once it passes liftWarn (the bed lifts or expands, which the
        //  model does not represent), prints it at write times
        void checkLift();

        //- Carried share of the grain mass at the last step (decision 0135)
        scalar liftShare() const noexcept { return liftShare_; }
};


} // End namespace Foam

#endif

// ************************************************************************* //
