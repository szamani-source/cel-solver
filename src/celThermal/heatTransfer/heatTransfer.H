/*---------------------------------------------------------------------------*\
Class
    Foam::heatTransfer

Description
    Temperature in the granular material (version 2 temperature work, T
    series, decision 0135). Library celThermal, called by the celFoam
    driver after the gas solve.

    WP T.1: conduction in the grains of a bed that may move but is not yet
    carried by the remap (T.2). Per cell the energy is C T with the heat
    capacity C = rho c V, rho the bulk density of the solid (alpha times
    the material density, phi_s rho_s). The explicit finite-volume update

        C_P (T_P' - T_P) = dt sum_f k_f |S_f| (n . grad T)_f

    uses the over-relaxed non-orthogonal split of OpenFOAM's corrected
    Laplacian:

        (n . grad T)_f = (T_N - T_P) delta_f + c_f . (grad T)_f,
        delta_f = 1/max(n . d, 0.05 |d|),  c_f = n - delta_f d,

    with grad T the Gauss gradient (linear face values), so that a linear
    temperature is reproduced exactly on non-orthogonal parallelogram cells.
    k_f = k min(a_P, a_N), a the material fraction (1 without a free
    surface): empty cells neither store nor conduct heat.

    Boundaries come from the patch types in 0/T:
      fixedValue (and any type that fixes the value): Dirichlet, with the
          non-orthogonal correction from the cell gradient;
      zeroGradient: adiabatic;
      fixedGradient: heat flux q = k gradient into the material (gradient =
          dT/dn on the outward normal);
      processor, cyclic: as interior faces.

    The step is sub-cycled if the solid step ever exceeds half the explicit
    conduction limit min_P C_P / sum_f k_f |S_f| delta_f (about 20 s on
    4 mm cells of sand; it does not happen in practice). The heat that
    crosses the boundaries is summed (heatIn) so that the energy balance can
    be checked to round-off.

    constant/thermalProperties:

        model        equilibrium;   // none | equilibrium (T.1);
                                    // twoTemperature from T.4
        c            830;           // specific heat of the grains [J/(kg K)]
        conductivity constant;      // constant | zehnerSchlunder (T.3)
        k            0.25;          // bed conductivity [W/(m K)] (constant)
        kParticle    3;             // zehnerSchlunder: grain material [W/(m K)]
        kGas         0.02569;       // zehnerSchlunder: gas [W/(m K)], air 20 C
        rhoP         2650;          // zehnerSchlunder: grain density [kg/m^3]
        walls                       // optional, overrides the 0/T patch type
        {
            hot  { h 50; Tinf 350; }   // convective (film + half cell)
            side { q 1000; }           // heat flux into the material [W/m^2]
        }

    WP T.3: zehnerSchlunder takes the packing of the filled part,
    phi_s = rho/(alpha rhoP), per cell (Zehner and Schluender 1970 as split
    by Kuipers et al. 1992; docs/spec/thermal-references.md); faces use the
    harmonic mean of the two cells. Convective walls put the film 1/h and
    the half cell 1/(k delta) in series (no non-orthogonal correction).
        limiter      superbee;      // remap limiter of T (T.2), as the
                                    // history fields (decision 0017)
        voidT        293.15;        // written in empty cells [K] (T.2)

    WP T.2: T moves with the material. The solid's remap carries it mass
    weighted (eulerianSolid::addPassenger), so sum m c T is conserved; an
    emptied cell holds no heat (voidT is only what is written there); open
    outlets carry out m c T of the cell, open inlets bring material at the
    inlet patch value of 0/T. Immersed STL walls are not yet supported (WP
    T.2b): the remap's small-cell mixing does not carry T.

    With no file, or model none, nothing is created and every v1 and v2
    result is unchanged bit for bit.

SourceFiles
    heatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef celThermal_heatTransfer_H
#define celThermal_heatTransfer_H

#include "eulerianSolid.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatTransfer Declaration
\*---------------------------------------------------------------------------*/

class heatTransfer
{
    // Private Data

        eulerianSolid& solid_;

        const fvMesh& mesh_;

        IOdictionary dict_;

        word model_;

        //- Specific heat of the grains [J/(kg K)]
        scalar c_;

        //- Bed conductivity model and value [W/(m K)]
        word conductivity_;
        scalar k_;

        //- Zehner-Schluender inputs (WP T.3): gas and particle
        //  conductivity [W/(m K)], particle density [kg/m^3]
        scalar kGas_;
        scalar kParticle_;
        scalar rhoP_;

        //- Temperature of the grains [K]
        volScalarField T_;

        //- Heat that entered through the boundaries since the start [J]
        scalar heatIn_;

        //- Thermal energy at the start of this run [J]
        scalar energy0_;

        //- Heat carried in by material through inlets less that carried
        //  out through outlets, since the start [J] (WP T.2)
        scalar carriedIn_;

        //- Temperature before the remap (outlet faces carry it out)
        scalarField TOld_;

        //- Remap limiter of T and the value written in empty cells
        word limiter_;
        scalar voidT_;

        //- Bed conductivity per cell [W/(m K)] (written as kBed with
        //  zehnerSchlunder)
        volScalarField kCell_;

        //- Wall heat transfer per patch (WP T.3)
        enum wallKind { NONE, CONVECTIVE, FLUX };
        List<wallKind> wallType_;
        scalarList wallH_;
        scalarList wallTinf_;
        scalarList wallQ_;

        //- Sub-cycles used in the last step
        label nSub_;


    // Private Member Functions

        //- Material fraction per cell (alpha, or 1 without a free surface)
        tmp<scalarField> fill() const;

        //- The same as a field with processor halos
        tmp<volScalarField> fillField() const;

        //- Bed conductivity per cell from the current packing
        void correctK();

        //- One explicit conduction step of length dt
        void conduct(const scalar dt);


public:

    // Constructors

        heatTransfer(eulerianSolid& solid, const IOobject& io);

        //- Null without constant/thermalProperties or with model none
        static autoPtr<heatTransfer> New(eulerianSolid& solid);

        heatTransfer(const heatTransfer&) = delete;
        void operator=(const heatTransfer&) = delete;


    // Member Functions

        const volScalarField& T() const noexcept { return T_; }

        scalar c() const noexcept { return c_; }

        //- Heat capacity per cell, rho c V [J/K]
        tmp<scalarField> capacity() const;

        //- Total thermal energy sum C T [J] (all processors)
        scalar energy() const;

        //- Heat that entered through the boundaries [J] (all processors)
        scalar heatIn() const noexcept { return heatIn_; }

        //- Explicit conduction limit [s] (all processors)
        scalar stableDeltaT() const;

        //- Call before the solid's remap: keeps T for the outlet faces
        void beforeRemap();

        //- Call after the solid's remap (which carries T, mass weighted):
        //  counts the heat carried through outlets and inlets and sets
        //  empty cells to voidT
        void afterRemap();

        //- Heat carried in through inlets less out through outlets [J]
        scalar carriedIn() const noexcept { return carriedIn_; }

        //- Advance the temperature over dt
        void solve(const scalar dt);

        //- Print energy, heat in and their balance
        void report() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
