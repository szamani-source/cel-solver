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
        conductivity constant;      // bed conductivity model (T.3 adds
                                    // zehnerSchlunder)
        k            0.25;          // bed conductivity [W/(m K)]

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

        const eulerianSolid& solid_;

        const fvMesh& mesh_;

        IOdictionary dict_;

        word model_;

        //- Specific heat of the grains [J/(kg K)]
        scalar c_;

        //- Bed conductivity model and value [W/(m K)]
        word conductivity_;
        scalar k_;

        //- Temperature of the grains [K]
        volScalarField T_;

        //- Heat that entered through the boundaries since the start [J]
        scalar heatIn_;

        //- Thermal energy at the start of this run [J]
        scalar energy0_;

        //- Sub-cycles used in the last step
        label nSub_;


    // Private Member Functions

        //- Material fraction per cell (alpha, or 1 without a free surface)
        tmp<scalarField> fill() const;

        //- The same as a field with processor halos
        tmp<volScalarField> fillField() const;

        //- One explicit conduction step of length dt
        void conduct(const scalar dt);


public:

    // Constructors

        heatTransfer(const eulerianSolid& solid, const IOobject& io);

        //- Null without constant/thermalProperties or with model none
        static autoPtr<heatTransfer> New(const eulerianSolid& solid);

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
