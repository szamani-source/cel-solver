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
    F.4). Nothing here acts on the solid: no pressure is solved (F.2).

    constant/gasProperties:
    \verbatim
        model   poreGas;        // none | poreGas (none: as if absent)
        rhoP    2650;           // grain (particle) density [kg/m^3]
        d       0.001;          // grain diameter [m]
        T       293.15;         // gas temperature [K], default 293.15
        p0      101325;         // reference (ambient) pressure [Pa], default
        drag    ergunWenYu;     // drag law; F.2 adds the run-time selector,
                                // F.7 the other laws (decision 0125 item 7)
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

        //- Drag law name (only ergunWenYu until F.2/F.7)
        word drag_;

        //- Gas density at p0 and T [kg/m^3], viscosity at T [Pa s]
        scalar rhoG_;
        scalar mu_;

        //- Solids and gas volume fractions, permeability of the filled part
        volScalarField phiS_;
        volScalarField phiG_;
        volScalarField k_;


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


    // Cycle

        //- Fractions and permeability from the solid's current rho and
        //  alpha; call after remap() (and after moving walls). The fields
        //  are written with the case's other fields.
        void correctFractions();
};


} // End namespace Foam

#endif

// ************************************************************************* //
