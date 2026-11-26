/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModels::diFelice

Description
    Di Felice (1994):

        K = 0.75 Cd rhoG (1-e) u/d e^(2 - beta)
          = 0.75 Cd Re mu (1-e)/d^2 e^(1 - beta),
        Cd = (0.63 + 4.8/sqrt(Re))^2   (Dallavalle),
        beta = 3.7 - 0.65 exp(-(1.5 - log10 Re)^2/2),
        Re = e rhoG u d/mu (superficial slip, Di Felice's definition)

    i.e. one sphere's drag at the superficial slip times e^-beta
    (El Geitani et al. 2022 eq. 30 write the same Cd e^(2-beta); they
    give Re without e). Cd Re = (0.63 sqrt(Re) + 4.8)^2 and beta -> 3.7
    as Re -> 0, so u = 0 is finite.

    A user option (WP F.7, decision 0125 item 7); every test uses the
    default ergunWenYu (Sia, 2026-11-24). Checked against
    tools/gas_reference.py by tests/test_drag_laws.sh.

    \verbatim
        drag    diFelice;        // constant/gasProperties
    \endverbatim

SourceFiles
    diFelice.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModels_diFelice_H
#define Foam_gasDragModels_diFelice_H

#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{

class diFelice
:
    public gasDragModel
{
public:

    //- Runtime type information
    TypeName("diFelice");


    // Constructors

        explicit diFelice(const dictionary& dict);


    //- Destructor
    virtual ~diFelice() = default;


    // Member Functions

        virtual scalar K
        (
            const scalar e,
            const scalar u,
            const scalar d,
            const scalar mu,
            const scalar rhoG
        ) const;
};


} // End namespace gasDragModels
} // End namespace Foam

#endif

// ************************************************************************* //
