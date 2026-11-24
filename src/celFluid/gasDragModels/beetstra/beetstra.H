/*---------------------------------------------------------------------------*\
Class
    Foam::gasDragModels::beetstra

Description
    Beetstra, van der Hoef and Kuipers (2007), lattice-Boltzmann fit of
    the normalised drag F = F_d/(3 pi mu d e u), with phi = 1 - e and
    Re = e rhoG u d/mu (superficial slip):

        F = 10 phi/e^2 + e^2 (1 + 1.5 sqrt(phi))
          + 0.413 Re/(24 e^2) (1/e + 3 e phi + 8.4 Re^-0.343)
           /(1 + 10^(3 phi) Re^-(1 + 4 phi)/2)

        K = 18 mu e (1-e) F/d^2

    (so a packed bed at rest gives Carman-Kozeny, 180 (1-e)^2 mu/(e^3 d^2)
    per superficial velocity, and e -> 1 one Stokes sphere). The inertial
    term vanishes at Re = 0.

    A user option (WP F.7, decision 0125 item 7); every test uses the
    default ergunWenYu (Sia, 2026-11-24). Checked against
    tools/gas_reference.py by tests/test_drag_laws.sh.

    \verbatim
        drag    beetstra;        // constant/gasProperties
    \endverbatim

SourceFiles
    beetstra.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gasDragModels_beetstra_H
#define Foam_gasDragModels_beetstra_H

#include "gasDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace gasDragModels
{

class beetstra
:
    public gasDragModel
{
public:

    //- Runtime type information
    TypeName("beetstra");


    // Constructors

        explicit beetstra(const dictionary& dict);


    //- Destructor
    virtual ~beetstra() = default;


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
