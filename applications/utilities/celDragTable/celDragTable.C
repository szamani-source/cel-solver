/*---------------------------------------------------------------------------*\
Application
    celDragTable

Description
    Prints the momentum exchange coefficient K [kg/(m^3 s)] of a gas-grain
    drag law (gasDragModel, library celFluid) for given gas fractions and
    slip speeds, 17 digits, one "e u K" line each. For checking the laws
    against tools/gas_reference.py (WP F.7, tests/test_drag_laws.sh).
    No case needed.

Usage
    celDragTable <law> <d> <mu> <rhoG>   then "e u" pairs on stdin

\*---------------------------------------------------------------------------*/

#include "gasDragModel.H"
#include "dictionary.H"
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace Foam;

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cerr << "usage: celDragTable <law> <d> <mu> <rhoG> < pairs\n";
        return 1;
    }
    const word law(argv[1]);
    const scalar d = std::strtod(argv[2], nullptr);
    const scalar mu = std::strtod(argv[3], nullptr);
    const scalar rhoG = std::strtod(argv[4], nullptr);

    const dictionary dict;
    autoPtr<gasDragModel> drag(gasDragModel::New(law, dict));

    double e, u;
    while (std::cin >> e >> u)
    {
        std::printf
        (
            "%.17g %.17g %.17g\n",
            e, u, drag->K(e, u, d, mu, rhoG)
        );
    }
    return 0;
}


// ************************************************************************* //
