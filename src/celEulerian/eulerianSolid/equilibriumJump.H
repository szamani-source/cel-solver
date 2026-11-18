/*---------------------------------------------------------------------------*\
    Hydrostatic stress jump used by the Lagrange step (decisions 0030, 0039)
    and by the immersed walls (WP 4.3a): one definition for both.
\*---------------------------------------------------------------------------*/

#ifndef Foam_eulerianEquilibriumJump_H
#define Foam_eulerianEquilibriumJump_H

#include "symmTensor.H"

namespace Foam
{
namespace eulerianDetail
{

//- Hydrostatic stress jump over the vector d in a cell of density rho
//  under gravity g: -rho (g . d) gHat gHat for the stress, zero for other
//  types (decision 0030). Gravity not along an axis (an inclined layer):
//  the jump of a layer uniform across the dominant axis e of g,
//  -rho (d . e) (g e + e g - (g . e) e e) (decision 0039)
template<class Type>
inline Type equilibriumJump(const scalar, const vector&, const vector&)
{
    return Zero;
}

template<>
inline symmTensor equilibriumJump<symmTensor>
(
    const scalar rho,
    const vector& d,
    const vector& g
)
{
    const scalar gg = magSqr(g);
    if (gg <= 0)
    {
        return Zero;
    }
    // Dominant axis e of g
    const vector ga(cmptMag(g));
    const direction a =
        ga.x() >= ga.y() && ga.x() >= ga.z() ? 0 : (ga.y() >= ga.z() ? 1 : 2);
    if (sqr(g[a]) == gg)
    {
        // g along an axis: hydrostatic column (decision 0030)
        return (-rho*(g & d)/gg)*sqr(g);
    }
    // Tilted g (WP 3.6b, decision 0039): a layer varying along e only,
    // div sigma + rho g = 0, so sigma & e jumps by -rho g (d . e); the
    // normal stresses across e do not jump
    vector e(Zero);
    e[a] = 1;
    return (-rho*(d & e))*(twoSymm(g*e) - (g & e)*sqr(e));
}

} // End namespace eulerianDetail
} // End namespace Foam

#endif
