/*---------------------------------------------------------------------------*\
    Nonlocal fluidity solve (WP N.1, decision 0111): k g - div(a_f grad g)
    = k g_loc for the stress sigma and material fraction a per cell. Used by
    eulerianSolid each step and by tests/nonlocal/fluidityCheck.
\*---------------------------------------------------------------------------*/

#ifndef nonlocalFluidity_H
#define nonlocalFluidity_H

#include "volFields.H"
#include "celNonlocal.H"

namespace Foam
{

void nonlocalFluidity
(
    const fvMesh& mesh,
    const cel::constitutive::Nonlocal& law,
    const volSymmTensorField& sigma,
    const scalarField& a,
    volScalarField& g
);

}

#endif
