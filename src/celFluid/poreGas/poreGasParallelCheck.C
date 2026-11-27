/*---------------------------------------------------------------------------*\
    celFoam: explicit Coupled Eulerian-Lagrangian solver on OpenFOAM
    Licence: GPL-3.0-or-later

    Parallel consistency check (WP F.9, debug; gasProperties parallelCheck).

    Two kinds of error make a parallel run differ from the serial one by
    more than round-off, and both show on the processor boundaries:

    - a stale halo: a cell field changed after its last
      correctBoundaryConditions, so the neighbour processor still sees the
      old values (patchNeighbourField, interpolation, gradients);
    - a one-sided face value: the two processors compute a face's value
      (a flux, a face pressure, a coefficient) from different data, so the
      flux leaving one side is not the flux entering the other.

    The check refreshes a copy of each cell field and compares its
    processor-boundary values with the stored ones, and swaps the face
    fields across the coupled patches and compares the two sides (fluxes
    with opposite signs, face values equal). Each mismatch is relative to
    the field's largest |value| over all processors.
\*---------------------------------------------------------------------------*/

#include "poreGas.H"
#include "syncTools.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

// Largest |value| over all processors (internal and boundary)
template<class Type>
scalar fieldScale(const GeometricField<Type, fvPatchField, volMesh>& f)
{
    return max(gMax(mag(f.primitiveField())()), VSMALL);
}

template<class Type>
scalar fieldScale(const GeometricField<Type, fvsPatchField, surfaceMesh>& f)
{
    // One reduction: processors have different numbers of patches
    scalar s = max(mag(f.primitiveField())());
    forAll(f.boundaryField(), patchi)
    {
        if (f.boundaryField()[patchi].size())
        {
            s = max(s, max(mag(f.boundaryField()[patchi])()));
        }
    }
    return max(returnReduce(s, maxOp<scalar>()), VSMALL);
}


// Stale halo: stored processor-boundary values against fresh ones
template<class Type>
scalar staleHalo(const GeometricField<Type, fvPatchField, volMesh>& f)
{
    GeometricField<Type, fvPatchField, volMesh> c(f.name() + "ParCheck", f);
    c.correctBoundaryConditions();
    scalar d = 0;
    forAll(f.boundaryField(), patchi)
    {
        if (!f.boundaryField()[patchi].coupled())
        {
            continue;
        }
        const Field<Type> a(f.boundaryField()[patchi].patchNeighbourField());
        const Field<Type> b(c.boundaryField()[patchi].patchNeighbourField());
        if (a.size())
        {
            d = max(d, max(mag(a - b)()));
        }
    }
    return returnReduce(d, maxOp<scalar>())/fieldScale(f);
}


// One-sided face values: this side against the other (sign -1 for fluxes)
template<class Type>
scalar faceMismatch
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const scalar sign
)
{
    const fvMesh& mesh = f.mesh();
    const label nInt = mesh.nInternalFaces();
    List<Type> bf(mesh.nBoundaryFaces(), Zero);
    forAll(f.boundaryField(), patchi)
    {
        const label start = mesh.boundaryMesh()[patchi].start() - nInt;
        const Field<Type>& pf = f.boundaryField()[patchi];
        forAll(pf, i)
        {
            bf[start + i] = pf[i];
        }
    }
    List<Type> nb(bf);
    syncTools::swapBoundaryFaceList(mesh, nb);
    scalar d = 0;
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        if (!pp.coupled())
        {
            continue;
        }
        const label start = pp.start() - nInt;
        forAll(pp, i)
        {
            d = max(d, mag(bf[start + i] - sign*nb[start + i]));
        }
    }
    return returnReduce(d, maxOp<scalar>())/fieldScale(f);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::poreGas::checkParallel(const word& where)
{
    if
    (
        !parallelCheck_ || !Pstream::parRun()
     || solid_.mesh().time().timeIndex() % parallelCheckEvery_ != 0
    )
    {
        return;
    }

    List<Tuple2<word, scalar>> r;
    auto add = [&](const word& n, const scalar v)
    {
        r.append(Tuple2<word, scalar>(n, v));
    };

    // Cell fields: stale halos
    add("p", staleHalo(*pPtr_));
    add("eM", staleHalo(eM_));
    add("k", staleHalo(k_));
    add("phiS", staleHalo(phiS_));
    add("phiG", staleHalo(phiG_));
    add("void", staleHalo(void_));
    add("qGas", staleHalo(qGas_));
    if (UgPtr_) add("Ug", staleHalo(*UgPtr_));
    if (rAUPtr_) add("rAUg", staleHalo(*rAUPtr_));
    if (HbyAPtr_) add("HbyAg", staleHalo(*HbyAPtr_));
    if (gasForce_) add("gasForce", staleHalo(*gasForce_));
    add("alpha", staleHalo(solid_.alpha()));
    add("rho", staleHalo(solid_.rho()));
    add("U", staleHalo(solid_.U()));
    add("sigma", staleHalo(solid_.sigma()));

    // Face fields: the two sides of each coupled face
    add("massFlux", faceMismatch(massFlux_, -1));
    add("pFace", faceMismatch(pf_, 1));
    if (rAUfPtr_) add("rAUgf", faceMismatch(*rAUfPtr_, 1));
    if (mHbyAPtr_) add("mHbyA", faceMismatch(*mHbyAPtr_, -1));
    if (phiOld_) add("phiUgOld", faceMismatch(*phiOld_, -1));
    if (qnOld_) add("qnOld", faceMismatch(*qnOld_, -1));
    if (voidFaceOld_) add("voidFaceOld", faceMismatch(*voidFaceOld_, 1));
    if (sInertia_) add("sInertia", faceMismatch(*sInertia_, -1));
    if (gInertia_) add("gInertia", faceMismatch(*gInertia_, 1));
    if (mInertia_) add("mInertia", faceMismatch(*mInertia_, -1));
    add("phiL", faceMismatch(solid_.phiL(), -1));
    add("remapMass", faceMismatch(solid_.remapMass(), -1));

    scalar worst = 0;
    word worstName("none");
    for (const auto& e : r)
    {
        if (e.second() > worst)
        {
            worst = e.second();
            worstName = e.first();
        }
    }
    parallelWorst_ = max(parallelWorst_, worst);

    const label oldPrecision = Info().precision(3);
    Info<< "Parallel check (" << where << "): largest " << worst << " ("
        << worstName << "), run " << parallelWorst_;
    if (worst > 1e-12)
    {
        Info<< "; over 1e-12:";
        for (const auto& e : r)
        {
            if (e.second() > 1e-12)
            {
                Info<< ' ' << e.first() << ' ' << e.second();
            }
        }
    }
    Info<< nl;
    Info().precision(oldPrecision);
}


// ************************************************************************* //
