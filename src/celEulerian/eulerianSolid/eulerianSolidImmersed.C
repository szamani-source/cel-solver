/*---------------------------------------------------------------------------*\
    eulerianSolid immersed walls (WP 4.3a, decision 0047): rigid STL bodies
    cut the fixed mesh (geometry from bodyCells, WP 4.2, decision 0043) and
    act on the material as cut-cell walls (decision 0042):

    - the material fills the open volume of each cell and crosses only the
      open part of each face;
    - a cell whose open faces do not close has one more face, the wall, with
      area vector minus the sum of the open face areas times their outward
      normals (so a uniform stress gives no force), solved in the Lagrange
      step as a one-sided Riemann problem against the body (no tension,
      Coulomb cap with the body's friction);
    - covered cells hold ghost values for the reconstructions;
    - small cut cells (open fraction below smallOpen) are mixed with their
      neighbour after each remap.

    This file: setting the geometry, ghost values, mixing, reports. The wall
    faces themselves are in eulerianSolidLagrange.C and the open volumes in
    eulerianSolidRemap.C.
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "equilibriumJump.H"
#include "mathematicalConstants.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * Local Helpers * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedDetail
{

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> work
(
    const fvMesh& mesh,
    const word& name,
    const Field<Type>& values
)
{
    auto t = GeometricField<Type, fvPatchField, volMesh>::New
    (
        name,
        IOobject::NO_REGISTER,
        mesh,
        dimensioned<Type>(dimless, Zero),
        fvPatchFieldBase::calculatedType()
    );
    t.ref().primitiveFieldRef() = values;
    t.ref().correctBoundaryConditions();
    return t;
}

} // End namespace immersedDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eulerianSolid::setImmersedWalls
(
    const volScalarField& alphaBody,
    const surfaceScalarField& faceOpen,
    const volVectorField& wallVelocity,
    const labelList& cellBody,
    const wordList& bodyNames,
    const scalarList& frictionDeg
)
{
    using namespace immersedDetail;

    if (bodyNames.size() != frictionDeg.size())
    {
        FatalErrorInFunction
            << "one friction angle per body" << exit(FatalError);
    }

    const fvMesh& mesh = mesh_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& magSf = mesh.magSf().primitiveField();
    const label nCells = mesh.nCells();

    immersed_ = bodyNames.size() > 0;
    bodyNames_ = bodyNames;
    bodyMu_.setSize(bodyNames.size());
    forAll(bodyMu_, b)
    {
        bodyMu_[b] = Foam::tan(frictionDeg[b]*constant::mathematical::pi/180);
    }
    if (bodyForce_.size() != bodyNames.size())
    {
        bodyForce_.setSize(bodyNames.size(), Zero);
        bodyMassIn_.setSize(bodyNames.size(), 0);
    }

    // Cell open fractions; covered cells hold no material
    cellOpen_.setSize(nCells);
    forAll(cellOpen_, celli)
    {
        const scalar o = 1 - alphaBody[celli];
        cellOpen_[celli] = o < coveredOpen ? 0 : min(o, scalar(1));
    }
    tmp<volScalarField> topen(work<scalar>(mesh, "cellOpen", cellOpen_));

    // Face open fractions: closed where either side is covered
    faceOpen_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "faceOpenEff", mesh.time().timeName(), mesh,
                IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
            ),
            faceOpen
        )
    );
    surfaceScalarField& fo = faceOpen_.ref();
    forAll(own, facei)
    {
        fo[facei] =
            (cellOpen_[own[facei]] > 0 && cellOpen_[nei[facei]] > 0)
          ? min(max(fo[facei], scalar(0)), scalar(1))
          : 0;
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        scalarField& pf = fo.boundaryFieldRef()[patchi];
        const scalarField on
        (
            p.coupled()
          ? topen().boundaryField()[patchi].patchNeighbourField()
          : tmp<scalarField>::New(p.size(), scalar(1))
        );
        forAll(fc, i)
        {
            pf[i] =
                (cellOpen_[fc[i]] > 0 && on[i] > 0)
              ? min(max(pf[i], scalar(0)), scalar(1))
              : 0;
        }
    }

    // Wall faces: minus the open face areas times outward normals
    wallA_.setSize(nCells);
    wallA_ = Zero;
    scalarField sumMag(nCells, Zero);
    forAll(own, facei)
    {
        const vector s(fo[facei]*Sf[facei]);
        wallA_[own[facei]] -= s;
        wallA_[nei[facei]] += s;
        sumMag[own[facei]] += magSf[facei];
        sumMag[nei[facei]] += magSf[facei];
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pS = p.Sf();
        const scalarField& pf = fo.boundaryField()[patchi];
        forAll(fc, i)
        {
            wallA_[fc[i]] -= pf[i]*pS[i];
            sumMag[fc[i]] += mag(pS[i]);
        }
    }

    // Directions without cells (2D) carry no wall
    vector mask(Zero);
    const Vector<label>& geomD = mesh.geometricD();
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        mask[d] = geomD[d] < 0 ? 0 : 1;
    }

    wallBody_.setSize(nCells);
    wallBody_ = -1;
    label nWall = 0, nCovered = 0, nOrphan = 0;
    forAll(wallA_, celli)
    {
        wallA_[celli] = cmptMultiply(mask, wallA_[celli]);
        if (cellOpen_[celli] <= 0)
        {
            wallA_[celli] = Zero;
            ++nCovered;
            continue;
        }
        if (mag(wallA_[celli]) <= 1e-8*sumMag[celli])
        {
            wallA_[celli] = Zero;
            continue;
        }
        ++nWall;
        label b = cellBody[celli];
        if (b < 0)
        {
            // Open cell whose face lies on a body: the neighbour's body
            scalar best = 0;
            for (const label facei : mesh.cells()[celli])
            {
                if (facei >= mesh.nInternalFaces()) continue;
                const label other =
                    own[facei] == celli ? nei[facei] : own[facei];
                const scalar closed = (1 - fo[facei])*magSf[facei];
                if (closed > best && cellBody[other] >= 0)
                {
                    best = closed;
                    b = cellBody[other];
                }
            }
        }
        if (b < 0)
        {
            ++nOrphan;
            b = 0;
        }
        wallBody_[celli] = b;
    }

    wallU_ = wallVelocity.primitiveField();

    if (freeSurface_)
    {
        // WP 4.3c (decision 0049): alpha is the material fraction of the
        // open volume. Open length per axis for the remap's slab rule: open
        // volume over the larger open area of the two faces across the axis
        wallPhi_.setSize(nCells);
        wallPhi_ = 0;
        vectorField Amax(nCells, Zero);
        auto addA = [&](const label c, const vector& S, const scalar f)
        {
            const vector a(f*cmptMag(S));
            Amax[c] = max(Amax[c], a);
        };
        forAll(own, facei)
        {
            addA(own[facei], Sf[facei], fo[facei]);
            addA(nei[facei], Sf[facei], fo[facei]);
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const labelUList& fc = p.faceCells();
            const vectorField& pS = p.Sf();
            const scalarField& pf = fo.boundaryField()[patchi];
            forAll(fc, i)
            {
                addA(fc[i], pS[i], pf[i]);
            }
        }
        const scalarField& V = mesh.V().field();
        const Vector<label>& geomD = mesh.geometricD();
        openLen_ = hCell_;
        forAll(openLen_, celli)
        {
            const scalar o = cellOpen_[celli];
            if (o >= 1 || o <= 0)
            {
                continue;
            }
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                if (geomD[d] < 0)
                {
                    continue;   // no cells along d (2D): the full length
                }
                openLen_[celli][d] =
                    Amax[celli][d] > 0
                  ? o*V[celli]/Amax[celli][d]
                  : o*hCell_[celli][d];
            }
        }

        // No material in covered cells (ghost values follow)
        scalarField& a = alpha_.primitiveFieldRef();
        scalarField& r = rho_.primitiveFieldRef();
        forAll(a, celli)
        {
            if (cellOpen_[celli] <= 0)
            {
                a[celli] = 0;
                r[celli] = 0;
            }
        }
        alphaL_ = a;
    }

    // Mixing partner of each small cut cell: the face neighbour on this
    // processor, not small itself, across the largest open area; if every
    // neighbour is small (a corner), step to the small neighbour across the
    // largest open area and look again (up to 3 steps)
    mixTo_.setSize(nCells);
    mixTo_ = -1;
    nUnmixed_ = 0;
    label nSmall = 0;
    forAll(cellOpen_, celli)
    {
        if (cellOpen_[celli] <= 0 || cellOpen_[celli] >= smallOpen)
        {
            continue;
        }
        ++nSmall;
        label cur = celli, prev = -1;
        for (label hop = 0; hop < 4 && mixTo_[celli] < 0; ++hop)
        {
            scalar bestBig = 0, bestSmall = 0;
            label big = -1, small = -1;
            for (const label facei : mesh.cells()[cur])
            {
                if (facei >= mesh.nInternalFaces()) continue;
                const label other =
                    own[facei] == cur ? nei[facei] : own[facei];
                const scalar a = fo[facei]*magSf[facei];
                if (a <= 0 || other == prev || other == celli) continue;
                if (cellOpen_[other] >= smallOpen)
                {
                    if (a > bestBig) { bestBig = a; big = other; }
                }
                else if (a > bestSmall)
                {
                    bestSmall = a;
                    small = other;
                }
            }
            if (big >= 0)
            {
                mixTo_[celli] = big;
            }
            else if (small >= 0)
            {
                prev = cur;
                cur = small;
            }
            else
            {
                break;
            }
        }
        if (mixTo_[celli] < 0)
        {
            ++nUnmixed_;
        }
    }

    Info<< "Immersed walls: " << returnReduce(nCovered, sumOp<label>())
        << " covered cells, " << returnReduce(nWall, sumOp<label>())
        << " wall cells, " << returnReduce(nSmall, sumOp<label>())
        << " small cut cells (" << nUnmixed() << " unmixed)";
    reduce(nOrphan, sumOp<label>());
    if (nOrphan)
    {
        Info<< ", " << nOrphan << " wall cells without a body (body 0)";
    }
    Info<< nl;
    forAll(bodyNames_, b)
    {
        Info<< "    body " << bodyNames_[b] << ": friction "
            << frictionDeg[b] << " deg" << nl;
    }
    if (nUnmixed() > 0)
    {
        WarningInFunction
            << nUnmixed() << " small cut cells have no neighbour to mix with "
            << "on their processor (mixing across processors: not done)"
            << endl;
    }

    fillCovered();
}


Foam::tmp<Foam::scalarField> Foam::eulerianSolid::cellOpen() const
{
    if (!immersed_)
    {
        return tmp<scalarField>::New(mesh_.nCells(), scalar(1));
    }
    return tmp<scalarField>::New(cellOpen_);
}


Foam::tmp<Foam::scalarField> Foam::eulerianSolid::openVolume() const
{
    if (!immersed_)
    {
        return tmp<scalarField>::New(mesh_.V().field());
    }
    return tmp<scalarField>::New(mesh_.V().field()*cellOpen_);
}


void Foam::eulerianSolid::fillCovered()
{
    using namespace immersedDetail;
    using eulerianDetail::equilibriumJump;

    const fvMesh& mesh = mesh_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& C = mesh.C().primitiveField();
    const label nCells = mesh.nCells();
    const vector g(gravityAt(mesh.time().value()));

    scalarField& rho = rho_.primitiveFieldRef();
    vectorField& U = U_.primitiveFieldRef();
    symmTensorField& sigma = sigma_.primitiveFieldRef();
    symmTensorField& epsP = epsP_.primitiveFieldRef();
    scalarField& epsPEq = epsPEq_.primitiveFieldRef();

    const bool fs = freeSurface_;
    scalarField& alpha = alpha_.primitiveFieldRef();

    scalarField n(nCells, Zero), r(nCells, Zero), q(nCells, Zero);
    scalarField al(nCells, Zero);
    vectorField u(nCells, Zero);
    symmTensorField s(nCells, Zero), e(nCells, Zero);

    // Ghost of covered cell c from open neighbour k: the neighbour's state
    // with its hydrostatic jump over C_c - C_k (a mirror ghost at a mesh
    // wall, decision 0030)
    auto add = [&]
    (
        const label c,
        const scalar rk, const vector& uk, const symmTensor& sk,
        const symmTensor& ek, const scalar qk, const vector& d,
        const scalar ak
    )
    {
        n[c] += 1;
        al[c] += ak;
        r[c] += rk;
        u[c] += uk;
        s[c] += sk + equilibriumJump<symmTensor>(rk, d, g);
        e[c] += ek;
        q[c] += qk;
    };

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const bool oP = cellOpen_[P] > 0, oN = cellOpen_[N] > 0;
        if (oP == oN) continue;
        if (oP)
        {
            add
            (
                N, rho[P], U[P], sigma[P], epsP[P], epsPEq[P], C[N] - C[P],
                fs ? alpha[P] : 1
            );
        }
        else
        {
            add
            (
                P, rho[N], U[N], sigma[N], epsP[N], epsPEq[N], C[P] - C[N],
                fs ? alpha[N] : 1
            );
        }
    }

    // Neighbours across processor (and cyclic) patches
    {
        tmp<volScalarField> to(work<scalar>(mesh, "o", cellOpen_));
        tmp<volScalarField> tr(work<scalar>(mesh, "r", rho));
        tmp<volVectorField> tu(work<vector>(mesh, "u", U));
        tmp<volSymmTensorField> ts(work<symmTensor>(mesh, "s", sigma));
        tmp<volSymmTensorField> te(work<symmTensor>(mesh, "e", epsP));
        tmp<volScalarField> tq(work<scalar>(mesh, "q", epsPEq));
        tmp<volScalarField> ta(work<scalar>(mesh, "a", alpha));
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!p.coupled()) continue;
            const labelUList& fc = p.faceCells();
            const vectorField d(p.delta());
            const scalarField on(to().boundaryField()[patchi].patchNeighbourField());
            const scalarField rn(tr().boundaryField()[patchi].patchNeighbourField());
            const vectorField un(tu().boundaryField()[patchi].patchNeighbourField());
            const symmTensorField sn(ts().boundaryField()[patchi].patchNeighbourField());
            const symmTensorField en(te().boundaryField()[patchi].patchNeighbourField());
            const scalarField qn(tq().boundaryField()[patchi].patchNeighbourField());
            const scalarField an(ta().boundaryField()[patchi].patchNeighbourField());
            forAll(fc, i)
            {
                if (cellOpen_[fc[i]] <= 0 && on[i] > 0)
                {
                    add
                    (
                        fc[i], rn[i], un[i], sn[i], en[i], qn[i], -d[i],
                        fs ? an[i] : 1
                    );
                }
            }
        }
    }

    forAll(cellOpen_, celli)
    {
        if (cellOpen_[celli] > 0) continue;
        if (n[celli] > 0)
        {
            rho[celli] = r[celli]/n[celli];
            U[celli] = u[celli]/n[celli];
            sigma[celli] = s[celli]/n[celli];
            epsP[celli] = e[celli]/n[celli];
            epsPEq[celli] = q[celli]/n[celli];
            if (fs)
            {
                alpha[celli] = al[celli]/n[celli];
            }
        }
        else if (fs)
        {
            // Away from the material: empty, at rest
            rho[celli] = 0;
            alpha[celli] = 0;
            U[celli] = wallU_[celli];
            sigma[celli] = Zero;
            epsP[celli] = Zero;
            epsPEq[celli] = 0;
        }
        else
        {
            rho[celli] = rho0_;
            U[celli] = wallU_[celli];
            sigma[celli] = Zero;
            epsP[celli] = Zero;
            epsPEq[celli] = 0;
        }
    }
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();
    epsP_.correctBoundaryConditions();
    epsPEq_.correctBoundaryConditions();
    if (fs)
    {
        alpha_.correctBoundaryConditions();
    }
}


void Foam::eulerianSolid::mixSmallCells()
{
    using eulerianDetail::equilibriumJump;

    const fvMesh& mesh = mesh_;
    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V().field();
    const label nCells = mesh.nCells();
    const vector g(gravityAt(mesh.time().value()));

    scalarField& rho = rho_.primitiveFieldRef();
    vectorField& U = U_.primitiveFieldRef();
    symmTensorField& sigma = sigma_.primitiveFieldRef();
    symmTensorField& epsP = epsP_.primitiveFieldRef();
    scalarField& epsPEq = epsPEq_.primitiveFieldRef();

    // Group of each target: the target and the small cells mixed into it.
    // Stress is mixed less its hydrostatic part relative to the target's
    // centre, so a column at rest stays at rest.
    const bool fs = freeSurface_;
    scalarField& alpha = alpha_.primitiveFieldRef();
    scalarField gm(nCells, Zero), gV(nCells, Zero), gq(nCells, Zero);
    scalarField ga(nCells, Zero);       // material volume (freeSurface)
    vectorField gp(nCells, Zero);
    symmTensorField gs(nCells, Zero), ge(nCells, Zero);
    boolList isTarget(nCells, false);

    auto add = [&](const label k, const label t)
    {
        const scalar mk = rho[k]*V[k]*cellOpen_[k];
        gm[t] += mk;
        gV[t] += V[k]*cellOpen_[k];
        if (fs)
        {
            ga[t] += alpha[k]*V[k]*cellOpen_[k];
        }
        gp[t] += mk*U[k];
        gs[t] += mk*(sigma[k] - equilibriumJump<symmTensor>(rho[k], C[k] - C[t], g));
        ge[t] += mk*epsP[k];
        gq[t] += mk*epsPEq[k];
    };

    forAll(mixTo_, celli)
    {
        const label t = mixTo_[celli];
        if (t < 0) continue;
        if (!isTarget[t])
        {
            isTarget[t] = true;
            add(t, t);
        }
        add(celli, t);
    }

    auto set = [&](const label k, const label t)
    {
        if (fs && (gm[t] <= 0 || ga[t] <= 0))
        {
            // Empty group (free surface, WP 4.3c)
            rho[k] = 0;
            alpha[k] = 0;
            U[k] = Zero;
            sigma[k] = Zero;
            epsP[k] = Zero;
            epsPEq[k] = 0;
            return;
        }
        const scalar r = gm[t]/gV[t];
        rho[k] = r;
        U[k] = gp[t]/gm[t];
        sigma[k] =
            gs[t]/gm[t] + equilibriumJump<symmTensor>(r, C[k] - C[t], g);
        epsP[k] = ge[t]/gm[t];
        epsPEq[k] = gq[t]/gm[t];
        if (fs)
        {
            const scalar a = ga[t]/gV[t];
            // more material than open volume: clipped (its mass stays,
            // compressed)
            alpha[k] = a >= 1 - 1e-12 ? 1 : max(a, scalar(0));

        }
    };
    forAll(mixTo_, celli)
    {
        const label t = mixTo_[celli];
        if (t >= 0 && (fs || gm[t] > 0))
        {
            set(celli, t);
        }
    }
    forAll(isTarget, t)
    {
        if (isTarget[t] && (fs || gm[t] > 0))
        {
            set(t, t);
        }
    }
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();
    epsP_.correctBoundaryConditions();
    epsPEq_.correctBoundaryConditions();
    if (fs)
    {
        alpha_.correctBoundaryConditions();
    }
}


Foam::vector Foam::eulerianSolid::bodyForce(const label b) const
{
    return returnReduce(bodyForce_[b], sumOp<vector>());
}


Foam::scalar Foam::eulerianSolid::bodyMassIn(const label b) const
{
    return returnReduce(bodyMassIn_[b], sumOp<scalar>());
}


Foam::label Foam::eulerianSolid::nSmallCells() const
{
    label n = 0;
    if (immersed_)
    {
        forAll(cellOpen_, celli)
        {
            if (cellOpen_[celli] > 0 && cellOpen_[celli] < smallOpen) ++n;
        }
    }
    return returnReduce(n, sumOp<label>());
}


Foam::label Foam::eulerianSolid::nUnmixed() const
{
    return returnReduce(nUnmixed_, sumOp<label>());
}


// ************************************************************************* //
