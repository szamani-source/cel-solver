/*---------------------------------------------------------------------------*\
    eulerianSolid::remap (WP 2.4): see eulerianSolid.H and
    docs/decisions/0014-lagrange-remap-contract.md (what remap receives) and
    docs/decisions/0017-remap-cpp.md (how it is done).

    One unsplit, face-based, second-order remap on any fvMesh:

    1. Cell "material velocity" uC: least-squares fit of the face fluxes,
       sum_f (uC.S_f - phiL_f)^2 -> min. Exact for a translation; in 1D it is
       the mean of the two face velocities, so the Lagrangian centroid
       C + dt*uC is the exact 1D centroid.
    2. Face displacement d_f = dt*(tangential part of the mean uC of the two
       cells + phiL/|S| along the normal). The swept region of face f has
       volume dt*phiL and centroid Cf + d_f/2 (exact for a translation).
    3. Per field component, a least-squares gradient g on the Lagrangian
       centroids, limited by one scalar factor beta per cell (Barth-Jespersen
       form): for each face, t_f = bound/(g.d_PN), with bound the largest
       (g.d > 0) or smallest (g.d < 0) neighbour difference max/min(0, q_N - q)
       over all the cell's neighbours, d_PN between Lagrangian centroids.
       MC: min(1, 2 t_min), superbee: max(min(2 t_min, t_max),
       min(t_min, 2 t_max)), minmod: min(1, t_min), upwind: 0. On a uniform
       1D mesh t is the one-sided over the central slope, t_l + t_r = 2, and
       these are exactly the 1D limiters of prototype/lr1d*.py. beta is also
       capped so the value at every swept centroid the cell donates to stays
       in the neighbour range (never active in 1D). Using the neighbour range
       rather than each face's own difference keeps round-off in one
       direction from switching off the gradient in another.
    4. Donor value at the swept centroid; mass flux dm = dt*phiL*rho_f;
       velocity, stress and history per unit mass fluxed with dm.
    5. New state written as increments of the Lagrangian state, so a cell
       with no flux through any face keeps every value bit for bit.

    Non-coupled boundary faces (wall, symmetry, empty, ...) carry no flux:
    their phiL is not used. Coupled faces (processor, cyclic) are handled
    like internal faces with neighbour data exchanged through temporary
    volFields, so both sides compute the same flux (to round-off).
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * Local helpers * * * * * * * * * * * * * * //

namespace Foam
{
namespace remapDetail
{

// Relative size below which g.d counts as zero (direction normal to the
// gradient): such faces do not constrain the limiter.
static const scalar perpTol = 1e-10;

//- Limiter types, see the file header
enum class limiterType { upwind, minmod, MC, superbee };

static limiterType limiterFromWord(const word& name, const dictionary& dict)
{
    if (name == "upwind") return limiterType::upwind;
    if (name == "minmod") return limiterType::minmod;
    if (name == "MC") return limiterType::MC;
    if (name == "superbee") return limiterType::superbee;

    FatalIOErrorInFunction(dict)
        << "Unknown remap limiter " << name
        << ". Valid: upwind minmod MC superbee" << exit(FatalIOError);
    return limiterType::upwind;
}


//- Limiter factor beta from the smallest and largest face ratio t (>= 0)
static inline scalar beta(limiterType lim, scalar tMin, scalar tMax)
{
    switch (lim)
    {
        case limiterType::upwind:
            return 0;
        case limiterType::minmod:
            return min(scalar(1), tMin);
        case limiterType::MC:
            return min(scalar(1), 2*tMin);
        case limiterType::superbee:
            return max(min(2*tMin, tMax), min(tMin, 2*tMax));
    }
    return 0;
}


//- Neighbour-cell values of a cell list across every coupled patch
//  (processor, cyclic). Empty fields for other patches.
template<class Type>
static List<Field<Type>> neighbourValues
(
    const fvMesh& mesh,
    const Field<Type>& cellValues
)
{
    GeometricField<Type, fvPatchField, volMesh> f
    (
        IOobject
        (
            "remapNeighbour",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensioned<Type>(dimless, Zero),
        fvPatchFieldBase::calculatedType()
    );
    f.primitiveFieldRef() = cellValues;
    f.correctBoundaryConditions();

    List<Field<Type>> nbr(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        if (mesh.boundary()[patchi].coupled())
        {
            nbr[patchi] = f.boundaryField()[patchi].patchNeighbourField();
        }
    }
    return nbr;
}


//- Geometry of one remap: everything that does not depend on the field
struct remapGeometry
{
    //- Lagrangian centroid differences P -> N, internal faces
    vectorField dInt;

    //- Same for each coupled patch (empty for others)
    List<vectorField> dPatch;

    //- Inverse least-squares matrix per cell
    symmTensorField DInv;

    //- Swept-region centroid minus the donor's Lagrangian centroid
    vectorField rInt;
    List<vectorField> rPatch;

    //- Directions without cells (2D): gradient component forced to 0
    vector emptyMask;
};


static remapGeometry makeGeometry
(
    const fvMesh& mesh,
    const surfaceScalarField& phi,
    const scalar dt
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& C = mesh.C().primitiveField();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& magSf = mesh.magSf().primitiveField();
    const vectorField& Cf = mesh.Cf().primitiveField();
    const scalarField& phiI = phi.primitiveField();
    const label nCells = mesh.nCells();
    const fvBoundaryMesh& bm = mesh.boundary();

    // Directions without cells (empty in 2D): regularise and mask
    const Vector<label>& geomD = mesh.geometricD();
    vector mask(Zero);
    symmTensor emptyId(Zero);
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        mask[cmpt] = (geomD[cmpt] < 0 ? 0 : 1);
    }
    if (geomD.x() < 0) emptyId.xx() = 1;
    if (geomD.y() < 0) emptyId.yy() = 1;
    if (geomD.z() < 0) emptyId.zz() = 1;

    // 1. Cell velocity from the face fluxes (least squares)
    symmTensorField M(nCells, emptyId);
    vectorField b(nCells, Zero);
    forAll(own, facei)
    {
        const symmTensor SS(sqr(Sf[facei]));
        M[own[facei]] += SS;
        M[nei[facei]] += SS;
        b[own[facei]] += Sf[facei]*phiI[facei];
        b[nei[facei]] += Sf[facei]*phiI[facei];  // (-S)(-phi) for the neighbour
    }
    forAll(bm, patchi)
    {
        const fvPatch& p = bm[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pSf = p.Sf();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        forAll(fc, i)
        {
            M[fc[i]] += sqr(pSf[i]);
            if (p.coupled())
            {
                b[fc[i]] += pSf[i]*pPhi[i];
            }
        }
    }
    vectorField uC(nCells);
    forAll(uC, celli)
    {
        uC[celli] = cmptMultiply(mask, inv(M[celli]) & b[celli]);
    }
    const List<vectorField> uNbr(neighbourValues(mesh, uC));

    remapGeometry geo;
    geo.emptyMask = mask;

    // 2. Lagrangian centroid differences and the least-squares matrices
    const vectorField deltaI(mesh.delta()().primitiveField());
    symmTensorField D(nCells, emptyId);

    geo.dInt.setSize(own.size());
    geo.rInt.setSize(own.size());
    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector d(deltaI[facei] + dt*(uC[N] - uC[P]));
        geo.dInt[facei] = d;
        const symmTensor wdd(sqr(d)/magSqr(d));
        D[P] += wdd;
        D[N] += wdd;

        // Swept centroid relative to the donor's Lagrangian centroid
        const vector n(Sf[facei]/magSf[facei]);
        const vector uMean(0.5*(uC[P] + uC[N]));
        const vector disp
        (
            dt*((uMean - (uMean & n)*n) + (phiI[facei]/magSf[facei])*n)
        );
        const vector xs(Cf[facei] + 0.5*disp);
        geo.rInt[facei] =
        (
            phiI[facei] > 0
          ? xs - (C[P] + dt*uC[P])
          : xs - (C[N] + dt*uC[N])
        );
    }

    geo.dPatch.setSize(bm.size());
    geo.rPatch.setSize(bm.size());
    forAll(bm, patchi)
    {
        const fvPatch& p = bm[patchi];
        const labelUList& fc = p.faceCells();

        if (p.coupled())
        {
            const vectorField pDelta(p.delta());
            const vectorField& pSf = p.Sf();
            const scalarField& pMagSf = p.magSf();
            const vectorField& pCf = p.Cf();
            const scalarField& pPhi = phi.boundaryField()[patchi];
            const vectorField& uN = uNbr[patchi];

            vectorField& dp = geo.dPatch[patchi];
            vectorField& rp = geo.rPatch[patchi];
            dp.setSize(fc.size());
            rp.setSize(fc.size());
            forAll(fc, i)
            {
                const label P = fc[i];
                const vector d(pDelta[i] + dt*(uN[i] - uC[P]));
                dp[i] = d;
                D[P] += sqr(d)/magSqr(d);

                const vector n(pSf[i]/pMagSf[i]);
                const vector uMean(0.5*(uC[P] + uN[i]));
                const vector disp
                (
                    dt*((uMean - (uMean & n)*n) + (pPhi[i]/pMagSf[i])*n)
                );
                const vector xs(pCf[i] + 0.5*disp);
                rp[i] =
                (
                    pPhi[i] > 0
                  ? xs - (C[P] + dt*uC[P])
                  : xs - (C[P] + pDelta[i] + dt*uN[i])
                );
            }
        }
        else
        {
            // Mirror ghost cell across the face (zero-gradient side):
            // regularises the fit next to walls and symmetry planes.
            const vectorField nf(p.nf());
            const vectorField& pCf = p.Cf();
            forAll(fc, i)
            {
                const vector d(2*(nf[i] & (pCf[i] - C[fc[i]]))*nf[i]);
                if (magSqr(d) > VSMALL)
                {
                    D[fc[i]] += sqr(d)/magSqr(d);
                }
            }
        }
    }

    geo.DInv.setSize(nCells);
    forAll(D, celli)
    {
        geo.DInv[celli] = inv(D[celli]);
    }

    return geo;
}


//- Limited gradient of one scalar component on the Lagrangian centroids
static vectorField limitedGradient
(
    const fvMesh& mesh,
    const remapGeometry& geo,
    const surfaceScalarField& phi,
    const scalarField& q,
    const List<scalarField>& qNbr,
    const limiterType lim
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const fvBoundaryMesh& bm = mesh.boundary();
    const label nCells = mesh.nCells();

    vectorField g(nCells, Zero);
    if (lim == limiterType::upwind)
    {
        return g;
    }

    // Least squares, weights 1/|d|^2 (ghost cells add nothing to b: dq = 0)
    vectorField b(nCells, Zero);
    forAll(own, facei)
    {
        const vector& d = geo.dInt[facei];
        const vector wdq(d*(q[nei[facei]] - q[own[facei]])/magSqr(d));
        b[own[facei]] += wdq;
        b[nei[facei]] += wdq;
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const vectorField& dp = geo.dPatch[patchi];
        const scalarField& qn = qNbr[patchi];
        forAll(fc, i)
        {
            b[fc[i]] += dp[i]*(qn[i] - q[fc[i]])/magSqr(dp[i]);
        }
    }
    forAll(g, celli)
    {
        g[celli] = cmptMultiply(geo.emptyMask, geo.DInv[celli] & b[celli]);
    }

    // Range of the neighbours: dMax = max(0, q_N - q), dMin = min(0, q_N - q)
    scalarField dMax(nCells, Zero);
    scalarField dMin(nCells, Zero);
    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar dq = q[N] - q[P];
        dMax[P] = max(dMax[P], dq);
        dMin[P] = min(dMin[P], dq);
        dMax[N] = max(dMax[N], -dq);
        dMin[N] = min(dMin[N], -dq);
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const scalarField& qn = qNbr[patchi];
        forAll(fc, i)
        {
            const scalar dq = qn[i] - q[fc[i]];
            dMax[fc[i]] = max(dMax[fc[i]], dq);
            dMin[fc[i]] = min(dMin[fc[i]], dq);
        }
    }

    // Allowed fraction r = (bound)/(g.x) of an extrapolation g.x, x the
    // offset from the cell's Lagrangian centroid
    auto ratio = [&](const label celli, const scalar e) -> scalar
    {
        return e > 0 ? dMax[celli]/e : dMin[celli]/e;
    };

    // Face ratios at the midpoints between Lagrangian centroids (in 1D on a
    // uniform mesh r = 2 t with t the one-sided over the central slope)
    scalarField rMin(nCells, GREAT);
    scalarField rMax(nCells, -GREAT);
    auto addFace = [&](const label celli, const vector& d)
    {
        const scalar gd = g[celli] & d;
        if (mag(gd) <= perpTol*mag(g[celli])*mag(d))
        {
            return;     // face (nearly) normal to the gradient
        }
        const scalar r = ratio(celli, 0.5*gd);
        rMin[celli] = min(rMin[celli], r);
        rMax[celli] = max(rMax[celli], r);
    };

    // Extra cap: the value at each swept centroid this cell donates to stays
    // within the neighbour range (weaker than the midpoint bound in 1D)
    scalarField cap(nCells, GREAT);
    auto addSwept = [&](const label celli, const vector& x)
    {
        const scalar e = g[celli] & x;
        if (mag(e) <= perpTol*mag(g[celli])*mag(x))
        {
            return;
        }
        cap[celli] = min(cap[celli], ratio(celli, e));
    };

    const scalarField& phiI = phi.primitiveField();
    forAll(own, facei)
    {
        addFace(own[facei], geo.dInt[facei]);
        addFace(nei[facei], -geo.dInt[facei]);
        if (phiI[facei] > 0)
        {
            addSwept(own[facei], geo.rInt[facei]);
        }
        else if (phiI[facei] < 0)
        {
            addSwept(nei[facei], geo.rInt[facei]);
        }
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const vectorField& dp = geo.dPatch[patchi];
        const scalarField& pPhi = phi.boundaryField()[patchi];
        forAll(fc, i)
        {
            addFace(fc[i], dp[i]);
            if (pPhi[i] > 0)
            {
                addSwept(fc[i], geo.rPatch[patchi][i]);
            }
        }
    }

    forAll(g, celli)
    {
        if (rMax[celli] < 0)
        {
            g[celli] = Zero;    // no face along the gradient
        }
        else
        {
            g[celli] *= min
            (
                beta(lim, 0.5*rMin[celli], 0.5*rMax[celli]),
                cap[celli]
            );
        }
    }

    return g;
}


//- Donor value of one component at the swept centroid of every face with
//  flux: internal faces and coupled patch faces
static void faceValues
(
    const fvMesh& mesh,
    const remapGeometry& geo,
    const surfaceScalarField& phi,
    const scalarField& q,
    const limiterType lim,
    scalarField& qfInt,
    List<scalarField>& qfPatch
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const fvBoundaryMesh& bm = mesh.boundary();
    const scalarField& phiI = phi.primitiveField();

    const List<scalarField> qNbr(neighbourValues(mesh, q));
    const vectorField g(limitedGradient(mesh, geo, phi, q, qNbr, lim));
    const List<vectorField> gNbr(neighbourValues(mesh, g));

    qfInt.setSize(own.size());
    forAll(own, facei)
    {
        const label donor = (phiI[facei] > 0 ? own[facei] : nei[facei]);
        qfInt[facei] = q[donor] + (g[donor] & geo.rInt[facei]);
    }

    qfPatch.setSize(bm.size());
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        const vectorField& rp = geo.rPatch[patchi];
        scalarField& qf = qfPatch[patchi];
        qf.setSize(fc.size());
        forAll(fc, i)
        {
            qf[i] =
            (
                pPhi[i] > 0
              ? q[fc[i]] + (g[fc[i]] & rp[i])
              : qNbr[patchi][i] + (gNbr[patchi][i] & rp[i])
            );
        }
    }
}


//- Remap one mass-weighted component: q <- (m q - sum dm q_f)/m_new
static void remapMassWeighted
(
    const fvMesh& mesh,
    const remapGeometry& geo,
    const surfaceScalarField& phi,
    const scalarField& dmInt,
    const List<scalarField>& dmPatch,
    const scalarField& mNew,
    const limiterType lim,
    scalarField& q
)
{
    scalarField qfInt;
    List<scalarField> qfPatch;
    faceValues(mesh, geo, phi, q, lim, qfInt, qfPatch);

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const fvBoundaryMesh& bm = mesh.boundary();

    // sum over the cell's faces of dm_out*(q_f - q), from the old q
    scalarField sumF(q.size(), Zero);
    forAll(own, facei)
    {
        const scalar dm = dmInt[facei];
        if (dm == 0) continue;
        sumF[own[facei]] += dm*(qfInt[facei] - q[own[facei]]);
        sumF[nei[facei]] -= dm*(qfInt[facei] - q[nei[facei]]);
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        forAll(fc, i)
        {
            const scalar dm = dmPatch[patchi][i];
            if (dm == 0) continue;
            sumF[fc[i]] += dm*(qfPatch[patchi][i] - q[fc[i]]);
        }
    }

    forAll(q, celli)
    {
        if (sumF[celli] != 0)
        {
            q[celli] -= sumF[celli]/mNew[celli];
        }
    }
}


//- Remap every component of a volField, mass-weighted
template<class Type>
static void remapField
(
    const fvMesh& mesh,
    const remapGeometry& geo,
    const surfaceScalarField& phi,
    const scalarField& dmInt,
    const List<scalarField>& dmPatch,
    const scalarField& mNew,
    const limiterType lim,
    GeometricField<Type, fvPatchField, volMesh>& f
)
{
    Field<Type>& fi = f.primitiveFieldRef();
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        scalarField q(fi.component(cmpt));
        remapMassWeighted(mesh, geo, phi, dmInt, dmPatch, mNew, lim, q);
        fi.replace(cmpt, q);
    }
    f.correctBoundaryConditions();
}

} // End namespace remapDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::HashTable<Foam::word> Foam::eulerianSolid::readRemapLimiters
(
    const dictionary& materialProperties
)
{
    // Defaults: decision 0013 item 5
    HashTable<word> lim;
    lim.set("rho", "MC");
    lim.set("U", "MC");
    lim.set("sigma", "MC");
    lim.set("epsP", "superbee");
    lim.set("epsPEq", "superbee");

    const dictionary* dictPtr = materialProperties.findDict("remapLimiters");
    if (dictPtr)
    {
        for (const entry& e : *dictPtr)
        {
            if (!lim.found(e.keyword()))
            {
                FatalIOErrorInFunction(*dictPtr)
                    << "Unknown field " << e.keyword()
                    << " in remapLimiters. Valid: " << lim.sortedToc()
                    << exit(FatalIOError);
            }
            const word name(e.stream());
            remapDetail::limiterFromWord(name, *dictPtr);   // validate
            lim.set(e.keyword(), name);
        }
    }
    return lim;
}


void Foam::eulerianSolid::remap(const scalar dt)
{
    using namespace remapDetail;

    const fvMesh& mesh = mesh_;
    const surfaceScalarField& phi = phiL_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const fvBoundaryMesh& bm = mesh.boundary();
    const scalarField& V = mesh.V().field();
    const scalarField& phiI = phi.primitiveField();

    auto lim = [&](const word& field)
    {
        return limiterFromWord(remapLimiter_[field], materialProperties_);
    };

    const remapGeometry geo(makeGeometry(mesh, phi, dt));

    // Lagrangian volume change dt*sum(phiL out), internal and coupled faces
    scalarField dVL(mesh.nCells(), Zero);
    forAll(own, facei)
    {
        dVL[own[facei]] += dt*phiI[facei];
        dVL[nei[facei]] -= dt*phiI[facei];
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        forAll(fc, i)
        {
            dVL[fc[i]] += dt*pPhi[i];
        }
    }

    // Mass fluxes dm = dt*phiL*rho_f (positive out of the owner / local cell)
    scalarField& rho = rho_.primitiveFieldRef();
    scalarField rhofInt;
    List<scalarField> rhofPatch;
    faceValues(mesh, geo, phi, rho, lim("rho"), rhofInt, rhofPatch);

    scalarField dmInt(dt*phiI*rhofInt);
    List<scalarField> dmPatch(bm.size());
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        dmPatch[patchi] = dt*phi.boundaryField()[patchi]*rhofPatch[patchi];
    }

    scalarField dmOut(mesh.nCells(), Zero);
    forAll(own, facei)
    {
        dmOut[own[facei]] += dmInt[facei];
        dmOut[nei[facei]] -= dmInt[facei];
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        forAll(fc, i)
        {
            dmOut[fc[i]] += dmPatch[patchi][i];
        }
    }

    // New density, written as an increment of the Lagrangian density:
    // (rho_L V_L - dmOut)/V = rho_L + (rho_L (V_L - V) - dmOut)/V
    scalar minRho = GREAT;
    forAll(rho, celli)
    {
        if (dVL[celli] != 0 || dmOut[celli] != 0)
        {
            rho[celli] += (rho[celli]*dVL[celli] - dmOut[celli])/V[celli];
        }
        minRho = min(minRho, rho[celli]);
    }
    reduce(minRho, minOp<scalar>());
    if (minRho <= 0)
    {
        FatalErrorInFunction
            << "Remap gave a non-positive density " << minRho
            << ": the material Courant number is too large for dt = " << dt
            << exit(FatalError);
    }
    rho_.correctBoundaryConditions();

    const scalarField mNew(rho*V);

    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("U"), U_);
    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("sigma"), sigma_);
    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("epsP"), epsP_);
    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("epsPEq"), epsPEq_);
}


// ************************************************************************* //
