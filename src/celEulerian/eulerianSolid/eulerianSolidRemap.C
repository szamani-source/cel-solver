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

    Non-coupled boundary faces (wall, symmetry, free surface, prescribed
    velocity): the Lagrange step's V_L includes their phiL, and the remap
    leaves that volume outside the fixed mesh, so material crosses them with
    the boundary cell's own (Lagrangian) values: outflow removes rho_L dt phiL
    of mass, inflow brings in the cell's own state (first order, zero
    gradient). Walls and symmetry planes have phiL = 0. Coupled faces (processor, cyclic) are handled
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

// Relative size below which a face flux counts as round-off noise
static const scalar fluxTol = 1e-10;

// With flowImpedance (WP 3.6c2, decision 0044): a cell whose void fraction
// is below this is full. In stiff material the top cell of a flowing layer
// opens 1e-9 to 1e-5 of void from round-off and wave noise; without the
// tolerance every such void switches it to the partly filled rules (first
// order, no flow viscosity, leaning slab) and back, which kicks its stress
static const scalar nearlyFull = 1e-4;

// ... but only while the snap keeps the material within this fraction of
// the density its stress implies, rho0 (1 + p/K) (WP 3.6d, decision 0045).
// Each snap adds up to nearlyFull of material volume. Above an open outlet
// the converging flow opens about 1e-6 of void per step; snapped every
// step, that material stayed full and stressed instead of loosening: the
// 1 MPa discharge slowed by a fifth and at 10 MPa the stream inflated to
// half its density
static const scalar nearlyFullRho = 1e-3;

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

    //- Faces whose flux is not round-off: |phiL| > fluxTol * sum|phiL| of
    //  the donor cell. Only these cap the donor's gradient (a face whose
    //  flux is round-off noise would switch the cap on and off at random)
    boolList fluxInt;
    List<boolList> fluxPatch;

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
            b[fc[i]] += pSf[i]*pPhi[i];
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

    // Faces carrying more than round-off flux, relative to the donor cell
    scalarField sumAbs(nCells, Zero);
    forAll(own, facei)
    {
        sumAbs[own[facei]] += mag(phiI[facei]);
        sumAbs[nei[facei]] += mag(phiI[facei]);
    }
    forAll(bm, patchi)
    {
        const labelUList& fc = bm[patchi].faceCells();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        forAll(fc, i)
        {
            sumAbs[fc[i]] += mag(pPhi[i]);
        }
    }
    geo.fluxInt.setSize(own.size());
    forAll(own, facei)
    {
        const label donor = (phiI[facei] > 0 ? own[facei] : nei[facei]);
        geo.fluxInt[facei] = mag(phiI[facei]) > fluxTol*sumAbs[donor];
    }
    geo.fluxPatch.setSize(bm.size());
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        boolList& fp = geo.fluxPatch[patchi];
        fp.setSize(fc.size());
        forAll(fc, i)
        {
            fp[i] = pPhi[i] > fluxTol*sumAbs[fc[i]];
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
        const scalar e = 0.5*(g[celli] & d);
        if (mag(e) <= 0.5*perpTol*mag(g[celli])*mag(d) || mag(e) < ROOTVSMALL)
        {
            // Face (nearly) normal to the gradient, or a vanishing gradient.
            // The absolute test catches a subnormal g, whose relative
            // tolerance underflows to 0 (WP 2.4b: 0.5 g.d rounded to -0 and
            // the ratio divided by zero)
            return;
        }
        const scalar r = ratio(celli, e);
        rMin[celli] = min(rMin[celli], r);
        rMax[celli] = max(rMax[celli], r);
    };

    // Extra cap: the value at each swept centroid this cell donates to stays
    // within the neighbour range (weaker than the midpoint bound in 1D)
    scalarField cap(nCells, GREAT);
    auto addSwept = [&](const label celli, const vector& x)
    {
        const scalar e = g[celli] & x;
        if (mag(e) <= perpTol*mag(g[celli])*mag(x) || mag(e) < ROOTVSMALL)
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
        if (geo.fluxInt[facei])
        {
            addSwept
            (
                phiI[facei] > 0 ? own[facei] : nei[facei],
                geo.rInt[facei]
            );
        }
    }
    forAll(bm, patchi)
    {
        if (!bm[patchi].coupled()) continue;
        const labelUList& fc = bm[patchi].faceCells();
        const vectorField& dp = geo.dPatch[patchi];
        forAll(fc, i)
        {
            addFace(fc[i], dp[i]);
            if (geo.fluxPatch[patchi][i])
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
    List<scalarField>& qfPatch,
    const scalarField* smooth = nullptr
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const fvBoundaryMesh& bm = mesh.boundary();
    const scalarField& phiI = phi.primitiveField();

    const List<scalarField> qNbr(neighbourValues(mesh, q));
    vectorField g(limitedGradient(mesh, geo, phi, q, qNbr, lim));
    if (smooth)
    {
        // Free surface: first order unless the cell and its neighbours
        // are full (decision 0021 rule 6)
        forAll(g, celli)
        {
            if ((*smooth)[celli] < 0.5)
            {
                g[celli] = Zero;
            }
        }
    }
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
    scalarField& q,
    const scalarField* smooth = nullptr
)
{
    scalarField qfInt;
    List<scalarField> qfPatch;
    faceValues(mesh, geo, phi, q, lim, qfInt, qfPatch, smooth);

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
        if (mNew[celli] <= 0)
        {
            q[celli] = 0;       // emptied (free surface only)
        }
        else if (sumF[celli] != 0)
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
    GeometricField<Type, fvPatchField, volMesh>& f,
    const scalarField* smooth = nullptr
)
{
    Field<Type>& fi = f.primitiveFieldRef();
    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        scalarField q(fi.component(cmpt));
        remapMassWeighted
        (
            mesh, geo, phi, dmInt, dmPatch, mNew, lim, q, smooth
        );
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

    if (freeSurface_)
    {
        remapFreeSurface(dt);
        return;
    }

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


void Foam::eulerianSolid::remapFreeSurface(const scalar dt)
{
    // Decision 0031 (decision 0021 rule 6 face by face on a Cartesian mesh).
    // The Lagrange step left rho = m/V_L (bulk), alphaL_ (material fraction
    // of V_L), lean_ and smooth_ (from the alpha at the start of the step).
    using namespace remapDetail;

    const fvMesh& mesh = mesh_;
    const surfaceScalarField& phi = phiL_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const fvBoundaryMesh& bm = mesh.boundary();
    const scalarField& V = mesh.V().field();
    const scalarField& phiI = phi.primitiveField();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& magSf = mesh.magSf().primitiveField();
    const scalarField& dc = mesh.deltaCoeffs().primitiveField();
    const label nCells = mesh.nCells();

    auto lim = [&](const word& field)
    {
        return limiterFromWord(remapLimiter_[field], materialProperties_);
    };
    auto axisOf = [](const vector& S) -> direction
    {
        const vector a(cmptMag(S));
        return a.x() >= a.y() && a.x() >= a.z() ? 0 : (a.y() >= a.z() ? 1 : 2);
    };

    const remapGeometry geo(makeGeometry(mesh, phi, dt));

    // Lagrangian volume and lengths H_a = h_a + dt sum_{faces along a} phi/A_a
    scalarField VL(V);
    vectorField dVa(nCells, Zero);
    forAll(own, facei)
    {
        const direction d = axisOf(Sf[facei]);
        VL[own[facei]] += dt*phiI[facei];
        VL[nei[facei]] -= dt*phiI[facei];
        dVa[own[facei]][d] += dt*phiI[facei];
        dVa[nei[facei]][d] -= dt*phiI[facei];
    }
    forAll(bm, patchi)
    {
        const labelUList& fc = bm[patchi].faceCells();
        const vectorField& pS = bm[patchi].Sf();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        forAll(fc, i)
        {
            VL[fc[i]] += dt*pPhi[i];
            dVa[fc[i]][axisOf(pS[i])] += dt*pPhi[i];
        }
    }
    vectorField H(hCell_);
    forAll(H, celli)
    {
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            H[celli][d] += dVa[celli][d]*hCell_[celli][d]/V[celli];
        }
    }

    scalarField& rho = rho_.primitiveFieldRef();
    const scalarField mL(rho*VL);
    const scalarField AL(alphaL_*VL);
    scalarField rhoM(nCells, Zero);
    forAll(rhoM, celli)
    {
        if (AL[celli] > 0)
        {
            rhoM[celli] = mL[celli]/AL[celli];
        }
    }

    // Material volume in the swept region of length ad at the end of a
    // donor cell (fraction aL, length Hd) along its outward normal, slab
    // side la (+1 towards that end, -1 away, 0 spread)
    auto overlap = []
    (
        const scalar aL, const scalar Hd, const scalar ad, const scalar la
    ) -> scalar
    {
        const scalar A = aL*Hd;
        if (la > 0.5)
        {
            return min(ad, A);
        }
        if (la < -0.5)
        {
            return max(scalar(0), ad - (Hd - A));
        }
        return aL*ad;
    };
    auto sideOf = [&](const vector& lean, const vector& Sout) -> scalar
    {
        const direction d = axisOf(Sout);
        return lean[d]*sign(Sout[d]);
    };

    const List<scalarField> aLNbr(neighbourValues(mesh, alphaL_));
    const List<vectorField> HNbr(neighbourValues(mesh, H));
    const List<vectorField> leanNbr(neighbourValues(mesh, lean_));
    const List<scalarField> rhoMNbr(neighbourValues(mesh, rhoM));
    const List<scalarField> smoothNbr(neighbourValues(mesh, smooth_));

    // Raw material outflow per face, positive, from its donor; total per
    // local donor cell
    scalarField vfInt(own.size(), Zero);
    List<scalarField> vfPatch(bm.size());
    scalarField outTot(nCells, Zero);
    forAll(own, facei)
    {
        if (phiI[facei] == 0) continue;
        const bool fromP = phiI[facei] > 0;
        const label D = fromP ? own[facei] : nei[facei];
        const vector Sout(fromP ? Sf[facei] : -Sf[facei]);
        const scalar ad = dt*mag(phiI[facei])/magSf[facei];
        const scalar v =
            overlap
            (
                alphaL_[D], H[D][axisOf(Sout)], ad, sideOf(lean_[D], Sout)
            )*magSf[facei];
        vfInt[facei] = v;
        outTot[D] += v;
    }
    forAll(bm, patchi)
    {
        const fvPatch& p = bm[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pS = p.Sf();
        const scalarField& pmS = p.magSf();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        scalarField& vp = vfPatch[patchi];
        vp.setSize(fc.size(), Zero);
        forAll(fc, i)
        {
            if (pPhi[i] == 0) continue;
            const scalar ad = dt*mag(pPhi[i])/pmS[i];
            if (pPhi[i] > 0)
            {
                const label D = fc[i];
                vp[i] =
                    overlap
                    (
                        alphaL_[D], H[D][axisOf(pS[i])], ad,
                        sideOf(lean_[D], pS[i])
                    )*pmS[i];
                outTot[D] += vp[i];
            }
            else if (p.coupled())
            {
                const vector Sout(-pS[i]);
                vp[i] =
                    overlap
                    (
                        aLNbr[patchi][i], HNbr[patchi][i][axisOf(Sout)], ad,
                        sideOf(leanNbr[patchi][i], Sout)
                    )*pmS[i];
            }
            // inflow through other boundary faces: nothing outside
        }
    }

    // A donor that would keep less than 1e-9 of its volume gives all of it
    scalarField scale(nCells, scalar(1));
    forAll(scale, celli)
    {
        if
        (
            AL[celli] > 0 && outTot[celli] > 0
         && AL[celli] - outTot[celli] < 1e-9*V[celli]
        )
        {
            scale[celli] = AL[celli]/outTot[celli];
        }
    }
    const List<scalarField> scaleNbr(neighbourValues(mesh, scale));

    // Signed material volume and mass fluxes (out of the owner / local
    // cell); fluxes below 1e-9 of a cell are dropped. Mass: MUSCL density
    // where the donor is smooth, else the donor's Lagrangian material
    // density.
    scalarField rhofInt;
    List<scalarField> rhofPatch;
    faceValues(mesh, geo, phi, rho, lim("rho"), rhofInt, rhofPatch);

    scalarField dmInt(own.size(), Zero);
    List<scalarField> dmPatch(bm.size());
    scalarField dmOut(nCells, Zero);
    scalarField vOut(nCells, Zero);
    forAll(own, facei)
    {
        if (phiI[facei] == 0) continue;
        const bool fromP = phiI[facei] > 0;
        const label D = fromP ? own[facei] : nei[facei];
        scalar v = vfInt[facei]*scale[D];
        if (v <= 1e-9*magSf[facei]/dc[facei])
        {
            continue;
        }
        const scalar dm =
            smooth_[D] > 0.5
          ? dt*mag(phiI[facei])*rhofInt[facei]
          : v*rhoM[D];
        const scalar sgn = fromP ? 1 : -1;
        dmInt[facei] = sgn*dm;
        dmOut[own[facei]] += sgn*dm;
        dmOut[nei[facei]] -= sgn*dm;
        vOut[own[facei]] += sgn*v;
        vOut[nei[facei]] -= sgn*v;
    }
    forAll(bm, patchi)
    {
        const fvPatch& p = bm[patchi];
        const labelUList& fc = p.faceCells();
        const scalarField& pPhi = phi.boundaryField()[patchi];
        const scalarField& pmS = p.magSf();
        const scalarField pdc(p.deltaCoeffs());
        dmPatch[patchi].setSize(fc.size(), Zero);
        forAll(fc, i)
        {
            if (pPhi[i] == 0) continue;
            const bool local = pPhi[i] > 0;
            const scalar v =
                vfPatch[patchi][i]
               *(local ? scale[fc[i]] : (p.coupled() ? scaleNbr[patchi][i] : 0));
            if (v <= 1e-9*pmS[i]/pdc[i])
            {
                continue;
            }
            scalar dm;
            if (local)
            {
                dm =
                    (p.coupled() && smooth_[fc[i]] > 0.5)
                  ? dt*pPhi[i]*rhofPatch[patchi][i]
                  : v*rhoM[fc[i]];
            }
            else
            {
                dm =
                    smoothNbr[patchi][i] > 0.5
                  ? dt*mag(pPhi[i])*rhofPatch[patchi][i]
                  : v*rhoMNbr[patchi][i];
            }
            const scalar sgn = local ? 1 : -1;
            if (p.coupled())
            {
                dmPatch[patchi][i] = sgn*dm;
            }
            else
            {
                // Leaves the domain (an open outlet, WP 3.3)
                discharged_ += dm;
            }
            dmOut[fc[i]] += sgn*dm;
            vOut[fc[i]] += sgn*v;
        }
    }

    // New mass and volume fraction; round-off dust (alpha < 1e-9 with
    // mass) gets its volume from its mass and no velocity or stress
    scalarField& alpha = alpha_.primitiveFieldRef();
    scalarField mNew(mL - dmOut);
    const scalar Kbulk =
        material().waveModulus() - 4.0*material().shearModulus()/3.0;
    boolList dust(nCells, false);
    forAll(mNew, celli)
    {
        if (mNew[celli] <= 0)
        {
            mNew[celli] = 0;
            alpha[celli] = 0;
            rho[celli] = 0;
            continue;
        }
        scalar a = (AL[celli] - vOut[celli])/V[celli];
        a = min(max(a, scalar(0)), scalar(1));
        if
        (
            a >= 1 - 1e-12
         || (
                flowImpedance_
             && a >= 1 - nearlyFull
             && mNew[celli]
             >= (1 - nearlyFullRho)*rho0_*V[celli]
               *(1 + max(-tr(sigma_[celli])/3, scalar(0))/Kbulk)
            )
        )
        {
            a = 1;
        }
        if (a < 1e-9)
        {
            a = mNew[celli]/(rho0_*V[celli]);
            dust[celli] = true;
        }
        alpha[celli] = a;
        rho[celli] = mNew[celli]/V[celli];
    }
    rho_.correctBoundaryConditions();
    alpha_.correctBoundaryConditions();

    // dmPatch for non-coupled patches stays zero: remapMassWeighted skips
    // them, and the boundary donor value is the cell's own (first order)
    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("U"), U_, &smooth_);
    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("sigma"), sigma_, &smooth_);
    remapField(mesh, geo, phi, dmInt, dmPatch, mNew, lim("epsP"), epsP_, &smooth_);
    remapField
    (
        mesh, geo, phi, dmInt, dmPatch, mNew, lim("epsPEq"), epsPEq_, &smooth_
    );

    vectorField& U = U_.primitiveFieldRef();
    symmTensorField& sigma = sigma_.primitiveFieldRef();
    forAll(U, celli)
    {
        if (dust[celli] || mNew[celli] <= 0)
        {
            U[celli] = Zero;
            sigma[celli] = Zero;
        }
    }
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();
}


// ************************************************************************* //
