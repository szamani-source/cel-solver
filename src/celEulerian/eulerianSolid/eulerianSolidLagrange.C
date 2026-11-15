/*---------------------------------------------------------------------------*\
    eulerianSolid::lagrangeStep (WP 2.3)

    The scheme of decision 0013, unsplit on any fvMesh:

    1. Face states. With secondOrder, each of the nine primitive components
       (U and sigma) gets a Gauss cell gradient; the state on each side of a
       face is the cell value plus half the MC-limited slope along the line
       joining the two cell centres (on a uniform mesh exactly the 1D MC
       MUSCL slope). The half-step (Hancock) predictor adds, per cell,
       dt/2 div(sigma)/rho to U and dt/2 C:sym(grad U) to sigma, both from
       the cell's own limited face states and the elastic moduli.
       Non-coupled boundary faces reconstruct towards the patch value of
       the field through a mirror ghost cell: exact for a linear velocity
       with fixedValue patches, the cell value on zeroGradient patches.
    2. Acoustic Riemann solver per face in the face-normal frame: normal
       velocity and normal traction with the P impedance rho c_p, tangential
       velocity and shear traction with the S impedance rho c_s. Result: face
       velocity vStar and the traction tStar = sigma* & n on the owner side.
       Boundary faces from the patch type of U (decision 0018). With
       solveMomentum no (prescribed motion) the face velocity is the mean of
       the two face states, without the stress-jump term (decision 0022).
    3. phiL = vStar & Sf (decision 0014); face force F = tStar |Sf|.
       Coupled faces (cyclic, processor, processorCyclic) are solved on the
       owner side only; the other side receives vStar and -F, -phiL, so both
       sides hold bit-identical values and momentum is exchanged exactly.
    4. Momentum: m dU = dt sum(F) (if solveMomentum).
    5. Velocity gradient by Gauss from vStar over the mid-step volume;
       strain increment dt sym(L), spin dt skew(L); sigma and epsP rotated
       with Hughes-Winget, one model call per cell (none if the strain
       increment is exactly zero).
    6. rho = m / V_L, V_L = V + dt sum(phiL out of the cell).

    Gravity and walls (WP 3.1, decision 0030). Gravity g(t) is balanced in
    the face states: across a face the stress reconstruction works on the
    deviation from each cell's own hydrostatic jump
    e = -rho (g . d) gHat gHat (d the vector between the cells), so a
    column in discrete equilibrium (cell stresses differing by the mean of
    the two cells' e) has equal face states on both sides and stays at
    rest to round-off; the momentum and Hancock predictor add g. Wall
    patches (materialProperties walls): v_n = 0 unless that gives tension
    (then zero traction); sticking shear kept if at most mu times the wall
    pressure, else clipped to it with v_t from the interior characteristic.
    Damping U <- U/(1 + damping dt) after the momentum update.

    Free surface (WP 3.2a, decision 0031; materialProperties freeSurface):
    one material plus void on a Cartesian hex mesh, decision 0021 rules 1
    to 7 face by face. Partly filled cells are homogenised (alpha sigma,
    alpha z) and first order; each leans along each axis towards its fuller
    neighbour; faces the material does not touch on both sides are
    traction-free (also a pulling contact face of material without tensile
    strength); compression closes void first, expansion unloads the
    material and then opens void. alphaL_, lean_ and smooth_ go to the remap.
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "processorPolyPatch.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "PstreamBuffers.H"
#include "calculatedFvPatchFields.H"

#include "transform.H"

#include <map>
#include <tuple>

// * * * * * * * * * * * * * * * Local Helpers * * * * * * * * * * * * * * //

namespace Foam
{
namespace lagrangeStepDetail
{

//- Minmod of three numbers: zero unless all have the same sign, else the
//  one of smallest magnitude
inline scalar minmod3(const scalar a, const scalar b, const scalar c)
{
    if (a > 0 && b > 0 && c > 0)
    {
        return min(a, min(b, c));
    }
    if (a < 0 && b < 0 && c < 0)
    {
        return max(a, max(b, c));
    }
    return 0;
}


//- MC-limited slope (change of q over the cell, along the line to the
//  other cell) from the jump a = q_other - q_this and the cell-gradient
//  increment c = (x_other - x_this) & grad(q)_this, with
//  the far-side jump 2c - a. Symmetric in the two one-sided jumps.
inline scalar mcSlope(const scalar a, const scalar c)
{
    return minmod3(2*a, c, 2*(2*c - a));
}


//- Hydrostatic stress jump over the vector d in a cell of density rho
//  under gravity g: -rho (g . d) gHat gHat for the stress, zero for other
//  types (decision 0030)
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
    return (-rho*(g & d)/gg)*sqr(g);
}


//- A vol field of Type with calculated patches (constraint patches keep
//  their constraint type, so coupled patches exchange neighbour values)
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> workField
(
    const fvMesh& mesh,
    const word& name,
    const dimensionSet& dims
)
{
    return GeometricField<Type, fvPatchField, volMesh>::New
    (
        name,
        IOobject::NO_REGISTER,
        mesh,
        dimensioned<Type>(dims, Zero),
        fvPatchFieldBase::calculatedType()
    );
}


//- Gauss gradient (linear interpolation to faces) of each component of q.
//  Non-coupled boundary faces take the patch value of q. Coupled patches of the
//  results are evaluated, so patchNeighbourField() gives the neighbour's
//  gradient.
template<class Type>
void componentGradients
(
    const GeometricField<Type, fvPatchField, volMesh>& q,
    PtrList<volVectorField>& g
)
{
    const fvMesh& mesh = q.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceScalarField& w = mesh.weights();
    const surfaceVectorField& Sf = mesh.Sf();
    const direction nCmpt = pTraits<Type>::nComponents;

    g.setSize(nCmpt);
    for (direction d = 0; d < nCmpt; ++d)
    {
        g.set
        (
            d,
            workField<vector>(mesh, "grad" + Foam::name(d), dimless/dimLength)
        );
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const Type qf = w[facei]*q[P] + (1 - w[facei])*q[N];
        for (direction d = 0; d < nCmpt; ++d)
        {
            const vector s = Sf[facei]*component(qf, d);
            g[d][P] += s;
            g[d][N] -= s;
        }
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const fvPatchField<Type>& pq = q.boundaryField()[patchi];

        if (p.coupled())
        {
            const scalarField& pw = w.boundaryField()[patchi];
            const Field<Type> qn(pq.patchNeighbourField());
            forAll(fc, i)
            {
                const Type qf = pw[i]*q[fc[i]] + (1 - pw[i])*qn[i];
                for (direction d = 0; d < nCmpt; ++d)
                {
                    g[d][fc[i]] += pSf[i]*component(qf, d);
                }
            }
        }
        else
        {
            forAll(fc, i)
            {
                for (direction d = 0; d < nCmpt; ++d)
                {
                    g[d][fc[i]] += pSf[i]*component(pq[i], d);
                }
            }
        }
    }

    for (direction d = 0; d < nCmpt; ++d)
    {
        g[d].primitiveFieldRef() /= mesh.V().field();
        g[d].correctBoundaryConditions();
    }
}


//- Limited face states of q on the owner side (qO) and neighbour side (qN)
//  of every internal and coupled face. Non-coupled boundary faces: qO from a
//  mirror ghost cell at 2 (Cf - C) holding 2 q_b - q_P (q_b the patch value),
//  so a linear field with exact patch values is reconstructed exactly and a
//  zeroGradient patch gives the cell value; qN = qO. With secondOrder false
//  all states are the cell values.
//  With gravity g nonzero and rho given, the reconstruction works on the
//  deviation from each side's hydrostatic jump (equilibriumJump): the
//  jump a and the gradient increments lose the equilibrium part, which
//  is added back when extrapolating to the face (decision 0030).
template<class Type>
void faceStates
(
    const GeometricField<Type, fvPatchField, volMesh>& q,
    const bool secondOrder,
    GeometricField<Type, fvsPatchField, surfaceMesh>& qO,
    GeometricField<Type, fvsPatchField, surfaceMesh>& qN,
    const volScalarField* rho = nullptr,
    const vector& gravity = Zero,
    const volScalarField* smooth = nullptr,
    const volScalarField* eFac = nullptr
)
{
    const bool grav = rho && magSqr(gravity) > 0;
    const Type zeroT(Zero);
    const fvMesh& mesh = q.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const direction nCmpt = pTraits<Type>::nComponents;

    PtrList<volVectorField> g;
    if (secondOrder)
    {
        componentGradients(q, g);
    }
    const tmp<surfaceVectorField> tdelta(mesh.delta());
    const surfaceVectorField& delta = tdelta();

    // Owner and neighbour side states from cell values qP, qQ, gradients
    // and the vector dPQ from owner to neighbour centre
    // eP, eQ: hydrostatic jumps of the two cells over dPQ (zero without
    // gravity, which leaves the original scheme unchanged)
    auto states0 = [&]
    (
        const Type& qP, const Type& qQ,
        const UList<vector>* gP, const UList<vector>* gQ,
        const vector& dPQ,
        Type& outP, Type& outQ,
        const Type& eP, const Type& eQ
    )
    {
        outP = qP;
        outQ = qQ;
        if (!secondOrder)
        {
            if (grav)
            {
                outP += 0.5*eP;
                outQ -= 0.5*eQ;
            }
            return;
        }
        if (!grav)
        {
            for (direction d = 0; d < nCmpt; ++d)
            {
                const scalar a = component(qQ, d) - component(qP, d);
                const scalar sP = mcSlope(a, dPQ & (*gP)[d]);
                const scalar sQ = mcSlope(a, dPQ & (*gQ)[d]);
                setComponent(outP, d) = component(qP, d) + 0.5*sP;
                setComponent(outQ, d) = component(qQ, d) - 0.5*sQ;
            }
            return;
        }
        for (direction d = 0; d < nCmpt; ++d)
        {
            const scalar ep = component(eP, d);
            const scalar eq = component(eQ, d);
            const scalar a =
                component(qQ, d) - component(qP, d) - 0.5*(ep + eq);
            const scalar sP = mcSlope(a, (dPQ & (*gP)[d]) - ep);
            const scalar sQ = mcSlope(a, (dPQ & (*gQ)[d]) - eq);
            setComponent(outP, d) = component(qP, d) + 0.5*(sP + ep);
            setComponent(outQ, d) = component(qQ, d) - 0.5*(sQ + eq);
        }
    };

    // With the free surface (decision 0031): a side whose cell or a
    // neighbour of it is not full is first order; under gravity it keeps
    // the fraction f = 1 - phi of its hydrostatic jump (rule 7 of decision
    // 0021). The other side as usual. Without, states0 unchanged.
    auto states = [&]
    (
        const Type& qP, const Type& qQ,
        const UList<vector>* gP, const UList<vector>* gQ,
        const vector& dPQ,
        Type& outP, Type& outQ,
        const Type& eP, const Type& eQ,
        const scalar smP = 1, const scalar smQ = 1,
        const scalar fP = 1, const scalar fQ = 1
    )
    {
        if (smP > 0.5 && smQ > 0.5)
        {
            states0(qP, qQ, gP, gQ, dPQ, outP, outQ, eP, eQ);
            return;
        }
        Type tP(qP), tQ(qQ);
        if (smP > 0.5 || smQ > 0.5)
        {
            states0(qP, qQ, gP, gQ, dPQ, tP, tQ, eP, eQ);
        }
        outP = smP > 0.5 ? tP : (grav ? Type(qP + (0.5*fP)*eP) : qP);
        outQ = smQ > 0.5 ? tQ : (grav ? Type(qQ - (0.5*fQ)*eQ) : qQ);
    };

    auto sm = [&](const label celli) -> scalar
    {
        return smooth ? (*smooth)[celli] : 1;
    };
    auto ef = [&](const label celli) -> scalar
    {
        return eFac ? (*eFac)[celli] : 1;
    };

    auto eqJump = [&](const scalar r, const vector& d) -> Type
    {
        return grav ? equilibriumJump<Type>(r, d, gravity) : zeroT;
    };

    List<vector> gP(nCmpt), gQ(nCmpt);

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        if (secondOrder)
        {
            for (direction d = 0; d < nCmpt; ++d)
            {
                gP[d] = g[d][P];
                gQ[d] = g[d][N];
            }
        }
        states
        (
            q[P], q[N], &gP, &gQ, delta[facei], qO[facei], qN[facei],
            eqJump(grav ? (*rho)[P] : 0, delta[facei]),
            eqJump(grav ? (*rho)[N] : 0, delta[facei]),
            sm(P), sm(N), ef(P), ef(N)
        );
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        Field<Type>& pO = qO.boundaryFieldRef()[patchi];
        Field<Type>& pN = qN.boundaryFieldRef()[patchi];

        if (p.coupled())
        {
            const Field<Type> qn(q.boundaryField()[patchi].patchNeighbourField());
            const vectorField& pd = delta.boundaryField()[patchi];
            const scalarField rn
            (
                grav
              ? rho->boundaryField()[patchi].patchNeighbourField()
              : tmp<scalarField>::New(p.size(), Zero)
            );
            const scalarField smN
            (
                smooth
              ? smooth->boundaryField()[patchi].patchNeighbourField()
              : tmp<scalarField>::New(p.size(), scalar(1))
            );
            const scalarField efN
            (
                eFac
              ? eFac->boundaryField()[patchi].patchNeighbourField()
              : tmp<scalarField>::New(p.size(), scalar(1))
            );

            PtrList<vectorField> gn(secondOrder ? nCmpt : 0);
            if (secondOrder)
            {
                for (direction d = 0; d < nCmpt; ++d)
                {
                    gn.set
                    (
                        d,
                        g[d].boundaryField()[patchi].patchNeighbourField()
                    );
                }
            }

            forAll(fc, i)
            {
                if (secondOrder)
                {
                    for (direction d = 0; d < nCmpt; ++d)
                    {
                        gP[d] = g[d][fc[i]];
                        gQ[d] = gn[d][i];
                    }
                }
                states
                (
                    q[fc[i]], qn[i], &gP, &gQ, pd[i], pO[i], pN[i],
                    eqJump(grav ? (*rho)[fc[i]] : 0, pd[i]),
                    eqJump(rn[i], pd[i]),
                    sm(fc[i]), smN[i], ef(fc[i]), efN[i]
                );
            }
        }
        else
        {
            const Field<Type>& qb = q.boundaryField()[patchi];
            const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
            const vectorField& C = mesh.C().primitiveField();
            forAll(fc, i)
            {
                if (secondOrder)
                {
                    for (direction d = 0; d < nCmpt; ++d)
                    {
                        gP[d] = g[d][fc[i]];
                    }
                }
                // The mirror ghost carries the cell's own hydrostatic
                // jump, so the deviation is mirrored (decision 0030)
                Type ghost;
                const vector dG(2*(Cf[i] - C[fc[i]]));
                const Type e(eqJump(grav ? (*rho)[fc[i]] : 0, dG));
                states
                (
                    q[fc[i]], 2*qb[i] - q[fc[i]] + e, &gP, &gP,
                    dG, pO[i], ghost, e, e,
                    sm(fc[i]), sm(fc[i]), ef(fc[i]), ef(fc[i])
                );
                pN[i] = pO[i];
            }
        }
    }
}


//- Acoustic Riemann solution at a face with unit normal n (from L to R):
//  face velocity vS and traction tS = sigma* & n
inline void riemann
(
    const vector& n,
    const vector& vL, const symmTensor& sL, const scalar zpL, const scalar zsL,
    const vector& vR, const symmTensor& sR, const scalar zpR, const scalar zsR,
    vector& vS,
    vector& tS
)
{
    const vector tL(sL & n);
    const vector tR(sR & n);

    const scalar vnL = vL & n;
    const scalar vnR = vR & n;
    const scalar tnL = tL & n;
    const scalar tnR = tR & n;
    const vector vtL(vL - vnL*n);
    const vector vtR(vR - vnR*n);
    const vector ttL(tL - tnL*n);
    const vector ttR(tR - tnR*n);

    const scalar zp = zpL + zpR;
    const scalar zs = zsL + zsR;

    const scalar vn = (zpL*vnL + zpR*vnR + (tnR - tnL))/zp;
    const scalar tn = (zpR*tnL + zpL*tnR + zpL*zpR*(vnR - vnL))/zp;
    const vector vt((zsL*vtL + zsR*vtR + (ttR - ttL))/zs);
    const vector tt((zsR*ttL + zsL*ttR + zsL*zsR*(vtR - vtL))/zs);

    vS = vn*n + vt;
    tS = tn*n + tt;
}


//- Boundary kinds (decision 0018)
enum class boundaryKind
{
    coupled, empty, velocity, symmetry, traction, wall
};


boundaryKind kindOf
(
    const fvPatch& p,
    const fvPatchVectorField& pU,
    const HashTable<scalar>& walls
)
{
    if (walls.found(p.name()))
    {
        if (p.coupled() || isA<emptyFvPatch>(p))
        {
            FatalErrorInFunction
                << "Patch " << p.name() << " of type " << p.type()
                << " cannot be a wall" << exit(FatalError);
        }
        return boundaryKind::wall;
    }
    if (p.coupled())
    {
        if (isA<processorFvPatch>(p) || isA<cyclicFvPatch>(p))
        {
            return boundaryKind::coupled;
        }
        FatalErrorInFunction
            << "Patch " << p.name() << " of type " << p.type()
            << ": coupled patches other than cyclic and processor are not "
            << "supported by the Lagrange step" << exit(FatalError);
    }
    if (isA<emptyFvPatch>(p))
    {
        return boundaryKind::empty;
    }
    if (isA<wedgeFvPatch>(p))
    {
        FatalErrorInFunction
            << "Patch " << p.name() << ": wedge is not supported"
            << exit(FatalError);
    }
    if (isA<symmetryFvPatch>(p) || isA<symmetryPlaneFvPatch>(p))
    {
        return boundaryKind::symmetry;
    }
    if (pU.fixesValue())
    {
        return boundaryKind::velocity;
    }
    return boundaryKind::traction;
}


//- Riemann solution at a rigid wall at rest with unit normal n out of the
//  cell, from the cell-side state (vL, sL): v_n = 0 unless that pulls on
//  the wall (then zero traction); sticking shear clipped to mu |t_n|
//  (decision 0030). Returns face velocity v and traction t = sigma* & n.
inline void wallFace
(
    const vector& n,
    const vector& vL, const symmTensor& sL,
    const scalar zp, const scalar zs, const scalar mu,
    vector& v, vector& t
)
{
    const vector tL(sL & n);
    const scalar vnL = vL & n;
    const scalar tnL = tL & n;
    const vector vtL(vL - vnL*n);
    const vector ttL(tL - tnL*n);

    const scalar tn = tnL - zp*vnL;
    if (tn >= 0)
    {
        v = (vnL - tnL/zp)*n + vtL - ttL/zs;
        t = Zero;
        return;
    }
    vector tt(ttL - zs*vtL);
    const scalar limit = -tn*mu;
    const scalar magTt = mag(tt);
    v = Zero;
    if (magTt > limit)
    {
        tt *= limit/magTt;
        v = vtL + (tt - ttL)/zs;
    }
    t = tn*n + tt;
}


//- Copy the owner-side coupled face values to the other side:
//  vStar as is, F and phi negated. Cyclic locally, processor by message.
//  With the shear cap (dO, dN not null): the slip velocities seen by the
//  owner and the neighbour cell swap sides.
void syncCoupled
(
    const fvMesh& mesh,
    surfaceVectorField& vS,
    surfaceVectorField& F,
    surfaceScalarField& phi,
    surfaceVectorField* dO = nullptr,
    surfaceVectorField* dN = nullptr
)
{
    auto& bv = vS.boundaryFieldRef();
    auto& bF = F.boundaryFieldRef();
    auto& bp = phi.boundaryFieldRef();
    const bool slip = dO && dN;

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (isA<cyclicFvPatch>(p))
        {
            const cyclicFvPatch& cp = refCast<const cyclicFvPatch>(p);
            if (!cp.owner())
            {
                const label nbr = cp.neighbPatchID();
                const vectorField v(bv[nbr]);
                const vectorField f(-bF[nbr]);
                const scalarField ph(-bp[nbr]);
                bv[patchi] = v;
                bF[patchi] = f;
                bp[patchi] = ph;
                if (slip)
                {
                    const vectorField so(dN->boundaryField()[nbr]);
                    const vectorField sn(dO->boundaryField()[nbr]);
                    dO->boundaryFieldRef()[patchi] = so;
                    dN->boundaryFieldRef()[patchi] = sn;
                }
            }
        }
    }

    if (!UPstream::parRun())
    {
        return;
    }

    // One message per neighbour processor holding every owner-side patch
    // to it, keyed by the patch tag (unique per processor pair: the plain
    // processor patch and one per cyclic for processorCyclic)
    labelHashSet nbrProcs;
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (isA<processorFvPatch>(p))
        {
            nbrProcs.insert
            (
                refCast<const processorFvPatch>(p).neighbProcNo()
            );
        }
    }

    PstreamBuffers pBufs(UPstream::commsTypes::nonBlocking);

    for (const label proci : nbrProcs.sortedToc())
    {
        label n = 0;
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (isA<processorFvPatch>(p))
            {
                const processorPolyPatch& pp =
                    refCast<const processorFvPatch>(p).procPolyPatch();
                if (pp.neighbProcNo() == proci && pp.owner())
                {
                    ++n;
                }
            }
        }

        UOPstream os(proci, pBufs);
        os << n;
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (isA<processorFvPatch>(p))
            {
                const processorPolyPatch& pp =
                    refCast<const processorFvPatch>(p).procPolyPatch();
                if (pp.neighbProcNo() == proci && pp.owner())
                {
                    os  << label(pp.tag())
                        << vectorField(bv[patchi])
                        << vectorField(bF[patchi])
                        << scalarField(bp[patchi]);
                    if (slip)
                    {
                        os  << vectorField(dO->boundaryField()[patchi])
                            << vectorField(dN->boundaryField()[patchi]);
                    }
                }
            }
        }
    }

    pBufs.finishedSends();

    for (const label proci : nbrProcs.sortedToc())
    {
        UIPstream is(proci, pBufs);
        label n;
        is >> n;

        std::map
        <
            label,
            std::tuple<vectorField, vectorField, scalarField, vectorField, vectorField>
        > got;
        for (label k = 0; k < n; ++k)
        {
            label tag;
            vectorField v, f, so, sn;
            scalarField ph;
            is >> tag >> v >> f >> ph;
            if (slip)
            {
                is >> so >> sn;
            }
            got[tag] = std::make_tuple(v, f, ph, so, sn);
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!isA<processorFvPatch>(p))
            {
                continue;
            }
            const processorPolyPatch& pp =
                refCast<const processorFvPatch>(p).procPolyPatch();
            if (pp.neighbProcNo() != proci || pp.owner())
            {
                continue;
            }

            const auto iter = got.find(label(pp.tag()));
            if (iter == got.end() || std::get<0>(iter->second).size() != p.size())
            {
                FatalErrorInFunction
                    << "No matching face values from processor " << proci
                    << " for patch " << p.name() << exit(FatalError);
            }
            const vectorField f(-std::get<1>(iter->second));
            const scalarField ph(-std::get<2>(iter->second));
            bv[patchi] = std::get<0>(iter->second);
            bF[patchi] = f;
            bp[patchi] = ph;
            if (slip)
            {
                // the sender's owner cell is this side's neighbour
                dO->boundaryFieldRef()[patchi] = std::get<4>(iter->second);
                dN->boundaryFieldRef()[patchi] = std::get<3>(iter->second);
            }
        }
    }
}


//- Hughes-Winget rotation for the spin increment dW:
//  Q = (I - dW/2)^-1 (I + dW/2)
inline tensor hughesWinget(const tensor& dW)
{
    return inv(I - 0.5*dW) & (I + 0.5*dW);
}


// Free surface (WP 3.2a, decision 0031) ---------------------------------------

//- Volume fraction thresholds (decision 0021 rule 1)
static const scalar alphaEmpty = 1e-12;
static const scalar alphaFull = 1 - 1e-12;


//- Axis of a face normal to an axis (Cartesian mesh, checked at start)
inline direction faceAxis(const vector& S)
{
    const vector a(cmptMag(S));
    return a.x() >= a.y() && a.x() >= a.z() ? 0 : (a.y() >= a.z() ? 1 : 2);
}


//- Slab side per axis (+1 towards +x_a, -1 towards -x_a, 0 spread or
//  full) and the smooth flag (cell and every neighbour full), decision
//  0021 rule 2 face by face: the slab leans towards the fuller of the two
//  neighbours along each axis; walls, symmetry planes and prescribed
//  velocity patches mirror the cell, traction-free patches count as empty.
//  alpha must have its coupled patches evaluated.
//  Rule 9 (decision 0024, C++ in decision 0032) when gAxis >= 0 (the axis
//  of gravity): a partly filled cell with alpha >= 0.1 is spread (lean 0)
//  along each other axis d where its Youngs (1-2-1) alpha gradient is
//  larger along gAxis than along d.
void slabSides
(
    const fvMesh& mesh,
    const volScalarField& alpha,
    const UList<boundaryKind>& kinds,
    vectorField& lean,
    scalarField& smooth,
    const label gAxis = -1
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf().primitiveField();
    const scalarField& a = alpha.primitiveField();

    // Fraction beyond the + and - face along each axis (default: mirror)
    vectorField aPlus(mesh.nCells()), aMinus(mesh.nCells());
    forAll(a, celli)
    {
        aPlus[celli] = a[celli]*vector::one;
        aMinus[celli] = aPlus[celli];
    }
    auto set = [&](const label c, const vector& Sout, const scalar ab)
    {
        const direction d = faceAxis(Sout);
        if (Sout[d] > 0)
        {
            aPlus[c][d] = ab;
        }
        else
        {
            aMinus[c][d] = ab;
        }
    };
    forAll(own, facei)
    {
        set(own[facei], Sf[facei], a[nei[facei]]);
        set(nei[facei], -Sf[facei], a[own[facei]]);
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pS = p.Sf();
        if (kinds[patchi] == boundaryKind::coupled)
        {
            const scalarField an
            (
                alpha.boundaryField()[patchi].patchNeighbourField()
            );
            forAll(fc, i)
            {
                set(fc[i], pS[i], an[i]);
            }
        }
        else if (kinds[patchi] == boundaryKind::traction)
        {
            forAll(fc, i)
            {
                set(fc[i], pS[i], 0);
            }
        }
    }

    lean.setSize(mesh.nCells());
    smooth.setSize(mesh.nCells());
    forAll(a, celli)
    {
        const bool full = a[celli] >= alphaFull;
        const vector dif(aPlus[celli] - aMinus[celli]);
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            lean[celli][d] =
                (full || mag(dif[d]) <= 1e-9) ? 0 : sign(dif[d]);
        }
        smooth[celli] =
            (
                full
             && cmptMin(aPlus[celli]) >= alphaFull
             && cmptMin(aMinus[celli]) >= alphaFull
            )
          ? 1 : 0;
    }

    if (gAxis < 0)
    {
        return;
    }

    // Rule 9: centred alpha differences D per axis (aPlus - aMinus), then
    // the 1-2-1 sum of D_d over the neighbours along e (e != d); the
    // neighbour beyond a non-coupled boundary is the cell itself (copy)
    tmp<volVectorField> tD = workField<vector>(mesh, "alphaDif", dimless);
    volVectorField& D = tD.ref();
    D.primitiveFieldRef() = aPlus - aMinus;
    D.correctBoundaryConditions();
    const vectorField& Dc = D.primitiveField();

    // Neighbour D along each axis and side: index 2 e + (0 for +, 1 for -)
    List<FixedList<vector, 6>> Dn(mesh.nCells());
    forAll(Dn, celli)
    {
        Dn[celli] = Dc[celli];
    }
    auto setN = [&](const label c, const vector& Sout, const vector& Db)
    {
        const direction e = faceAxis(Sout);
        Dn[c][2*e + (Sout[e] > 0 ? 0 : 1)] = Db;
    };
    forAll(own, facei)
    {
        setN(own[facei], Sf[facei], Dc[nei[facei]]);
        setN(nei[facei], -Sf[facei], Dc[own[facei]]);
    }
    forAll(mesh.boundary(), patchi)
    {
        if (kinds[patchi] != boundaryKind::coupled)
        {
            continue;
        }
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pS = p.Sf();
        const vectorField dn(D.boundaryField()[patchi].patchNeighbourField());
        forAll(fc, i)
        {
            setN(fc[i], pS[i], dn[i]);
        }
    }

    const direction g = direction(gAxis);
    const Vector<label>& geomD = mesh.geometricD();
    forAll(a, celli)
    {
        if (a[celli] < 0.1 || a[celli] >= alphaFull)
        {
            continue;
        }
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            if (d == g || geomD[d] < 0 || lean[celli][d] == 0)
            {
                continue;
            }
            // gradient along d, smoothed along g; along g, smoothed along d
            const scalar Gd =
                2*Dc[celli][d] + Dn[celli][2*g][d] + Dn[celli][2*g + 1][d];
            const scalar Gg =
                2*Dc[celli][g] + Dn[celli][2*d][g] + Dn[celli][2*d + 1][g];
            if (mag(Gg) > mag(Gd))
            {
                lean[celli][d] = 0;
            }
        }
    }
}


//- Does the material of a cell touch its face with outward area vector
//  Sout (non-empty and its slab not leaning away)?
inline bool touches(const scalar a, const vector& lean, const vector& Sout)
{
    if (a <= alphaEmpty)
    {
        return false;
    }
    if (a >= alphaFull)
    {
        return true;
    }
    const direction d = faceAxis(Sout);
    return lean[d]*sign(Sout[d]) > -0.5;
}


//- Velocity of a traction-free face from the one-sided (characteristic)
//  relation of one side. n from L to R; side +1: the L state (vL, sL)
//  with impedances of L; side -1: the R state.
inline vector oneSided
(
    const vector& n, const label side,
    const vector& v, const symmTensor& s, const scalar zp, const scalar zs
)
{
    const vector t(s & n);
    const scalar vn = v & n;
    const scalar tn = t & n;
    const vector vt(v - vn*n);
    const vector tt(t - tn*n);
    return (vn - side*tn/zp)*n + vt - side*tt/zs;
}


//- Free-surface face (decision 0021 rules 3 and 4): contact face (both
//  sides touch) gets the Riemann solution unless the material has no
//  tensile strength and the solution pulls; every other face is
//  traction-free with the one-sided velocity of the side that touches it
//  (owner first), else of the side(s) holding material, else zero.
//  Returns true for a contact face (Riemann solution kept).
inline bool voidFace
(
    const vector& n,
    const vector& vL, const symmTensor& sL, const scalar zpL, const scalar zsL,
    const bool hasL, const bool touchL,
    const vector& vR, const symmTensor& sR, const scalar zpR, const scalar zsR,
    const bool hasR, const bool touchR,
    const bool noTension,
    vector& vS,
    vector& tS
)
{
    if (touchL && touchR)
    {
        riemann(n, vL, sL, zpL, zsL, vR, sR, zpR, zsR, vS, tS);
        if (!noTension || (tS & n) <= 0)
        {
            return true;
        }
    }
    tS = Zero;
    if (touchL)
    {
        vS = oneSided(n, 1, vL, sL, zpL, zsL);
    }
    else if (touchR)
    {
        vS = oneSided(n, -1, vR, sR, zpR, zsR);
    }
    else if (hasL && hasR)
    {
        vS =
            0.5
           *(
                oneSided(n, 1, vL, sL, zpL, zsL)
              + oneSided(n, -1, vR, sR, zpR, zsR)
            );
    }
    else if (hasL)
    {
        vS = oneSided(n, 1, vL, sL, zpL, zsL);
    }
    else if (hasR)
    {
        vS = oneSided(n, -1, vR, sR, zpR, zsR);
    }
    else
    {
        vS = Zero;
    }
    return false;
}


//- Mohr-Coulomb yield function (s1 - s3)/2 + (s1 + s3)/2 sin(phi)
//  - c cos(phi) of a stress (tension positive)
inline scalar mcYield(const symmTensor& s, const scalar sinPhi, const scalar cCos)
{
    const vector ev(eigenValues(s));
    const scalar s1 = cmptMax(ev);
    const scalar s3 = cmptMin(ev);
    return 0.5*(s1 - s3) + 0.5*(s1 + s3)*sinPhi - cCos;
}


//- Rule 8 (decision 0023; C++ decision 0032): shear cap at a contact face
//  with unit normal n (L to R). The face stress is the mean of the two
//  cells' stresses sMean with its n row replaced by the Riemann traction
//  tS. The shear jump term zsL zsR (vtR - vtL)/(zsL + zsR) is scaled back
//  (20 bisection steps) until the yield function is no larger than
//  without it (or zero). A capped face slips: L sees vtL + (tt - ttL)/zsL,
//  R sees vtR - (tt - ttR)/zsR. Returns true if capped, with tS changed
//  and the tangential velocities dL, dR that L and R see, relative to vS.
inline bool capShear
(
    const vector& n,
    const vector& vL, const symmTensor& sL, const scalar zsL,
    const vector& vR, const symmTensor& sR, const scalar zsR,
    const symmTensor& sMean,
    const scalar sinPhi, const scalar cCos,
    const vector& vS,
    vector& tS,
    vector& dL,
    vector& dR
)
{
    if (zsL <= 0 || zsR <= 0)
    {
        return false;
    }
    const vector tL(sL & n);
    const vector tR(sR & n);
    const vector ttL(tL - (tL & n)*n);
    const vector ttR(tR - (tR & n)*n);
    const vector vtL(vL - (vL & n)*n);
    const vector vtR(vR - (vR & n)*n);
    const scalar zs = zsL + zsR;
    const vector tb((zsR*ttL + zsL*ttR)/zs);
    const vector tj(zsL*zsR*(vtR - vtL)/zs);
    const scalar tn = tS & n;

    // face stress: tangential block of sMean, n row from the traction
    const tensor P(I - sqr(n));
    const tensor A(P & tensor(sMean) & P);
    auto faceYield = [&](const vector& tt)
    {
        return mcYield
        (
            symm(A + tn*sqr(n) + (n*tt) + (tt*n)),
            sinPhi,
            cCos
        );
    };

    const scalar lim = max(faceYield(tb), scalar(0));
    if (faceYield(tb + tj) <= lim)
    {
        return false;
    }
    scalar lo = 0, hi = 1;
    for (label k = 0; k < 20; ++k)
    {
        const scalar mid = 0.5*(lo + hi);
        if (faceYield(tb + mid*tj) <= lim)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    const vector tt(tb + lo*tj);
    tS = tn*n + tt;
    const vector vtS(vS - (vS & n)*n);
    dL = vtL + (tt - ttL)/zsL - vtS;
    dR = vtR - (tt - ttR)/zsR - vtS;
    return true;
}

} // End namespace lagrangeStepDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eulerianSolid::lagrangeStep(const scalar dt)
{
    using namespace lagrangeStepDetail;

    const fvMesh& mesh = mesh_;
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField& V = mesh.V().field();

    const scalar M = material_->waveModulus();
    const scalar G = material_->shearModulus();
    const scalar K = M - 4.0*G/3.0;

    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();

    // Cell mass, constant in the step
    const scalarField m(rho_.primitiveField()*V);

    // Gravity at mid-step (the driver advances the time before the step)
    const vector g(gravityAt(mesh.time().value() - 0.5*dt));

    // Impedances, with neighbour values on coupled patches
    tmp<volScalarField> tzp = workField<scalar>(mesh, "zp", dimless);
    tmp<volScalarField> tzs = workField<scalar>(mesh, "zs", dimless);
    volScalarField& zp = tzp.ref();
    volScalarField& zs = tzs.ref();
    zp.primitiveFieldRef() = sqrt(M*rho_.primitiveField());
    zs.primitiveFieldRef() = sqrt(G*rho_.primitiveField());

    // Free surface (decision 0031): homogenised partly filled cells
    // (decision 0021 rule 3): traction alpha sigma, impedance alpha z of
    // the material (density rho/alpha, bounded), mass rho V. The face
    // states and forces use seff = alpha sigma.
    const bool fs = freeSurface_;
    List<boundaryKind> kinds(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        kinds[patchi] =
            kindOf(mesh.boundary()[patchi], U_.boundaryField()[patchi], wallFriction_);
    }
    tmp<volScalarField> tsmooth, teFac;
    tmp<volSymmTensorField> tseff;
    if (fs)
    {
        alpha_.correctBoundaryConditions();
        const scalarField& a = alpha_.primitiveField();
        const scalarField& r = rho_.primitiveField();
        forAll(a, celli)
        {
            if (a[celli] > alphaEmpty)
            {
                const scalar rm =
                    min(max(r[celli]/a[celli], 0.5*rho0_), 2*rho0_);
                zp[celli] = a[celli]*sqrt(M*rm);
                zs[celli] = a[celli]*sqrt(G*rm);
            }
            else
            {
                zp[celli] = 0;
                zs[celli] = 0;
            }
        }

        label gAxis = -1;
        if (spreadLayers_ && mag(g_) > 0)
        {
            gAxis = label(faceAxis(g_));
        }
        slabSides(mesh, alpha_, kinds, lean_, smooth_, gAxis);

        tsmooth = workField<scalar>(mesh, "smooth", dimless);
        tsmooth.ref().primitiveFieldRef() = smooth_;
        tsmooth.ref().correctBoundaryConditions();

        // copy of sigma (zero-gradient patches, as sigma), times alpha
        tseff.reset
        (
            new volSymmTensorField
            (
                IOobject
                (
                    "seff",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                sigma_
            )
        );
        tseff.ref().primitiveFieldRef() = a*sigma_.primitiveField();
        tseff.ref().correctBoundaryConditions();

        // Rule 7 of decision 0021: a first-order cell keeps 1 - phi of its
        // hydrostatic jump, phi = 1 - |s_gg|/(w/2) in [0, 1], w its own
        // weight per area along g (phi 0: supported, at rest; 1: free fall)
        teFac = workField<scalar>(mesh, "eFac", dimless);
        scalarField& f = teFac.ref().primitiveFieldRef();
        f = 1;
        const scalar gm = mag(g);
        if (gm > 0)
        {
            const vector gHat(g/gm);
            const symmTensorField& se = tseff().primitiveField();
            forAll(f, celli)
            {
                const scalar w =
                    gm*r[celli]*(cmptMag(gHat) & hCell_[celli]);
                if (w > 0)
                {
                    const scalar sgg = gHat & (se[celli] & gHat);
                    f[celli] = 1 - min(max(1 - mag(sgg)/(0.5*w), 0), 1);
                }
            }
        }
        teFac.ref().correctBoundaryConditions();
    }
    zp.correctBoundaryConditions();
    zs.correctBoundaryConditions();
    const volSymmTensorField& sref = fs ? tseff() : sigma_;
    const volScalarField* smoothPtr = fs ? &tsmooth() : nullptr;
    const volScalarField* eFacPtr = fs ? &teFac() : nullptr;

    // Neighbour alpha and slab sides across coupled patches
    List<scalarField> aNbr(mesh.boundary().size());
    List<vectorField> leanNbr(mesh.boundary().size());
    if (fs)
    {
        tmp<volVectorField> tl = workField<vector>(mesh, "lean", dimless);
        tl.ref().primitiveFieldRef() = lean_;
        tl.ref().correctBoundaryConditions();
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                aNbr[patchi] =
                    alpha_.boundaryField()[patchi].patchNeighbourField();
                leanNbr[patchi] =
                    tl().boundaryField()[patchi].patchNeighbourField();
            }
        }
    }


    // 1. Face states ---------------------------------------------------------

    auto UO = surfaceVectorField::New
    (
        "UO", IOobject::NO_REGISTER, mesh, dimensionedVector(dimVelocity, Zero)
    );
    auto UN = surfaceVectorField::New
    (
        "UN", IOobject::NO_REGISTER, mesh, dimensionedVector(dimVelocity, Zero)
    );
    auto sO = surfaceSymmTensorField::New
    (
        "sO", IOobject::NO_REGISTER, mesh, dimensionedSymmTensor(dimPressure, Zero)
    );
    auto sN = surfaceSymmTensorField::New
    (
        "sN", IOobject::NO_REGISTER, mesh, dimensionedSymmTensor(dimPressure, Zero)
    );

    faceStates
    (
        U_, secondOrder_, UO.ref(), UN.ref(), nullptr, Zero, smoothPtr
    );
    faceStates
    (
        sref, secondOrder_, sO.ref(), sN.ref(), &rho_, g, smoothPtr, eFacPtr
    );

    if (secondOrder_)
    {
        // Hancock predictor from each cell's own face states
        tmp<volVectorField> tdU = workField<vector>(mesh, "dU", dimless);
        tmp<volSymmTensorField> tdS = workField<symmTensor>(mesh, "dS", dimless);
        volVectorField& dU = tdU.ref();
        volSymmTensorField& dS = tdS.ref();

        vectorField divS(mesh.nCells(), Zero);
        tensorField gradU(mesh.nCells(), Zero);

        forAll(own, facei)
        {
            const label P = own[facei];
            const label N = nei[facei];
            const vector& S = Sf[facei];
            divS[P] += (sO()[facei] - sref[P]) & S;
            divS[N] -= (sN()[facei] - sref[N]) & S;
            gradU[P] += S*(UO()[facei] - U_[P]);
            gradU[N] -= S*(UN()[facei] - U_[N]);
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const labelUList& fc = p.faceCells();
            const vectorField& pS = Sf.boundaryField()[patchi];
            const symmTensorField& ps = sO().boundaryField()[patchi];
            const vectorField& pu = UO().boundaryField()[patchi];
            if (kindOf(p, U_.boundaryField()[patchi], wallFriction_)
             == boundaryKind::wall)
            {
                // The predictor sees the wall's reaction: wall Riemann
                // solution from the face states (decision 0030)
                const scalar mu = wallFriction_[p.name()];
                const scalarField& pmS = magSf.boundaryField()[patchi];
                forAll(fc, i)
                {
                    const label P = fc[i];
                    vector v(U_[P]), t(Zero);
                    if (!fs || touches(alpha_[P], lean_[P], pS[i]))
                    {
                        wallFace
                        (
                            pS[i]/pmS[i], pu[i], ps[i], zp[P], zs[P], mu, v, t
                        );
                    }
                    divS[P] += t*pmS[i] - (sref[P] & pS[i]);
                    gradU[P] += pS[i]*(v - U_[P]);
                }
                continue;
            }
            forAll(fc, i)
            {
                divS[fc[i]] += (ps[i] - sref[fc[i]]) & pS[i];
                gradU[fc[i]] += pS[i]*(pu[i] - U_[fc[i]]);
            }
        }

        forAll(V, celli)
        {
            if (fs && smooth_[celli] < 0.5)
            {
                // First-order cell (decision 0021): no stress predictor;
                // empty cells do not move
                if (solveMomentum_ && m[celli] > 0)
                {
                    dU[celli] = 0.5*dt*(divS[celli]/m[celli] + g);
                }
                continue;
            }
            if (solveMomentum_)
            {
                dU[celli] = 0.5*dt*(divS[celli]/m[celli] + g);
            }
            const symmTensor D(symm(gradU[celli])/V[celli]);
            dS[celli] = 0.5*dt*(K*tr(D)*symmTensor::I + 2*G*dev(D));
        }
        dU.correctBoundaryConditions();
        dS.correctBoundaryConditions();

        surfaceVectorField& uo = UO.ref();
        surfaceVectorField& un = UN.ref();
        surfaceSymmTensorField& so = sO.ref();
        surfaceSymmTensorField& sn = sN.ref();
        forAll(own, facei)
        {
            uo[facei] += dU[own[facei]];
            un[facei] += dU[nei[facei]];
            so[facei] += dS[own[facei]];
            sn[facei] += dS[nei[facei]];
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const labelUList& fc = p.faceCells();
            vectorField& puo = uo.boundaryFieldRef()[patchi];
            vectorField& pun = un.boundaryFieldRef()[patchi];
            symmTensorField& pso = so.boundaryFieldRef()[patchi];
            symmTensorField& psn = sn.boundaryFieldRef()[patchi];
            const vectorField dUn
            (
                p.coupled()
              ? dU.boundaryField()[patchi].patchNeighbourField()
              : dU.boundaryField()[patchi].patchInternalField()
            );
            const symmTensorField dSn
            (
                p.coupled()
              ? dS.boundaryField()[patchi].patchNeighbourField()
              : dS.boundaryField()[patchi].patchInternalField()
            );
            forAll(fc, i)
            {
                puo[i] += dU[fc[i]];
                pso[i] += dS[fc[i]];
                pun[i] += dUn[i];
                psn[i] += dSn[i];
            }
        }
    }


    // 2. and 3. Riemann solver, face velocity, force and phiL ------------------

    auto tvS = surfaceVectorField::New
    (
        "vStar", IOobject::NO_REGISTER, mesh, dimensionedVector(dimVelocity, Zero)
    );
    auto tF = surfaceVectorField::New
    (
        "Fstar", IOobject::NO_REGISTER, mesh, dimensionedVector(dimForce, Zero)
    );
    surfaceVectorField& vS = tvS.ref();
    surfaceVectorField& F = tF.ref();
    nCapped_ = 0;

    // Shear cap (rule 8): tangential velocity seen by the owner and the
    // neighbour cell of a capped face, relative to vStar (zero elsewhere)
    tmp<surfaceVectorField> tdO, tdN;
    if (shearCap_)
    {
        tdO = surfaceVectorField::New
        (
            "slipO", IOobject::NO_REGISTER, mesh,
            dimensionedVector(dimVelocity, Zero)
        );
        tdN = surfaceVectorField::New
        (
            "slipN", IOobject::NO_REGISTER, mesh,
            dimensionedVector(dimVelocity, Zero)
        );
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector n(Sf[facei]/magSf[facei]);
        vector t;
        bool contact = true;
        if (fs)
        {
            const scalar aP = alpha_[P], aN = alpha_[N];
            contact = voidFace
            (
                n,
                UO()[facei], sO()[facei], zp[P], zs[P], aP > alphaEmpty,
                touches(aP, lean_[P], Sf[facei]),
                UN()[facei], sN()[facei], zp[N], zs[N], aN > alphaEmpty,
                touches(aN, lean_[N], -Sf[facei]),
                noTension_, vS[facei], t
            );
        }
        else
        {
            riemann
            (
                n,
                UO()[facei], sO()[facei], zp[P], zs[P],
                UN()[facei], sN()[facei], zp[N], zs[N],
                vS[facei], t
            );
        }
        if (shearCap_ && contact && solveMomentum_)
        {
            nCapped_ += capShear
            (
                n,
                UO()[facei], sO()[facei], zs[P],
                UN()[facei], sN()[facei], zs[N],
                0.5*(sref[P] + sref[N]),
                capSinPhi_, capCCos_, vS[facei], t,
                tdO.ref()[facei], tdN.ref()[facei]
            );
        }
        if (!solveMomentum_)
        {
            // Prescribed motion: the face moves with the material on both
            // sides, no acoustic correction from stress jumps (decision 0022)
            vS[facei] = 0.5*(UO()[facei] + UN()[facei]);
        }
        F[facei] = t*magSf[facei];
        phiL_[facei] = vS[facei] & Sf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const fvPatchVectorField& pU = U_.boundaryField()[patchi];
        const boundaryKind kind = kindOf(p, pU, wallFriction_);
        if (kind == boundaryKind::empty)
        {
            continue;
        }

        const labelUList& fc = p.faceCells();
        const vectorField& pS = Sf.boundaryField()[patchi];
        const scalarField& pmS = magSf.boundaryField()[patchi];
        const vectorField& puo = UO().boundaryField()[patchi];
        const symmTensorField& pso = sO().boundaryField()[patchi];
        vectorField& pvS = vS.boundaryFieldRef()[patchi];
        vectorField& pF = F.boundaryFieldRef()[patchi];
        scalarField& pphi = phiL_.boundaryFieldRef()[patchi];

        if (kind == boundaryKind::coupled)
        {
            // Solved on the owner side; syncCoupled copies to the other
            const bool owner =
                isA<cyclicFvPatch>(p)
              ? refCast<const cyclicFvPatch>(p).owner()
              : refCast<const processorFvPatch>(p).procPolyPatch().owner();
            if (!owner)
            {
                continue;
            }
            const vectorField& pun = UN().boundaryField()[patchi];
            const symmTensorField& psn = sN().boundaryField()[patchi];
            const scalarField zpn(zp.boundaryField()[patchi].patchNeighbourField());
            const scalarField zsn(zs.boundaryField()[patchi].patchNeighbourField());
            symmTensorField srn;
            if (shearCap_)
            {
                srn = sref.boundaryField()[patchi].patchNeighbourField();
            }
            forAll(fc, i)
            {
                const vector n(pS[i]/pmS[i]);
                vector t;
                bool contact = true;
                if (fs)
                {
                    const scalar aP = alpha_[fc[i]];
                    const scalar aN = aNbr[patchi][i];
                    contact = voidFace
                    (
                        n,
                        puo[i], pso[i], zp[fc[i]], zs[fc[i]], aP > alphaEmpty,
                        touches(aP, lean_[fc[i]], pS[i]),
                        pun[i], psn[i], zpn[i], zsn[i], aN > alphaEmpty,
                        touches(aN, leanNbr[patchi][i], -pS[i]),
                        noTension_, pvS[i], t
                    );
                }
                else
                {
                    riemann
                    (
                        n,
                        puo[i], pso[i], zp[fc[i]], zs[fc[i]],
                        pun[i], psn[i], zpn[i], zsn[i],
                        pvS[i], t
                    );
                }
                if (shearCap_ && contact && solveMomentum_)
                {
                    nCapped_ += capShear
                    (
                        n,
                        puo[i], pso[i], zs[fc[i]],
                        pun[i], psn[i], zsn[i],
                        0.5*(sref[fc[i]] + srn[i]),
                        capSinPhi_, capCCos_, pvS[i], t,
                        tdO.ref().boundaryFieldRef()[patchi][i],
                        tdN.ref().boundaryFieldRef()[patchi][i]
                    );
                }
                if (!solveMomentum_)
                {
                    pvS[i] = 0.5*(puo[i] + pun[i]);
                }
                pF[i] = t*pmS[i];
                pphi[i] = pvS[i] & pS[i];
            }
            continue;
        }

        forAll(fc, i)
        {
            const label P = fc[i];
            const vector n(pS[i]/pmS[i]);
            const vector& vL = puo[i];
            const vector tL(pso[i] & n);
            const scalar vnL = vL & n;
            const scalar tnL = tL & n;
            const vector vtL(vL - vnL*n);
            const vector ttL(tL - tnL*n);

            vector v, t;
            if (fs && !touches(alpha_[P], lean_[P], pS[i]))
            {
                // Material not at the face (decision 0021 rule 3)
                v = Zero;
                t = Zero;
            }
            else if (kind == boundaryKind::velocity)
            {
                // Prescribed velocity; traction from the interior
                // characteristic s* - z v* = sL - z vL
                v = pU[i];
                const scalar vn = v & n;
                const vector vt(v - vn*n);
                t = (tnL + zp[P]*(vn - vnL))*n + ttL + zs[P]*(vt - vtL);
            }
            else if (kind == boundaryKind::wall)
            {
                wallFace
                (
                    n, vL, pso[i], zp[P], zs[P], wallFriction_[p.name()],
                    v, t
                );
            }
            else if (kind == boundaryKind::symmetry)
            {
                // No normal velocity, no shear traction
                v = vtL - ttL/zs[P];
                t = (tnL - zp[P]*vnL)*n;
            }
            else
            {
                // Free surface: zero traction
                v = (vnL - tnL/zp[P])*n + vtL - ttL/zs[P];
                t = Zero;
            }
            pvS[i] = v;
            pF[i] = t*pmS[i];
            pphi[i] = v & pS[i];
        }
    }

    if (shearCap_)
    {
        syncCoupled(mesh, vS, F, phiL_, &tdO.ref(), &tdN.ref());
    }
    else
    {
        syncCoupled(mesh, vS, F, phiL_);
    }
    vStar_ = vS;
    Fstar_ = F;


    // 4. Momentum, velocity gradient and Lagrangian volume ----------------------

    vectorField force(mesh.nCells(), Zero);
    tensorField T(mesh.nCells(), Zero);      // sum Sf * vStar
    scalarField dV(mesh.nCells(), Zero);     // sum phiL out of the cell
    vectorField dVa;                         // the same per face axis (fs)
    if (fs)
    {
        dVa.setSize(mesh.nCells(), Zero);
        forAll(own, facei)
        {
            const direction d = faceAxis(Sf[facei]);
            dVa[own[facei]][d] += phiL_[facei];
            dVa[nei[facei]][d] -= phiL_[facei];
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const labelUList& fc = p.faceCells();
            const vectorField& pS = Sf.boundaryField()[patchi];
            const scalarField& pphi = phiL_.boundaryField()[patchi];
            forAll(fc, i)
            {
                dVa[fc[i]][faceAxis(pS[i])] += pphi[i];
            }
        }
    }

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        force[P] += F[facei];
        force[N] -= F[facei];
        if (shearCap_)
        {
            // a capped face slips: each side sees its own velocity
            T[P] += Sf[facei]*(vS[facei] + tdO()[facei]);
            T[N] -= Sf[facei]*(vS[facei] + tdN()[facei]);
        }
        else
        {
            const tensor t(Sf[facei]*vS[facei]);
            T[P] += t;
            T[N] -= t;
        }
        dV[P] += phiL_[facei];
        dV[N] -= phiL_[facei];
    }
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pS = Sf.boundaryField()[patchi];
        const vectorField& pF = F.boundaryField()[patchi];
        const vectorField& pvS = vS.boundaryField()[patchi];
        const scalarField& pphi = phiL_.boundaryField()[patchi];
        forAll(fc, i)
        {
            force[fc[i]] += pF[i];
            T[fc[i]] += pS[i]*pvS[i];
            dV[fc[i]] += pphi[i];
        }
        if (shearCap_)
        {
            const vectorField& pdO = tdO().boundaryField()[patchi];
            forAll(fc, i)
            {
                T[fc[i]] += pS[i]*pdO[i];
            }
        }
    }

    if (solveMomentum_)
    {
        vectorField& U = U_.primitiveFieldRef();
        forAll(U, celli)
        {
            if (fs && m[celli] <= 0)
            {
                U[celli] = Zero;
                continue;
            }
            U[celli] += dt*(force[celli]/m[celli] + g);
            if (damping_ > 0)
            {
                U[celli] /= 1 + damping_*dt;
            }
        }
    }


    // 5. Rotation and constitutive update --------------------------------------

    symmTensorField& sigma = sigma_.primitiveFieldRef();
    symmTensorField& epsP = epsP_.primitiveFieldRef();
    scalarField& epsPEq = epsPEq_.primitiveFieldRef();
    scalarField& rho = rho_.primitiveFieldRef();

    forAll(sigma, celli)
    {
        const scalar VL = V[celli] + dt*dV[celli];
        const scalar Vmid = 0.5*(V[celli] + VL);

        // T/V = grad(v) in OpenFOAM order (d_i v_j); L = grad(v)^T
        const tensor L(dt*T[celli].T()/Vmid);
        symmTensor dEps(symm(L));
        const tensor dW(skew(L));

        if (fs)
        {
            // Void closure (decision 0021 rule 5), per axis of the
            // Cartesian cell: compression closes the void first; expansion
            // first unloads the material to zero normal stress (no cap for
            // a full cell of a material with tensile strength), the rest
            // opens void. Shear goes to the material. A full cell that
            // opens no void keeps the strain of the full scheme.
            const scalar a = alpha_[celli];
            if (a <= alphaEmpty)
            {
                sigma[celli] = Zero;
                rho[celli] = 0;
                alphaL_[celli] = 0;
                continue;
            }
            const bool full = a >= alphaFull;
            const Vector<label>& geomD = mesh.geometricD();
            static const direction diag[3] =
                {symmTensor::XX, symmTensor::YY, symmTensor::ZZ};
            bool same = full;
            scalar matV = full ? 1 : a;     // material volume / V
            symmTensor dE(dEps);
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                if (geomD[d] < 0)
                {
                    continue;
                }
                const scalar r = dt*dVa[celli][d]/V[celli];
                scalar dm;
                if (r < 0)
                {
                    dm = full ? r : min(scalar(0), r + (1 - a));
                }
                else
                {
                    const scalar cap =
                        (full && !noTension_)
                      ? GREAT
                      : a*max(scalar(0), -sigma[celli].component(diag[d]))/M;
                    dm = min(r, cap);
                }
                if (dm != r)
                {
                    same = false;
                }
                matV += dm;
                dE.component(diag[d]) = dm/(a + 0.5*dm);
            }
            if (!same)
            {
                dEps = dE;
                alphaL_[celli] =
                    min(max(matV*V[celli]/VL, scalar(0)), scalar(1));
            }
            else
            {
                alphaL_[celli] = 1;
            }
        }

        if (dW != tensor::zero)
        {
            const tensor Q(hughesWinget(dW));
            sigma[celli] = transform(Q, sigma[celli]);
            epsP[celli] = transform(Q, epsP[celli]);
        }

        if (dEps != symmTensor::zero)
        {
            constitutiveModel::history hOld, hNew;
            constitutiveModel::toVoigt(epsP[celli], hOld.epsP, 2);
            hOld.epsPEq = epsPEq[celli];

            symmTensor sNew;
            material_->update(dEps, sigma[celli], hOld, sNew, hNew);

            sigma[celli] = sNew;
            symmTensor e(constitutiveModel::fromVoigt(hNew.epsP));
            e.xy() *= 0.5;
            e.yz() *= 0.5;
            e.xz() *= 0.5;
            epsP[celli] = e;
            epsPEq[celli] = hNew.epsPEq;
        }

        // 6. Lagrangian density
        rho[celli] = m[celli]/VL;
    }

    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();
    epsP_.correctBoundaryConditions();
    epsPEq_.correctBoundaryConditions();
}


// ************************************************************************* //
