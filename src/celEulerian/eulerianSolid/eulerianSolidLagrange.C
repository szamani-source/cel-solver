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
       Boundary faces from the patch type of U (decision 0018).
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
template<class Type>
void faceStates
(
    const GeometricField<Type, fvPatchField, volMesh>& q,
    const bool secondOrder,
    GeometricField<Type, fvsPatchField, surfaceMesh>& qO,
    GeometricField<Type, fvsPatchField, surfaceMesh>& qN
)
{
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
    auto states = [&]
    (
        const Type& qP, const Type& qQ,
        const UList<vector>* gP, const UList<vector>* gQ,
        const vector& dPQ,
        Type& outP, Type& outQ
    )
    {
        outP = qP;
        outQ = qQ;
        if (!secondOrder)
        {
            return;
        }
        for (direction d = 0; d < nCmpt; ++d)
        {
            const scalar a = component(qQ, d) - component(qP, d);
            const scalar sP = mcSlope(a, dPQ & (*gP)[d]);
            const scalar sQ = mcSlope(a, dPQ & (*gQ)[d]);
            setComponent(outP, d) = component(qP, d) + 0.5*sP;
            setComponent(outQ, d) = component(qQ, d) - 0.5*sQ;
        }
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
            q[P], q[N], &gP, &gQ, delta[facei], qO[facei], qN[facei]
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
                states(q[fc[i]], qn[i], &gP, &gQ, pd[i], pO[i], pN[i]);
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
                Type ghost;
                states
                (
                    q[fc[i]], 2*qb[i] - q[fc[i]], &gP, &gP,
                    2*(Cf[i] - C[fc[i]]), pO[i], ghost
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
enum class boundaryKind { coupled, empty, velocity, symmetry, traction };


boundaryKind kindOf(const fvPatch& p, const fvPatchVectorField& pU)
{
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


//- Copy the owner-side coupled face values to the other side:
//  vStar as is, F and phi negated. Cyclic locally, processor by message.
void syncCoupled
(
    const fvMesh& mesh,
    surfaceVectorField& vS,
    surfaceVectorField& F,
    surfaceScalarField& phi
)
{
    auto& bv = vS.boundaryFieldRef();
    auto& bF = F.boundaryFieldRef();
    auto& bp = phi.boundaryFieldRef();

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

        std::map<label, std::tuple<vectorField, vectorField, scalarField>> got;
        for (label k = 0; k < n; ++k)
        {
            label tag;
            vectorField v, f;
            scalarField ph;
            is >> tag >> v >> f >> ph;
            got[tag] = std::make_tuple(v, f, ph);
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
        }
    }
}


//- Hughes-Winget rotation for the spin increment dW:
//  Q = (I - dW/2)^-1 (I + dW/2)
inline tensor hughesWinget(const tensor& dW)
{
    return inv(I - 0.5*dW) & (I + 0.5*dW);
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

    // Impedances, with neighbour values on coupled patches
    tmp<volScalarField> tzp = workField<scalar>(mesh, "zp", dimless);
    tmp<volScalarField> tzs = workField<scalar>(mesh, "zs", dimless);
    volScalarField& zp = tzp.ref();
    volScalarField& zs = tzs.ref();
    zp.primitiveFieldRef() = sqrt(M*rho_.primitiveField());
    zs.primitiveFieldRef() = sqrt(G*rho_.primitiveField());
    zp.correctBoundaryConditions();
    zs.correctBoundaryConditions();


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

    faceStates(U_, secondOrder_, UO.ref(), UN.ref());
    faceStates(sigma_, secondOrder_, sO.ref(), sN.ref());

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
            divS[P] += (sO()[facei] - sigma_[P]) & S;
            divS[N] -= (sN()[facei] - sigma_[N]) & S;
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
            forAll(fc, i)
            {
                divS[fc[i]] += (ps[i] - sigma_[fc[i]]) & pS[i];
                gradU[fc[i]] += pS[i]*(pu[i] - U_[fc[i]]);
            }
        }

        forAll(V, celli)
        {
            if (solveMomentum_)
            {
                dU[celli] = 0.5*dt*divS[celli]/m[celli];
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

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const vector n(Sf[facei]/magSf[facei]);
        vector t;
        riemann
        (
            n,
            UO()[facei], sO()[facei], zp[P], zs[P],
            UN()[facei], sN()[facei], zp[N], zs[N],
            vS[facei], t
        );
        F[facei] = t*magSf[facei];
        phiL_[facei] = vS[facei] & Sf[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        const fvPatchVectorField& pU = U_.boundaryField()[patchi];
        const boundaryKind kind = kindOf(p, pU);
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
            forAll(fc, i)
            {
                const vector n(pS[i]/pmS[i]);
                vector t;
                riemann
                (
                    n,
                    puo[i], pso[i], zp[fc[i]], zs[fc[i]],
                    pun[i], psn[i], zpn[i], zsn[i],
                    pvS[i], t
                );
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
            if (kind == boundaryKind::velocity)
            {
                // Prescribed velocity; traction from the interior
                // characteristic s* - z v* = sL - z vL
                v = pU[i];
                const scalar vn = v & n;
                const vector vt(v - vn*n);
                t = (tnL + zp[P]*(vn - vnL))*n + ttL + zs[P]*(vt - vtL);
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

    syncCoupled(mesh, vS, F, phiL_);


    // 4. Momentum, velocity gradient and Lagrangian volume ----------------------

    vectorField force(mesh.nCells(), Zero);
    tensorField T(mesh.nCells(), Zero);      // sum Sf * vStar
    scalarField dV(mesh.nCells(), Zero);     // sum phiL out of the cell

    forAll(own, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        force[P] += F[facei];
        force[N] -= F[facei];
        const tensor t(Sf[facei]*vS[facei]);
        T[P] += t;
        T[N] -= t;
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
    }

    if (solveMomentum_)
    {
        vectorField& U = U_.primitiveFieldRef();
        forAll(U, celli)
        {
            U[celli] += dt*force[celli]/m[celli];
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
        const symmTensor dEps(symm(L));
        const tensor dW(skew(L));

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
