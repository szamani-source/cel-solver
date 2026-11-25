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
    rest to round-off; the momentum and Hancock predictor add g. With g
    not along an axis the jump is that of a layer uniform across g's
    dominant axis (inclined layer, decision 0039). Wall
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

#include "cubicEqn.H"
#include <type_traits>
#include <algorithm>
#include "eulerianSolid.H"
#include "celNonlocal.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "processorPolyPatch.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "PstreamBuffers.H"
#include "calculatedFvPatchFields.H"
#include "equilibriumJump.H"
#include "surfaceInterpolate.H"

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


using eulerianDetail::equilibriumJump;


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


//- Axisymmetric wedge meshes (WP 3.8a, decision 0061): a copy of q whose
//  values on non-coupled, non-wedge boundary faces normal to gravity (all
//  of them without gravity), such as the side wall of a vertical cylinder,
//  are linear extrapolations of the cell value instead of the cell value.
//  The cell gradient is the Gauss gradient made consistent with those
//  extrapolated face values, g_a = g0_a/(1 - K_aa), K = sum S_b d_b / V
//  over the extrapolated faces (diagonal on axis-aligned faces). On a wedge the hoop stress of a cell beside such a
//  wall is what balances the radial face forces over one cell, so a first-
//  order wall traction (the cell value) leaves an O(1) hoop stress error
//  (Lame cylinder); extrapolated, it is O(h). Faces along gravity keep the
//  hydrostatic mirror (decision 0030).
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> extrapolatedWalls
(
    const GeometricField<Type, fvPatchField, volMesh>& q,
    const vector& gravity
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    const fvMesh& mesh = q.mesh();
    const direction nCmpt = pTraits<Type>::nComponents;
    const vector gHat(mag(gravity) > 0 ? gravity/mag(gravity) : vector::zero);

    wordList types(q.boundaryField().types());
    List<boolList> extrap(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (p.coupled() || isA<emptyFvPatch>(p) || isA<wedgeFvPatch>(p))
        {
            continue;
        }
        types[patchi] = fvPatchFieldBase::calculatedType();
        const vectorField nf(p.nf());
        extrap[patchi].setSize(p.size());
        forAll(nf, i)
        {
            extrap[patchi][i] = mag(nf[i] & gHat) < 1e-6;
        }
    }
    tmp<fieldType> tx
    (
        new fieldType
        (
            IOobject
            (
                q.name() + "Walls",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            q,
            types
        )
    );
    fieldType& x = tx.ref();

    PtrList<volVectorField> g0;
    componentGradients(q, g0);
    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V().field();

    tensorField K(mesh.nCells(), Zero);
    boolList touched(mesh.nCells(), false);
    forAll(mesh.boundary(), patchi)
    {
        if (extrap[patchi].empty()) continue;
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& Sf = p.Sf();
        const vectorField& Cf = p.Cf();
        forAll(fc, i)
        {
            if (!extrap[patchi][i]) continue;
            K[fc[i]] += (Sf[i]*(Cf[i] - C[fc[i]]))/V[fc[i]];
            touched[fc[i]] = true;
        }
    }
    forAll(mesh.boundary(), patchi)
    {
        if (extrap[patchi].empty()) continue;
        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& Cf = p.Cf();
        Field<Type>& xb = x.boundaryFieldRef()[patchi];
        forAll(fc, i)
        {
            const label P = fc[i];
            Type v(q[P]);
            if (extrap[patchi][i] && touched[P])
            {
                // Faces are along the mesh axes, so K is diagonal: along
                // axis a, g_a = g0_a/(1 - K_aa). A cell with extrapolated
                // faces on both sides of an axis (one cell across) has no
                // gradient there: it keeps the cell value
                const vector d(Cf[i] - C[P]);
                const vector Kd(K[P].diag());
                for (direction c = 0; c < nCmpt; ++c)
                {
                    scalar dq = 0;
                    for (direction a = 0; a < vector::nComponents; ++a)
                    {
                        if (1 - Kd[a] > 0.25)
                        {
                            dq += g0[c][P][a]/(1 - Kd[a])*d[a];
                        }
                    }
                    setComponent(v, c) = component(q[P], c) + dq;
                }
            }
            xb[i] = v;
        }
    }
    return tx;
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
    const volScalarField* eFac = nullptr,
    const List<tensor>* rot = nullptr
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
    auto statesCore = [&]
    (
        const Type& qP, const Type& qQ,
        const Type& DP, const Type& DQ,
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
                const scalar sP = mcSlope(a, component(DP, d));
                const scalar sQ = mcSlope(a, component(DQ, d));
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
            const scalar sP = mcSlope(a, (component(DP, d)) - ep);
            const scalar sQ = mcSlope(a, (component(DQ, d)) - eq);
            setComponent(outP, d) = component(qP, d) + 0.5*(sP + ep);
            setComponent(outQ, d) = component(qQ, d) - 0.5*(sQ + eq);
        }
    };


    // Turned cells (WP G.3): the face states are reconstructed and limited
    // in the mesh lines' frame R of the face (faceR, set per face), so a
    // turned mesh gives what a Cartesian one gives with gravity turned:
    // states, jumps and the directional derivatives along dPQ turned into
    // R, the usual component-wise rules, the result turned back
    const tensor* faceR = nullptr;
    auto rotT = [](const tensor& R, const Type& x) -> Type
    {
        if constexpr (std::is_same<Type, symmTensor>::value)
        {
            return transform(R, x);
        }
        else if constexpr (std::is_same<Type, vector>::value)
        {
            return R & x;
        }
        else
        {
            return x;
        }
    };
    auto states0 = [&]
    (
        const Type& qP, const Type& qQ,
        const UList<vector>* gP, const UList<vector>* gQ,
        const vector& dPQ,
        Type& outP, Type& outQ,
        const Type& eP, const Type& eQ
    )
    {
        Type DP(Zero), DQ(Zero);
        if (secondOrder)
        {
            for (direction d = 0; d < nCmpt; ++d)
            {
                setComponent(DP, d) = dPQ & (*gP)[d];
                setComponent(DQ, d) = dPQ & (*gQ)[d];
            }
        }
        if (!faceR)
        {
            statesCore(qP, qQ, DP, DQ, outP, outQ, eP, eQ);
            return;
        }
        const tensor& R = *faceR;
        Type oP, oQ;
        statesCore
        (
            rotT(R, qP), rotT(R, qQ), rotT(R, DP), rotT(R, DQ), oP, oQ,
            rotT(R, eP), rotT(R, eQ)
        );
        outP = rotT(R.T(), oP);
        outQ = rotT(R.T(), oQ);
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
        if (smP == 1 && smQ == 1)
        {
            states0(qP, qQ, gP, gQ, dPQ, outP, outQ, eP, eQ);
            return;
        }
        Type tP(qP), tQ(qQ);
        if (smP > 0 || smQ > 0)
        {
            states0(qP, qQ, gP, gQ, dPQ, tP, tQ, eP, eQ);
        }
        const Type foP(grav ? Type(qP + (0.5*fP)*eP) : qP);
        const Type foQ(grav ? Type(qQ - (0.5*fQ)*eQ) : qQ);
        // smooth is 0 or 1 without the blend (WP 3.6h): one of the two
        outP =
            smP == 1 ? tP : (smP == 0 ? foP : Type(smP*tP + (1 - smP)*foP));
        outQ =
            smQ == 1 ? tQ : (smQ == 0 ? foQ : Type(smQ*tQ + (1 - smQ)*foQ));
    };

    auto sm = [&](const label celli) -> scalar
    {
        return smooth ? (*smooth)[celli] : 1;
    };
    auto ef = [&](const label celli) -> scalar
    {
        return eFac ? (*eFac)[celli] : 1;
    };

    // In a turned cell (square mesh lines, not on the axes; WP G.3) the
    // jump is taken in the lines' frame R, as on a Cartesian mesh with
    // gravity turned, and turned back: R^T J(R d, R g) R
    auto eqJump = [&](const scalar r, const vector& d, const label c) -> Type
    {
        if (!grav)
        {
            return zeroT;
        }
        if constexpr (std::is_same<Type, symmTensor>::value)
        {
            if (rot && !rot->empty() && (*rot)[c] != tensor::I)
            {
                const tensor& R = (*rot)[c];
                return transform
                (
                    R.T(), equilibriumJump<Type>(r, R & d, R & gravity)
                );
            }
        }
        return equilibriumJump<Type>(r, d, gravity);
    };

    List<vector> gP(nCmpt), gQ(nCmpt);

    auto turnedR = [&](const label c) -> const tensor*
    {
        return rot && !rot->empty() && (*rot)[c] != tensor::I
            ? &(*rot)[c] : nullptr;
    };

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
        faceR = turnedR(P);
        if (!faceR)
        {
            faceR = turnedR(N);
        }
        states
        (
            q[P], q[N], &gP, &gQ, delta[facei], qO[facei], qN[facei],
            eqJump(grav ? (*rho)[P] : 0, delta[facei], P),
            eqJump(grav ? (*rho)[N] : 0, delta[facei], N),
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
                faceR = turnedR(fc[i]);
                states
                (
                    q[fc[i]], qn[i], &gP, &gQ, pd[i], pO[i], pN[i],
                    eqJump(grav ? (*rho)[fc[i]] : 0, pd[i], fc[i]),
                    eqJump(rn[i], pd[i], fc[i]),
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
                faceR = turnedR(fc[i]);
                const vector dG(2*(Cf[i] - C[fc[i]]));
                const Type e(eqJump(grav ? (*rho)[fc[i]] : 0, dG, fc[i]));
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
//  face velocity vS and traction tS = sigma* & n. kj >= 0 replaces the
//  shear jump coefficient zsL zsR/(zsL + zsR) of the traction (flow
//  impedance, WP 3.6c); the velocity keeps the elastic impedances.
inline void riemann
(
    const vector& n,
    const vector& vL, const symmTensor& sL, const scalar zpL, const scalar zsL,
    const vector& vR, const symmTensor& sR, const scalar zpR, const scalar zsR,
    vector& vS,
    vector& tS,
    const scalar kj = -1
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
    const vector tt
    (
        kj < 0
      ? (zsR*ttL + zsL*ttR + zsL*zsR*(vtR - vtL))/zs
      : (zsR*ttL + zsL*ttR)/zs + kj*(vtR - vtL)
    );

    vS = vn*n + vt;
    tS = tn*n + tt;
}


//- Boundary kinds (decision 0018)
enum class boundaryKind
{
    coupled, empty, velocity, symmetry, traction, wall, wedge
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
        return boundaryKind::wedge;     // WP 3.8a, decision 0061
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
    vector& v, vector& t,
    const scalar kj = -1,
    const scalar adh = 0
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
    vector tt(ttL - (kj < 0 ? zs : kj)*vtL);
    // Coulomb with adhesion (rough walls, WP 4.3e): |tt| <= adh - tn mu
    const scalar limit = adh - tn*mu;
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


//- Weight of the full rules for a cell of volume fraction a, with the
//  blend (WP 3.6h, decision 0049): 1 for void <= v1, 0 for void >= v2
//  (defaults 1e-4, 1e-2), linear in log(void) between. Without the blend
//  (v2 = 0): 1 if full, else 0.
inline scalar fullWeight(const scalar a, const scalar v1, const scalar v2)
{
    if (a >= alphaFull)
    {
        return 1;
    }
    if (!(v2 > 0))
    {
        return 0;
    }
    const scalar v = 1 - a;
    if (v <= v1)
    {
        return 1;
    }
    if (v >= v2)
    {
        return 0;
    }
    return Foam::log(v2/v)/Foam::log(v2/v1);
}


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
    const label gAxis = -1,
    const scalar v1 = 0,
    const scalar v2 = 0,
    const vector& gVec = Zero,
    const surfaceScalarField* faceOpen = nullptr,
    const List<tensor>* frame = nullptr,
    const vectorField* hLines = nullptr
)
{
    // Angled meshes (WP G.2, decision 0104): a face's area vector in the
    // cell's mesh-line frame; the vector itself on axis-aligned meshes
    auto L = [frame](const label c, const vector& S) -> vector
    {
        return frame ? ((*frame)[c] & S) : S;
    };
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
    auto set = [&](const label c, const vector& Sout0, const scalar ab)
    {
        const vector Sout(L(c, Sout0));
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
    // Immersed walls (WP 4.3c): a face closed by a body mirrors the cell,
    // as a wall patch
    auto closed = [&](const label patchi, const label i) -> bool
    {
        if (!faceOpen)
        {
            return false;
        }
        return
            patchi < 0
          ? (*faceOpen)[i] <= 0
          : faceOpen->boundaryField()[patchi][i] <= 0;
    };
    forAll(own, facei)
    {
        if (closed(-1, facei)) continue;
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
                if (closed(patchi, i)) continue;
                set(fc[i], pS[i], an[i]);
            }
        }
        else if (kinds[patchi] == boundaryKind::traction)
        {
            forAll(fc, i)
            {
                if (closed(patchi, i)) continue;
                set(fc[i], pS[i], 0);
            }
        }
    }

    // A wall above the cell (outward normal against gravity, a ceiling)
    // counts as empty for the slab side only (WP 3.6e, decision 0050):
    // mirrored, it made splashed material lean against the ceiling, never
    // touch the face below and fall in place forever. The smooth flag keeps
    // the mirror (a full cell under a lid stays second order).
    vectorField aPlusL(aPlus), aMinusL(aMinus);
    if (mag(gVec) > 0)
    {
        forAll(mesh.boundary(), patchi)
        {
            if
            (
                kinds[patchi] != boundaryKind::wall
             && kinds[patchi] != boundaryKind::velocity
            )
            {
                continue;
            }
            const fvPatch& p = mesh.boundary()[patchi];
            const labelUList& fc = p.faceCells();
            const vectorField& pS = p.Sf();
            forAll(fc, i)
            {
                if ((pS[i] & gVec) < -1e-6*mag(pS[i])*mag(gVec))
                {
                    const vector Sl(L(fc[i], pS[i]));
                    const direction d = faceAxis(Sl);
                    if (Sl[d] > 0)
                    {
                        aPlusL[fc[i]][d] = 0;
                    }
                    else
                    {
                        aMinusL[fc[i]][d] = 0;
                    }
                }
            }
        }
    }

    lean.setSize(mesh.nCells());
    smooth.setSize(mesh.nCells());
    forAll(a, celli)
    {
        // With the blend a nearly full cell (weight >= 1/2, void below
        // sqrt(v1 v2)) has an upright slab: it touches all its faces
        const bool blend = v2 > 0;
        const bool full =
            blend
          ? fullWeight(a[celli], v1, v2) >= 0.5
          : a[celli] >= alphaFull;
        const vector dif(aPlusL[celli] - aMinusL[celli]);
        for (direction d = 0; d < vector::nComponents; ++d)
        {
            lean[celli][d] =
                (full || mag(dif[d]) <= 1e-9) ? 0 : sign(dif[d]);
        }
        if (blend)
        {
            // Smallest weight of the cell and its neighbours
            scalar w = fullWeight(a[celli], v1, v2);
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                w = min(w, fullWeight(aPlus[celli][d], v1, v2));
                w = min(w, fullWeight(aMinus[celli][d], v1, v2));
            }
            smooth[celli] = w;
        }
        else
        {
            smooth[celli] =
                (
                    full
                 && cmptMin(aPlus[celli]) >= alphaFull
                 && cmptMin(aMinus[celli]) >= alphaFull
                )
              ? 1 : 0;
        }
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
    if (frame)
    {
        // Angled cells (WP G.2, decision 0104): the differences along the
        // mesh lines, D_a ~ 2 h_a t_a.grad(alpha), become the gradient's
        // x, y, z components (times 2 h, h the mean line length), so the
        // rule compares the gradient along gravity with the horizontal
        // one, as on an axis-aligned mesh. Along the mesh lines a leaning
        // column mixed the surface slope into the gradient along gravity
        vectorField& Dp = D.primitiveFieldRef();
        forAll(Dp, celli)
        {
            const vector& h = (*hLines)[celli];
            const scalar hm = cmptAv(h);
            vector q(Dp[celli]);
            for (direction k = 0; k < vector::nComponents; ++k)
            {
                q[k] *= h[k] > VSMALL ? hm/h[k] : 0;
            }
            Dp[celli] = inv((*frame)[celli]) & q;
        }
    }
    D.correctBoundaryConditions();
    const vectorField& Dc = D.primitiveField();

    // Neighbour D along each axis and side: index 2 e + (0 for +, 1 for -)
    List<FixedList<vector, 6>> Dn(mesh.nCells());
    forAll(Dn, celli)
    {
        Dn[celli] = Dc[celli];
    }
    auto setN = [&](const label c, const vector& Sout0, const vector& Db)
    {
        const vector Sout(L(c, Sout0));
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
//  sides touch) gets the Riemann solution unless its normal traction pulls
//  harder than the tensile strength sT (0: no tensile strength, decision
//  0054); every other face is
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
    const scalar sT,
    vector& vS,
    vector& tS,
    const scalar kj = -1
)
{
    if (touchL && touchR)
    {
        riemann(n, vL, sL, zpL, zsL, vR, sR, zpR, zsR, vS, tS, kj);
        if ((tS & n) <= sT)
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


//- Principal stresses, ascending: OpenFOAM's eigenValues (v2512,
//  symmTensor.cxx) step for step, so every result is unchanged bit for bit,
//  without its warning. On wedge stresses whose xz, yz shears are round-off
//  the cubic solver reports a non-real root and printed "Eigenvalue
//  computation fails" for each: millions of lines per Couto silo run, which
//  filled the disk and slowed the runs. Such a root is used as OpenFOAM
//  uses it.
inline vector principal(const symmTensor& T)
{
    if ((sqr(T.xy()) + sqr(T.xz()) + sqr(T.yz())) < ROOTSMALL)
    {
        vector eVals(T.diag());
        std::sort(eVals.begin(), eVals.end());
        return eVals;
    }
    const scalar b = - T.xx() - T.yy() - T.zz();
    const scalar c =
        T.xx()*T.yy() + T.xx()*T.zz() + T.yy()*T.zz()
      - T.xy()*T.yx() - T.yz()*T.zy() - T.zx()*T.xz();
    const scalar d =
      - T.xx()*T.yy()*T.zz()
      - T.xy()*T.yz()*T.zx() - T.xz()*T.zy()*T.yx()
      + T.xx()*T.yz()*T.zy() + T.yy()*T.zx()*T.xz() + T.zz()*T.xy()*T.yx();
    const Roots<3> roots(cubicEqn(1, b, c, d).roots());
    vector eVals(roots[0], roots[1], roots[2]);
    std::sort(eVals.begin(), eVals.end());
    return eVals;
}


//- Mohr-Coulomb yield function (s1 - s3)/2 + (s1 + s3)/2 sin(phi)
//  - c cos(phi) of a stress (tension positive)
inline scalar mcYield(const symmTensor& s, const scalar sinPhi, const scalar cCos)
{
    const vector ev(principal(s));
    const scalar s1 = cmptMax(ev);
    const scalar s3 = cmptMin(ev);
    return 0.5*(s1 - s3) + 0.5*(s1 + s3)*sinPhi - cCos;
}


//- Flow impedance (WP 3.6c, decision 0040). Tangent viscosity
//  d tau/d gamma_dot of the Jop mu(I) law at a cell's stress s, for
//  Mohr-Coulomb: tau = sin(phi_m) p_m with p_m = c cot(phi) - (s1 + s3)/2,
//  mu = tan(phi_m), I = I0 (mu - mu_s)/(mu_2 - mu), P = -tr(s)/3:
//  eta = p_m (1 + mu^2)^-1.5 (mu_2 - mu_s) I0/(I0 + I)^2 d sqrt(rho_s/P).
//  law = (mu_s, mu_2, I0, d, rho_s). Returns 0 for a cell more than
//  flowBand below the static surface in sin(phi) (not flowing; the band
//  keeps cells that hover at the static surface from switching every step,
//  which locked stiff layers) and sets sinMob = sin(phi_m) (0 if p_m <= 0).
static const scalar flowBand = 0.03;

inline scalar flowViscosity
(
    const symmTensor& s,
    const FixedList<scalar, 5>& law,
    const scalar cCot,
    scalar& sinMob
)
{
    sinMob = 0;
    const vector ev(principal(s));
    const scalar s1 = cmptMax(ev);
    const scalar s3 = cmptMin(ev);
    const scalar pm = cCot - 0.5*(s1 + s3);
    const scalar P = -tr(s)/3;
    // P so small that rhoS/P would overflow (a floating-point trap; the
    // stress-free silo fill at its first step, WP 3.6d): no flow viscosity,
    // as at P <= 0 (same guard as FrictionLaw::inertialNumber)
    if (pm <= 0 || !(P > 1e-300*law[4]))
    {
        return 0;
    }
    sinMob = min(0.5*(s1 - s3)/pm, 1 - 1e-12);
    const scalar muS = law[0], mu2 = law[1], I0 = law[2];
    const scalar sinS = muS/Foam::sqrt(1 + sqr(muS));
    if (sinMob < sinS*(1 - flowBand))
    {
        return 0;
    }
    const scalar mu =
        min
        (
            max(sinMob/Foam::sqrt(1 - sqr(sinMob)), muS),
            muS + (1 - 1e-6)*(mu2 - muS)
        );
    const scalar In = I0*(mu - muS)/(mu2 - mu);
    return
        pm*Foam::pow(1 + sqr(mu), -1.5)*(mu2 - muS)*I0/sqr(I0 + In)
       *law[3]*Foam::sqrt(law[4]/P);
}


//- Shear jump coefficient of a face between two cells (flow impedance):
//  each side contributes zs, or 2 eta/delta if flowing (whichever is
//  smaller), in series. delta is the distance between the cell centres.
inline scalar flowJump
(
    const scalar zsL, const scalar etaL,
    const scalar zsR, const scalar etaR,
    const scalar delta
)
{
    const scalar kL = etaL > 0 ? min(zsL, 2*etaL/delta) : zsL;
    const scalar kR = etaR > 0 ? min(zsR, 2*etaR/delta) : zsR;
    return (kL > 0 && kR > 0) ? kL*kR/(kL + kR) : 0;
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


//- Rule 12 (decision 0028; C++ WP 3.4b, decision 0035): normal jump cap at
//  a contact face with unit normal n (L to R), after the shear cap. The
//  face stress is the tangential block of sMean with the n row from the
//  traction tS (tangential part as left by rule 8). The normal jump term
//  zpL zpR (vnR - vnL)/(zpL + zpR) is scaled back (20 bisection steps)
//  until the yield function is no larger than without it (or zero). The
//  face keeps one normal velocity (geometry, remap); each side's own
//  characteristic at the capped traction gives the normal velocity it
//  sees, returned as dL, dR (along n, relative to vS). Returns true if
//  capped, with tS changed.
inline bool capNormal
(
    const vector& n,
    const vector& vL, const symmTensor& sL, const scalar zpL,
    const vector& vR, const symmTensor& sR, const scalar zpR,
    const symmTensor& sMean,
    const scalar sinPhi, const scalar cCos,
    const vector& vS,
    vector& tS,
    vector& dL,
    vector& dR
)
{
    if (zpL <= 0 || zpR <= 0)
    {
        return false;
    }
    const scalar tnL = (sL & n) & n;
    const scalar tnR = (sR & n) & n;
    const scalar vnL = vL & n;
    const scalar vnR = vR & n;
    const scalar zp = zpL + zpR;
    const scalar nb = (zpR*tnL + zpL*tnR)/zp;
    const scalar nj = zpL*zpR*(vnR - vnL)/zp;
    const vector tt(tS - (tS & n)*n);

    const tensor P(I - sqr(n));
    const tensor A(P & tensor(sMean) & P);
    auto faceYield = [&](const scalar tn)
    {
        return mcYield
        (
            symm(A + tn*sqr(n) + (n*tt) + (tt*n)),
            sinPhi,
            cCos
        );
    };

    const scalar lim = max(faceYield(nb), scalar(0));
    if (faceYield(nb + nj) <= lim)
    {
        return false;
    }
    scalar lo = 0, hi = 1;
    for (label k = 0; k < 20; ++k)
    {
        const scalar mid = 0.5*(lo + hi);
        if (faceYield(nb + mid*nj) <= lim)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    const scalar tn = nb + lo*nj;
    tS = tn*n + tt;
    const scalar vnS = vS & n;
    dL = (vnL + (tn - tnL)/zpL - vnS)*n;
    dR = (vnR - (tn - tnR)/zpR - vnS)*n;
    return true;
}

//- Immersed walls (WP 4.3a): a small cut cell and the cell it is mixed
//  with act as one cell. Sums f over each group (target plus its small
//  cells) and gives every member the sum. mixTo: target of each small cell
//  (-1 for others)
template<class Type>
void groupSum(const labelList& mixTo, Field<Type>& f)
{
    const Field<Type> f0(f);
    forAll(mixTo, celli)
    {
        if (mixTo[celli] >= 0)
        {
            f[mixTo[celli]] += f0[celli];
        }
    }
    forAll(mixTo, celli)
    {
        if (mixTo[celli] >= 0)
        {
            f[celli] = f[mixTo[celli]];
        }
    }
}

} // End namespace lagrangeStepDetail
} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eulerianSolid::lagrangeStep(const scalar dt)
{
    if (fluidity_.valid())
    {
        solveFluidity();
    }
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

    // WP F.8g (decision 0133): support by a pore pressure that only lowers
    // the friction: s = min(p, P) on the total stress of the step start
    if (supportP_)
    {
        const scalarField& ps = supportP_->primitiveField();
        const symmTensorField& sg = sigma_.primitiveField();
        forAll(ps, celli)
        {
            const scalar P = -tr(sg[celli])/3;
            const scalar sh = (P > 0 ? min(max(ps[celli], scalar(0)), P) : 0);
            supportShift_[celli] = sh;
            frictionScale_[celli] =
                P > 0 ? max(1 - sh/P, scalar(0)) : (ps[celli] > 0 ? 0 : 1);
        }
    }

    // Outlets open at the start of the first step at or after their
    // opening time (WP 3.3, decision 0033)
    {
        const scalar t0 = mesh.time().value() - dt;
        forAllConstIters(outletOpenTime_, iter)
        {
            if
            (
                wallFriction_.found(iter.key())
             && t0 >= iter.val() - 1e-9*dt
            )
            {
                openOutlets();
                break;
            }
        }
    }

    // Inlets open and close at the start of a step (WP 3.7a, decision 0060)
    updateInlets(mesh.time().value() - dt, dt);

    // Immersed walls (WP 4.3a, decision 0047): covered cells hold ghost
    // values for the reconstruction; the material fills the open volume Vo
    // of each cell and crosses only the open part of each face
    // (Sfo = open fraction times Sf); a wall cell has one more face, the
    // wall, with area vector wallA_ into the body
    const bool imm = immersed_;
    if (imm)
    {
        fillCovered();
    }
    const scalarField Vo(openVolume());
    tmp<surfaceVectorField> tSfo(Sf);
    if (imm)
    {
        tSfo = Sf*faceOpen_();
    }
    const surfaceVectorField& Sfo = tSfo();

    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();

    // Cell mass, constant in the step
    const scalarField m(rho_.primitiveField()*Vo);

    // Inertia of each cell (decision 0094): its mass; but a residue (with free
    // surface, less than 1e-6 of a full cell's mass) gets at least half the
    // rest density times its material volume alpha Vo, as the impedances use
    // the material density rho/alpha bounded to 0.5 to 2 rho0. A trace whose
    // mass and fill disagree (1e-18 kg at alpha 1e-6) cannot then be driven
    // to runaway speeds by its own face forces. Other cells are unchanged.
    scalarField mI(m);
    if (freeSurface_)
    {
        const scalarField& a = alpha_.primitiveField();
        forAll(mI, celli)
        {
            if (a[celli] > alphaEmpty && m[celli] < 1e-6*rho0_*Vo[celli])
            {
                mI[celli] = max(m[celli], 0.5*rho0_*a[celli]*Vo[celli]);
            }
        }
    }

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
        // Open inlets count as empty for the slab side, as traction-free
        // patches (WP 3.7a, decision 0060): counted as full, the stream
        // under an inlet leaned towards it and fell in place
        slabSides
        (
            mesh, alpha_, kinds, lean_, smooth_, gAxis,
            blendVoid_[0], fullBlend_ ? blendVoid_[1] : 0, g_,
            imm ? &faceOpen_() : nullptr,
            lineDirs_.empty() ? nullptr : &lineDirs_,
            &hCell_
        );

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
                    gm*r[celli]*(cmptMag(inLines(celli, gHat)) & hCell_[celli]);
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

    // Flow impedance (WP 3.6c, decision 0040): mu(I) tangent viscosity of
    // flowing full cells (0 elsewhere) and sin(phi_m), with neighbour values
    // on coupled patches
    tmp<volScalarField> teta, tsinMob;
    const bool fi = flowImpedance_ && solveMomentum_;
    if (fi)
    {
        teta = workField<scalar>(mesh, "etaFlow", dimless);
        tsinMob = workField<scalar>(mesh, "sinMob", dimless);
        scalarField& e = teta.ref().primitiveFieldRef();
        scalarField& sm = tsinMob.ref().primitiveFieldRef();
        const scalar cCot =
            capSinPhi_ > 0 ? capCCos_/capSinPhi_ : 0;
        forAll(e, celli)
        {
            const scalar w = 
                fs
              ? fullWeight
                (
                    alpha_[celli], blendVoid_[0],
                    fullBlend_ ? blendVoid_[1] : 0
                )
              : 1;
            if (w > 0)
            {
                // WP F.8g: with frictionOnly the law sees the effective
                // stress sigma + lambda P I
                const symmTensor sEff
                (
                    supportP_
                  ? sigma_[celli]
                  + (1 - frictionScale_[celli])
                   *max(-tr(sigma_[celli])/3, scalar(0))*symmTensor::I
                  : sigma_[celli]
                );
                e[celli] = w*flowViscosity(sEff, flowLaw_, cCot, sm[celli]);
            }
        }
        // Nonlocal fluidity (WP N.2, decision 0111; C.3, decision 0112): a
        // cell flows where the law's secant viscosity P*/g (tau = mu P*,
        // gamma_dot = g mu, P* = P + c cot(phi)) is below 10 zs over its
        // closest neighbour (2 eta dc < 10 zs), including cells below mu_s
        // that creep through the nonlocal term; its face viscosity is the
        // law's tangent d tau/d gamma_dot = b d sqrt(rho_s P*) (as the mu(I)
        // flow impedance), or the secant if smaller. Elsewhere 0, so the
        // shear cap is kept
        if (fluidity_.valid())
        {
            const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();
            scalarField dcMax(e.size(), 0);
            forAll(own, facei)
            {
                dcMax[own[facei]] = max(dcMax[own[facei]], deltaCoeffs[facei]);
                dcMax[nei[facei]] = max(dcMax[nei[facei]], deltaCoeffs[facei]);
            }
            forAll(mesh.boundary(), patchi)
            {
                const fvPatch& p = mesh.boundary()[patchi];
                if (p.coupled())
                {
                    const labelUList& fc = p.faceCells();
                    forAll(fc, i)
                    {
                        dcMax[fc[i]] =
                            max
                            (
                                dcMax[fc[i]],
                                deltaCoeffs.boundaryField()[patchi][i]
                            );
                    }
                }
            }
            const scalarField& g = fluidity_().primitiveField();
            const auto& law = *material_->nonlocal();
            const scalar bdRho = law.b()*law.d()*Foam::sqrt(law.rhoS());
            // Flowing where the secant viscosity P*/g is below 10 z_s over
            // the closest neighbour (the margin is not sensitive: 10 and
            // 100 agree to 0.2 % in the sheared layer, decision 0112)
            const scalar nlK = 10;
            forAll(e, celli)
            {
                const scalar w =
                    fs
                  ? fullWeight
                    (
                        alpha_[celli], blendVoid_[0],
                        fullBlend_ ? blendVoid_[1] : 0
                    )
                  : 1;
                // P* = P + c cot(phi), as the law reads it (decision 0112)
                const scalar P = -tr(sigma_[celli])/3 + cCot;
                e[celli] = 0;
                // (compared without dividing by g, which can be tiny)
                const scalar zsc = Foam::sqrt(G*rho_[celli]);
                if (w > 0 && P > 0 && 2*w*P*dcMax[celli] < nlK*zsc*g[celli])
                {
                    e[celli] = w*min(P/g[celli], bdRho*Foam::sqrt(P));
                }
            }
        }
        teta.ref().correctBoundaryConditions();
        tsinMob.ref().correctBoundaryConditions();
    }
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();
    const volScalarField* smoothPtr = fs ? &tsmooth() : nullptr;
    const volScalarField* eFacPtr = fs ? &teFac() : nullptr;

    // Neighbour alpha and slab sides across coupled patches
    List<scalarField> aNbr(mesh.boundary().size());
    List<vectorField> leanNbr(mesh.boundary().size());
    // Angled meshes: the neighbour cell's mesh-line frame across coupled
    // patches (WP G.6), so a processor face reads the neighbour's side as
    // an internal face does
    List<tensorField> dirsNbr(mesh.boundary().size());
    if (fs)
    {
        tmp<volVectorField> tl = workField<vector>(mesh, "lean", dimless);
        tl.ref().primitiveFieldRef() = lean_;
        tl.ref().correctBoundaryConditions();
        tmp<volTensorField> td;
        if (!lineDirs_.empty())
        {
            td = workField<tensor>(mesh, "lineDirs", dimless);
            td.ref().primitiveFieldRef() = lineDirs_;
            td.ref().correctBoundaryConditions();
        }
        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                aNbr[patchi] =
                    alpha_.boundaryField()[patchi].patchNeighbourField();
                leanNbr[patchi] =
                    tl().boundaryField()[patchi].patchNeighbourField();
                if (td.valid())
                {
                    dirsNbr[patchi] =
                        td().boundaryField()[patchi].patchNeighbourField();
                }
            }
        }
    }
    auto inLinesNbr = [&](const label patchi, const label i, const vector& S)
    {
        return dirsNbr[patchi].empty() ? S : (dirsNbr[patchi][i] & S);
    };


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

    // Wall ghost of a flowing cell (WP 3.6c2, decision 0044): the velocity
    // reconstruction mirrors about the wall's last face velocity (0 where it
    // sticks) instead of the cell value, so a flowing cell's face state at
    // the wall is its extrapolated wall velocity and the face shear jump
    // (eta/distance) does not act on the whole cell velocity a second time
    tmp<volVectorField> tUw;
    if (fi)
    {
        tUw.reset(new volVectorField("Uwall", U_));
        volVectorField& Uw = tUw.ref();
        forAll(mesh.boundary(), patchi)
        {
            if (kinds[patchi] != boundaryKind::wall)
            {
                continue;
            }
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            fvPatchVectorField& pU = Uw.boundaryFieldRef()[patchi];
            const vectorField& pv = vStar_.boundaryField()[patchi];
            forAll(fc, i)
            {
                if (teta()[fc[i]] > 0)
                {
                    pU[i] = pv[i];
                }
            }
        }
    }
    faceStates
    (
        tUw.valid() ? tUw() : U_, secondOrder_, UO.ref(), UN.ref(), nullptr, Zero,
        smoothPtr, nullptr, angled_ ? &lineRot_ : nullptr
    );
    // Axisymmetric wedge (WP 3.8a, decision 0061): wall stress extrapolated
    tmp<volSymmTensorField> tsWalls;
    if (wedge_ && secondOrder_)
    {
        tsWalls = extrapolatedWalls(sref, g);
    }
    faceStates
    (
        tsWalls.valid() ? tsWalls() : sref, secondOrder_, sO.ref(), sN.ref(),
        &rho_, g, smoothPtr, eFacPtr, angled_ ? &lineRot_ : nullptr
    );

    // Cell state the wall face of an immersed wall cell sees: the cell
    // value plus the Hancock predictor increment (second order)
    vectorField Uwall;
    symmTensorField Swall;
    if (imm)
    {
        Uwall = U_.primitiveField();
        Swall = sref.primitiveField();
    }

    // Predicted cell state (Hancock half step; the cell state when first
    // order): what a wedge side carries (WP 3.8b, decision 0061), so the
    // wedge faces move with the other faces of the cell
    vectorField Upred;
    symmTensorField Spred;
    if (wedge_)
    {
        Upred = U_.primitiveField();
        Spred = sref.primitiveField();
    }

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
            const vector& S = Sfo[facei];
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
                    if (!fs || touches(alpha_[P], lean_[P], inLines(P, pS[i])))
                    {
                        wallFace
                        (
                            pS[i]/pmS[i], pu[i], ps[i], zp[P], zs[P],
                            mu*(supportP_ ? frictionScale_[P] : scalar(1)),
                            v, t,
                            fi && teta()[P] > 0
                          ? min
                            (
                                zs[P],
                                teta()[P]
                               *deltaCoeffs.boundaryField()[patchi][i]
                            )
                          : -1
                        );
                    }
                    const scalar fo =
                        imm ? faceOpen_().boundaryField()[patchi][i] : 1;
                    divS[P] += fo*(t*pmS[i] - (sref[P] & pS[i]));
                    gradU[P] += fo*pS[i]*(v - U_[P]);
                }
                continue;
            }
            const vectorField& pSo = Sfo.boundaryField()[patchi];
            if (isA<wedgeFvPatch>(p))
            {
                // Wedge side (decision 0061): the rotated cell state, as in
                // the corrector, so a cell at rest stays at rest
                const tensor& Tw = refCast<const wedgeFvPatch>(p).faceT();
                forAll(fc, i)
                {
                    const label P = fc[i];
                    divS[P] += (transform(Tw, sref[P]) - sref[P]) & pSo[i];
                    gradU[P] += pSo[i]*(transform(Tw, U_[P]) - U_[P]);
                }
                continue;
            }
            forAll(fc, i)
            {
                divS[fc[i]] += (ps[i] - sref[fc[i]]) & pSo[i];
                gradU[fc[i]] += pSo[i]*(pu[i] - U_[fc[i]]);
            }
        }
        if (imm)
        {
            // Immersed wall faces: wall Riemann solution from the cell
            // state (WP 4.3a)
            forAll(wallA_, celli)
            {
                if (wallBody_[celli] < 0 || m[celli] <= 0)
                {
                    continue;
                }
                const vector& A = wallA_[celli];
                const scalar mA = mag(A);
                vector v(U_[celli]), t(Zero);
                if (!fs || touches(alpha_[celli], lean_[celli], inLines(celli, A)))
                {
                    wallFace
                    (
                        A/mA, U_[celli] - wallU_[celli], sref[celli],
                        zp[celli], zs[celli], bodyMu_[wallBody_[celli]], v, t,
                        -1, bodyAdhesion(wallBody_[celli])
                    );
                    v += wallU_[celli];
                }
                divS[celli] += t*mA - (sref[celli] & A);
                gradU[celli] += A*(v - U_[celli]);
            }
        }

        // Small cut cells move and strain with the cell they are mixed
        // with (one cell for the predictor)
        scalarField mP(mI), VoP(Vo);
        if (imm)
        {
            groupSum(mixTo_, divS);
            groupSum(mixTo_, gradU);
            groupSum(mixTo_, mP);
            groupSum(mixTo_, VoP);
        }
        forAll(V, celli)
        {
            if (imm && cellOpen_[celli] <= 0)
            {
                continue;       // covered by a body
            }
            if (fs && smooth_[celli] == 0)
            {
                // First-order cell (decision 0021): no stress predictor;
                // empty cells do not move
                if (solveMomentum_ && m[celli] > 0)
                {
                    dU[celli] = 0.5*dt*(divS[celli]/mP[celli] + g);
                }
                continue;
            }
            if (solveMomentum_)
            {
                dU[celli] = 0.5*dt*(divS[celli]/mP[celli] + g);
            }
            const symmTensor D(symm(gradU[celli])/VoP[celli]);
            symmTensor Dd(dev(D));
            if (fi && teta()[celli] > 0)
            {
                // Flowing cell: no elastic stress rate along the plastic
                // flow direction m = (e1 e1 - e3 e3)/sqrt(2) while loading
                const vector ev(principal(sigma_[celli]));
                const tensor E(eigenVectors(sigma_[celli], ev));
                direction i1 = 0, i3 = 0;
                for (direction k = 1; k < 3; ++k)
                {
                    if (ev[k] > ev[i1]) i1 = k;
                    if (ev[k] < ev[i3]) i3 = k;
                }
                if (i1 != i3)
                {
                    const vector e1(E.row(i1)), e3(E.row(i3));
                    const symmTensor mf
                    (
                        (sqr(e1) - sqr(e3))/Foam::sqrt(2.0)
                    );
                    const scalar dm = mf && Dd;
                    if (dm > 0)
                    {
                        Dd -= dm*mf;
                    }
                }
            }
            dS[celli] = 0.5*dt*(K*tr(D)*symmTensor::I + 2*G*Dd);
            if (fs && smooth_[celli] < 1)
            {
                // Blended cell (WP 3.6h): part of the stress predictor
                dS[celli] *= smooth_[celli];
            }
        }
        dU.correctBoundaryConditions();
        dS.correctBoundaryConditions();
        if (wedge_)
        {
            Upred += dU.primitiveField();
            Spred += dS.primitiveField();
        }
        if (imm)
        {
            Uwall += dU.primitiveField();
            Swall += dS.primitiveField();
        }

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

    // Normal cap (rule 12): normal velocity seen by the owner and the
    // neighbour cell of a capped face, relative to vStar (zero elsewhere)
    tmp<surfaceVectorField> tnO, tnN;
    nNormalCapped_ = 0;
    if (normalCap_)
    {
        tnO = surfaceVectorField::New
        (
            "nslipO", IOobject::NO_REGISTER, mesh,
            dimensionedVector(dimVelocity, Zero)
        );
        tnN = surfaceVectorField::New
        (
            "nslipN", IOobject::NO_REGISTER, mesh,
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
        scalar kj = -1;
        bool bothFlow = false;
        scalar capSin = capSinPhi_;
        if (fi)
        {
            const scalar eP = teta()[P], eN = teta()[N];
            if (eP > 0 || eN > 0)
            {
                kj = flowJump(zs[P], eP, zs[N], eN, 1/deltaCoeffs[facei]);
                bothFlow = eP > 0 && eN > 0;
                if (eP > 0) capSin = max(capSin, tsinMob()[P]);
                if (eN > 0) capSin = max(capSin, tsinMob()[N]);
            }
        }
        if (fs)
        {
            const scalar aP = alpha_[P], aN = alpha_[N];
            contact = voidFace
            (
                n,
                UO()[facei], sO()[facei], zp[P], zs[P], aP > alphaEmpty,
                touches(aP, lean_[P], inLines(P, Sf[facei])),
                UN()[facei], sN()[facei], zp[N], zs[N], aN > alphaEmpty,
                touches(aN, lean_[N], inLines(N, -Sf[facei])),
                tensileStrength_, vS[facei], t, kj
            );
        }
        else
        {
            riemann
            (
                n,
                UO()[facei], sO()[facei], zp[P], zs[P],
                UN()[facei], sN()[facei], zp[N], zs[N],
                vS[facei], t, kj
            );
        }
        if (shearCap_ && contact && solveMomentum_ && !bothFlow)
        {
            nCapped_ += capShear
            (
                n,
                UO()[facei], sO()[facei], zs[P],
                UN()[facei], sN()[facei], zs[N],
                0.5*(sref[P] + sref[N]),
                capSin, capCCos_, vS[facei], t,
                tdO.ref()[facei], tdN.ref()[facei]
            );
        }
        if (normalCap_ && contact && solveMomentum_)
        {
            nNormalCapped_ += capNormal
            (
                n,
                UO()[facei], sO()[facei], zp[P],
                UN()[facei], sN()[facei], zp[N],
                0.5*(sref[P] + sref[N]),
                capSinPhi_, capCCos_, vS[facei], t,
                tnO.ref()[facei], tnN.ref()[facei]
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

    // Surface loads at the start of the step (WP F.3, decision 0128)
    scalarList qLoad(mesh.boundary().size(), Zero);
    if (loadPressure_.size())
    {
        forAll(mesh.boundary(), patchi)
        {
            qLoad[patchi] = loadPressure
            (
                mesh.boundary()[patchi].name(), mesh.time().value() - dt
            );
        }
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
            if (shearCap_ || normalCap_)
            {
                srn = sref.boundaryField()[patchi].patchNeighbourField();
            }
            scalarField etan, smn;
            if (fi)
            {
                etan = teta().boundaryField()[patchi].patchNeighbourField();
                smn = tsinMob().boundaryField()[patchi].patchNeighbourField();
            }
            forAll(fc, i)
            {
                const vector n(pS[i]/pmS[i]);
                vector t;
                bool contact = true;
                scalar kj = -1;
                bool bothFlow = false;
                scalar capSin = capSinPhi_;
                if (fi)
                {
                    const scalar eP = teta()[fc[i]], eN = etan[i];
                    if (eP > 0 || eN > 0)
                    {
                        kj = flowJump
                        (
                            zs[fc[i]], eP, zsn[i], eN,
                            1/deltaCoeffs.boundaryField()[patchi][i]
                        );
                        bothFlow = eP > 0 && eN > 0;
                        if (eP > 0) capSin = max(capSin, tsinMob()[fc[i]]);
                        if (eN > 0) capSin = max(capSin, smn[i]);
                    }
                }
                if (fs)
                {
                    const scalar aP = alpha_[fc[i]];
                    const scalar aN = aNbr[patchi][i];
                    contact = voidFace
                    (
                        n,
                        puo[i], pso[i], zp[fc[i]], zs[fc[i]], aP > alphaEmpty,
                        touches(aP, lean_[fc[i]], inLines(fc[i], pS[i])),
                        pun[i], psn[i], zpn[i], zsn[i], aN > alphaEmpty,
                        touches(aN, leanNbr[patchi][i], inLinesNbr(patchi, i, -pS[i])),
                        tensileStrength_, pvS[i], t, kj
                    );
                }
                else
                {
                    riemann
                    (
                        n,
                        puo[i], pso[i], zp[fc[i]], zs[fc[i]],
                        pun[i], psn[i], zpn[i], zsn[i],
                        pvS[i], t, kj
                    );
                }
                if (shearCap_ && contact && solveMomentum_ && !bothFlow)
                {
                    nCapped_ += capShear
                    (
                        n,
                        puo[i], pso[i], zs[fc[i]],
                        pun[i], psn[i], zsn[i],
                        0.5*(sref[fc[i]] + srn[i]),
                        capSin, capCCos_, pvS[i], t,
                        tdO.ref().boundaryFieldRef()[patchi][i],
                        tdN.ref().boundaryFieldRef()[patchi][i]
                    );
                }
                if (normalCap_ && contact && solveMomentum_)
                {
                    nNormalCapped_ += capNormal
                    (
                        n,
                        puo[i], pso[i], zp[fc[i]],
                        pun[i], psn[i], zpn[i],
                        0.5*(sref[fc[i]] + srn[i]),
                        capSinPhi_, capCCos_, pvS[i], t,
                        tnO.ref().boundaryFieldRef()[patchi][i],
                        tnN.ref().boundaryFieldRef()[patchi][i]
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

        const label ink = inletOfPatch_[patchi];
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
            if (kind == boundaryKind::wedge)
            {
                // Axisymmetric wedge side (WP 3.8a, decision 0061): the
                // state beyond is the cell's own (predicted), turned onto the face's
                // plane (faceT, as OpenFOAM's wedge patch field). No jump,
                // so the face carries the rotated cell velocity and stress:
                // with the other side, the hoop strain rate and the hoop
                // stress. A swirl velocity crosses the face, but what
                // leaves one side enters the other, so no volume flux
                // (phiL = 0 below)
                const tensor& T = refCast<const wedgeFvPatch>(p).faceT();
                v = transform(T, Upred[P]);
                t = transform(T, Spred[P]) & n;
            }
            else if (ink >= 0 && inletOpen_[ink])
            {
                // Open inlet (WP 3.7a, decision 0060): the face moves into
                // the domain at the inlet speed; the remap fills the volume
                // it sweeps with fresh material. Traction from the interior
                // characteristic where the cell's material is at the face
                // and pressed against the fed material; else (no material
                // at the face, or it pulls away: no tension across the
                // feed) the fed material's own traction -p n
                v = -inletNow_[ink]*n;
                t = -inletP_[ink]*n;
                if (touches(alpha_[P], lean_[P], inLines(P, pS[i])))
                {
                    const scalar vn = v & n;
                    const scalar tn = tnL + zp[P]*(vn - vnL);
                    if (tn < -inletP_[ink])
                    {
                        t = tn*n + ttL - zs[P]*vtL;
                    }
                }
            }
            else if (fs && !touches(alpha_[P], lean_[P], inLines(P, pS[i])))
            {
                // Material not at the face (decision 0021 rule 3). An open
                // outlet takes the free one-sided velocity anyway, so the
                // void under the slab passes out first (rule 10b, decisions
                // 0027 and 0033); other faces stay still.
                t = Zero;
                v = Zero;
                if
                (
                    alpha_[P] > alphaEmpty
                 && kind == boundaryKind::traction
                 && outletOpenTime_.found(p.name())
                )
                {
                    v = (vnL - tnL/zp[P])*n + vtL - ttL/zs[P];
                }
            }
            else if (kind == boundaryKind::velocity)
            {
                // Prescribed velocity; traction from the interior
                // characteristic s* - z v* = sL - z vL
                v = pU[i];
                const scalar vn = v & n;
                const vector vt(v - vn*n);
                // Flowing cell with nonlocal fluidity: its viscosity over
                // the half cell, if smaller than zs (WP N.2; unchanged
                // for mu(I) flow impedance)
                const scalar ks =
                    fi && fluidity_.valid() && teta()[P] > 0
                  ? min(zs[P], teta()[P]*deltaCoeffs.boundaryField()[patchi][i])
                  : zs[P];
                t = (tnL + zp[P]*(vn - vnL))*n + ttL + ks*(vt - vtL);
            }
            else if (kind == boundaryKind::wall)
            {
                wallFace
                (
                    n, vL, pso[i], zp[P], zs[P],
                    wallFriction_[p.name()]
                   *(supportP_ ? frictionScale_[P] : scalar(1)),
                    v, t,
                    fi && teta()[P] > 0
                  ? min(zs[P], teta()[P]*deltaCoeffs.boundaryField()[patchi][i])
                  : -1
                );
            }
            else if (kind == boundaryKind::symmetry)
            {
                // No normal velocity, no shear traction
                v = vtL - ttL/zs[P];
                t = (tnL - zp[P]*vnL)*n;
            }
            else if (qLoad[patchi] != 0)
            {
                // Surface load (WP F.3, decision 0128): normal traction
                // -q n, the face moves as a free surface under it
                const scalar q = qLoad[patchi];
                v = (vnL - (tnL + q)/zp[P])*n + vtL - ttL/zs[P];
                t = -q*n;
            }
            else
            {
                // Free surface: zero traction
                v = (vnL - tnL/zp[P])*n + vtL - ttL/zs[P];
                t = Zero;
            }
            pvS[i] = v;
            pF[i] = t*pmS[i];
            pphi[i] = kind == boundaryKind::wedge ? 0 : (v & pS[i]);
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
    if (normalCap_)
    {
        // copies vS, F and phiL again (unchanged) and the normal slips
        syncCoupled(mesh, vS, F, phiL_, &tnO.ref(), &tnN.ref());
    }
    if (imm)
    {
        // Only the open part of a face carries force and volume (WP 4.3a)
        F *= faceOpen_();
        phiL_ *= faceOpen_();
    }
    vStar_ = vS;
    Fstar_ = F;


    // 4. Momentum, velocity gradient and Lagrangian volume ----------------------

    vectorField force(mesh.nCells(), Zero);
    tensorField T(mesh.nCells(), Zero);      // sum Sf * vStar
    scalarField dV(mesh.nCells(), Zero);     // sum phiL out of the cell
    tensorField Tn;                          // sum Sf * normal slip (rule 12)
    if (normalCap_)
    {
        Tn.setSize(mesh.nCells(), Zero);
        const surfaceVectorField& nO = tnO();
        const surfaceVectorField& nN = tnN();
        forAll(own, facei)
        {
            Tn[own[facei]] += Sfo[facei]*nO[facei];
            Tn[nei[facei]] -= Sfo[facei]*nN[facei];
        }
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const vectorField& pS = Sfo.boundaryField()[patchi];
            const vectorField& pnO = nO.boundaryField()[patchi];
            forAll(fc, i)
            {
                Tn[fc[i]] += pS[i]*pnO[i];
            }
        }
    }
    vectorField dVa;                         // the same per face axis (fs)
    if (fs)
    {
        dVa.setSize(mesh.nCells(), Zero);
        forAll(own, facei)
        {
            dVa[own[facei]][faceAxis(inLines(own[facei], Sf[facei]))] +=
                phiL_[facei];
            dVa[nei[facei]][faceAxis(inLines(nei[facei], Sf[facei]))] -=
                phiL_[facei];
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            const labelUList& fc = p.faceCells();
            const vectorField& pS = Sf.boundaryField()[patchi];
            const scalarField& pphi = phiL_.boundaryField()[patchi];
            forAll(fc, i)
            {
                dVa[fc[i]][faceAxis(inLines(fc[i], pS[i]))] += pphi[i];
            }
        }
    }
    if (wedge_)
    {
        // Hoop volume rate from the wedge sides' face velocities (decision
        // 0061): the radial faces' flux (r u)' holds both u' and u/r; the
        // u/r part belongs to the hoop axis
        hoopRate_ = 0;
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!isA<wedgeFvPatch>(p)) continue;
            const labelUList& fc = p.faceCells();
            const vectorField& pS = Sfo.boundaryField()[patchi];
            const vectorField& pvS = vS.boundaryField()[patchi];
            forAll(fc, i)
            {
                hoopRate_[fc[i]] += pvS[i] & pS[i];
            }
        }
        if (fs)
        {
            forAll(dVa, celli)
            {
                dVa[celli][wedgeR_] -= hoopRate_[celli];
                dVa[celli][wedgeH_] += hoopRate_[celli];
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
            T[P] += Sfo[facei]*(vS[facei] + tdO()[facei]);
            T[N] -= Sfo[facei]*(vS[facei] + tdN()[facei]);
        }
        else
        {
            const tensor t(Sfo[facei]*vS[facei]);
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
        const vectorField& pS = Sfo.boundaryField()[patchi];
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

    if (imm)
    {
        // Immersed wall faces (WP 4.3a, decision 0047): one-sided Riemann
        // solution against the body (no tension, Coulomb shear cap with the
        // body's friction), from the predicted cell state, in the body's
        // frame. The force on the body is minus the force on the material.
        bodyForce_ = Zero;
        wallForce_.setSize(wallA_.size());
        wallForce_ = Zero;
        if (fs)
        {
            wallPhi_ = 0;
        }
        forAll(wallA_, celli)
        {
            if (wallBody_[celli] < 0 || m[celli] <= 0)
            {
                continue;
            }
            const label b = wallBody_[celli];
            const vector& A = wallA_[celli];
            const scalar mA = mag(A);
            if (fs && !touches(alpha_[celli], lean_[celli], inLines(celli, A)))
            {
                // Material not at the wall (WP 4.3c): the wall face stays
                // with the body, no force
                T[celli] += A*wallU_[celli];
                wallPhi_[celli] = wallU_[celli] & A;
                dV[celli] += wallPhi_[celli];
                dVa[celli][faceAxis(inLines(celli, A))] += wallPhi_[celli];
                continue;
            }
            vector v, t;
            wallFace
            (
                A/mA, Uwall[celli] - wallU_[celli], Swall[celli],
                zp[celli], zs[celli], bodyMu_[b], v, t, -1, bodyAdhesion(b)
            );
            v += wallU_[celli];
            force[celli] += t*mA;
            T[celli] += A*v;
            dV[celli] += v & A;
            bodyForce_[b] += t*mA;
            wallForce_[celli] = t*mA;
            if (fs)
            {
                // With void (WP 4.3c, decision 0051) material that pulls
                // away from the wall opens void in the remap: the wall
                // face's volume flux counts in the Lagrangian volume, along
                // the axis of the face's largest component (as a staircase
                // face), and no mass enters
                wallPhi_[celli] = v & A;
                dVa[celli][faceAxis(inLines(celli, A))] += wallPhi_[celli];
                continue;
            }
            // Material that enters through the wall face (where it pulls
            // away from the body: the face is then a free surface, as an
            // open top), counted in the leak report
            const scalar VL = Vo[celli] + dt*dV[celli];
            if (VL > 0)
            {
                bodyMassIn_[b] -= dt*(v & A)*m[celli]/VL;
            }
        }
    }
    // Pore gas (WP F.3, decision 0128): sigma is the effective stress; the
    // gas pressure acts on each cell as -sum(p_f Sf) over its open faces
    // (and its immersed wall face), so a uniform pressure exerts no force
    // and the total stress sigma - p I is in equilibrium with the loads.
    // In a partly filled cell it acts on the material part, alpha.
    if (poreP_)
    {
        // Face values from the gas (WP F.4: the void pressure on faces
        // between material and void), else linear (F.3)
        tmp<surfaceScalarField> tpf
        (
            porePf_ ? tmp<surfaceScalarField>(*porePf_)
                    : fvc::interpolate(*poreP_)
        );
        const surfaceScalarField& pf = tpf();
        vectorField fGas(mesh.nCells(), Zero);
        forAll(own, facei)
        {
            const vector f(pf[facei]*Sfo[facei]);
            fGas[own[facei]] -= f;
            fGas[nei[facei]] += f;
        }
        forAll(mesh.boundary(), patchi)
        {
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            const vectorField& pS = Sfo.boundaryField()[patchi];
            const scalarField& pp = pf.boundaryField()[patchi];
            forAll(fc, i)
            {
                fGas[fc[i]] -= pp[i]*pS[i];
            }
        }
        if (imm)
        {
            const scalarField& pc = poreP_->primitiveField();
            forAll(wallA_, celli)
            {
                if (wallBody_[celli] >= 0)
                {
                    fGas[celli] -= pc[celli]*wallA_[celli];
                }
            }
        }
        if (gasForce_)
        {
            // Reaction of the pore gas inertia (WP F.7)
            fGas += gasForce_->primitiveField();
        }
        forAll(fGas, celli)
        {
            if (m[celli] > 0)
            {
                force[celli] += (fs ? alpha_[celli] : 1)*fGas[celli];
            }
        }
    }

    // Small cut cells: one cell with the cell they are mixed with for the
    // momentum and the strain (decision 0047); each keeps its own volume
    // change for its density
    scalarField mG(mI), VoG(Vo), dVG(dV);
    tensorField TG(T);
    vectorField dVaG(dVa);
    if (imm)
    {
        if (fs)
        {
            groupSum(mixTo_, dVaG);
        }
        groupSum(mixTo_, force);
        groupSum(mixTo_, mG);
        groupSum(mixTo_, VoG);
        groupSum(mixTo_, dVG);
        groupSum(mixTo_, TG);
    }
    if (solveMomentum_)
    {
        vectorField& U = U_.primitiveFieldRef();
        const bool damp = mesh_.time().value() >= dampingStart_ - SMALL;
        forAll(U, celli)
        {
            if (imm && cellOpen_[celli] <= 0)
            {
                U[celli] = wallU_[celli];
                continue;
            }
            if (fs && m[celli] <= 0)
            {
                U[celli] = Zero;
                continue;
            }
            U[celli] += dt*(force[celli]/mG[celli] + g);
            if (damping_ > 0 && damp)
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

    // Axes of a partly filled cell open to plastic flow into its void
    // (voidPlasticFlow, decision 0095)
    static const direction diag3[3] =
        {symmTensor::XX, symmTensor::YY, symmTensor::ZZ};
    label nVoidAxes = 0;
    direction voidAxis[3];
    scalar voidDm[3], voidDm0[3], voidDmMax[3];
    scalar voidA = 1, voidMatV = 0;
    // Angled meshes (WP G.3, decision 0106): in a cell whose mesh lines
    // are square to each other but turned, the void closure works on the
    // normal components along the lines (identity elsewhere: x, y, z)
    bool voidRot = false;
    tensor voidR(tensor::I);

    forAll(sigma, celli)
    {
        nVoidAxes = 0;
        if (imm && cellOpen_[celli] <= 0)
        {
            continue;       // covered by a body: ghost values stay
        }
        scalar VL = Vo[celli] + dt*dV[celli];
        if (imm && fs && VoG[celli] != Vo[celli])
        {
            // Mixing group with void (WP 4.3c, decision 0051): each member
            // takes the group's volume ratio, so a sliver cell's own
            // volume change cannot clip its material volume; the remap
            // uses the same Lagrangian volume
            VL = Vo[celli]*(VoG[celli] + dt*dVG[celli])/VoG[celli];
        }
        const scalar Vmid =
            imm
          ? VoG[celli] + 0.5*dt*dVG[celli]
          : 0.5*(Vo[celli] + VL);

        // T/V = grad(v) in OpenFOAM order (d_i v_j); L = grad(v)^T
        const tensor L(dt*TG[celli].T()/Vmid);
        symmTensor dEps(symm(L));
        const tensor dW(skew(L));

        // Rule 12: a full cell next to capped faces strains with its own
        // characteristic velocities there, deviatoric part only; the
        // volume stays tied to the geometry (decision 0028)
        symmTensor dDev(Zero);
        if (normalCap_ && (!fs || alpha_[celli] >= alphaFull))
        {
            const symmTensor e(symm(dt*Tn[celli].T()/Vmid));
            const Vector<label>& gD = mesh.geometricD();
            symmTensor Pg(Zero);
            if (gD[0] >= 0) Pg.xx() = 1;
            if (gD[1] >= 0) Pg.yy() = 1;
            if (gD[2] >= 0) Pg.zz() = 1;
            const scalar nD = tr(Pg);
            dDev = e - (tr(e)/nD)*Pg;
        }

        if (fs)
        {
            // Void closure (decision 0021 rule 5), per axis of the
            // Cartesian cell: compression closes the void first; expansion
            // first unloads the material to its tensile strength (zero
            // when cohesionless; no cap for a full cell of a material
            // without a strength limit, decision 0054), the rest
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
            voidRot = !lineRot_.empty() && lineRot_[celli] != tensor::I;
            if (voidRot)
            {
                voidR = lineRot_[celli];
            }
            symmTensor dE(voidRot ? transform(voidR, dEps) : dEps);
            const symmTensor sF
            (
                voidRot ? transform(voidR, sigma[celli]) : sigma[celli]
            );
            nVoidAxes = 0;
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                if (geomD[d] < 0 && label(d) != wedgeH_)
                {
                    continue;
                }
                // relative to the open volume; small cut cells: their
                // mixing group's (WP 4.3c)
                const scalar r =
                    imm
                  ? dt*dVaG[celli][d]/VoG[celli]
                  : dt*dVa[celli][d]/V[celli];
                scalar dm;
                if (r < 0)
                {
                    dm = full ? r : min(scalar(0), r + (1 - a));
                }
                else
                {
                    // the material carries tension up to its strength
                    // sT; a full cell of a material without a strength
                    // limit (elastic, J2) never opens void (decision 0054)
                    const scalar sT =
                        tensileStrength_ < GREAT ? tensileStrength_ : 0;
                    const scalar cap =
                        (full && tensileStrength_ >= GREAT)
                      ? GREAT
                      : a*max(scalar(0), sT - sF.component(diag[d]))/M;
                    dm = min(r, cap);
                }
                if (dm != r)
                {
                    same = false;
                }
                if (voidPlasticFlow_ && !full)
                {
                    // the material may take more normal strain on this
                    // axis, up to closing the void left after r
                    voidAxis[nVoidAxes] = d;
                    voidDm[nVoidAxes] = dm;
                    voidDm0[nVoidAxes] = dm;
                    voidDmMax[nVoidAxes] = max(dm, r + (1 - a));
                    ++nVoidAxes;
                }
                matV += dm;
                dE.component(diag[d]) = dm/(a + 0.5*dm);
            }
            voidA = a;
            voidMatV = matV - (full ? 1 : a);
            if (!same)
            {
                dEps = voidRot ? transform(voidR.T(), dE) : dE;
                alphaL_[celli] =
                    min(max(matV*Vo[celli]/VL, scalar(0)), scalar(1));
            }
            else
            {
                alphaL_[celli] = 1;
            }
        }

        dEps += dDev;

        if (dW != tensor::zero)
        {
            const tensor Q(hughesWinget(dW));
            sigma[celli] = transform(Q, sigma[celli]);
            epsP[celli] = transform(Q, epsP[celli]);
        }

        if (dEps != symmTensor::zero || fluidity_.valid())
        {
            constitutiveModel::history hOld, hNew;
            constitutiveModel::toVoigt(epsP[celli], hOld.epsP, 2);
            hOld.epsPEq = epsPEq[celli];
            hOld.pc = pc_.valid() ? pc_()[celli] : 0;

            // WP F.8g: the constitutive update sees sigma + s I with
            // s = lambda P_trial (lambda = min(p/P, 1) of the step start,
            // P_trial the elastic trial mean pressure), so a fully
            // supported cell keeps an elastic, isotropic total stress
            const scalar shift = supportP_
              ? (1 - frictionScale_[celli])
               *max(-tr(sigma[celli])/3 - K*tr(dEps), scalar(0))
              : 0;
            if (shift != 0)
            {
                sigma[celli] += shift*symmTensor::I;
            }

            symmTensor sNew;
            if (fluidity_.valid())
            {
                material_->updateFluidity
                (
                    dEps, dt, fluidity_()[celli], sigma[celli], hOld, sNew,
                    hNew
                );
            }
            else
            {
                material_->update(dEps, dt, sigma[celli], hOld, sNew, hNew);
            }

            // A.4 (decision 0095): a partly filled cell lets the plastic
            // normal strain of its material go into the void (the elastic
            // rules above only see unloading). A few fixed-point passes:
            // each adds the plastic extension of the last trial on the
            // cell's axes, within the void the cell has.
            if (nVoidAxes > 0)
            {
                const symmTensor eOld(constitutiveModel::fromVoigt(hOld.epsP));
                // strain along the turned lines (angled meshes)
                symmTensor dEF(voidRot ? transform(voidR, dEps) : dEps);
                for (label pass = 0; pass < 3; ++pass)
                {
                    const symmTensor dPg
                    (
                        constitutiveModel::fromVoigt(hNew.epsP) - eOld
                    );
                    const symmTensor dP
                    (
                        voidRot ? transform(voidR, dPg) : dPg
                    );
                    bool changed = false;
                    for (label k = 0; k < nVoidAxes; ++k)
                    {
                        const direction d = voidAxis[k];
                        const scalar ext = dP.component(diag3[d]);
                        if (ext <= 0)
                        {
                            continue;
                        }
                        const scalar dm = min
                        (
                            voidDmMax[k],
                            voidDm[k] + voidA*ext
                        );
                        if (dm > voidDm[k])
                        {
                            voidDm[k] = dm;
                            dEF.component(diag3[d]) = dm/(voidA + 0.5*dm);
                            changed = true;
                        }
                    }
                    if (!changed)
                    {
                        break;
                    }
                    // the axes share the cell's void: the material volume
                    // stays within the Lagrangian cell volume
                    scalar extra = 0, room = VL/Vo[celli] - voidA;
                    for (label k = 0; k < nVoidAxes; ++k)
                    {
                        extra += voidDm[k] - voidDm0[k];
                        room -= voidDm0[k];
                    }
                    if (extra > max(room, scalar(0)))
                    {
                        const scalar f = max(room, scalar(0))/extra;
                        for (label k = 0; k < nVoidAxes; ++k)
                        {
                            voidDm[k] =
                                voidDm0[k] + f*(voidDm[k] - voidDm0[k]);
                            dEF.component(diag3[voidAxis[k]]) =
                                voidDm[k]/(voidA + 0.5*voidDm[k]);
                        }
                    }
                    dEps = voidRot ? transform(voidR.T(), dEF) : dEF;
                    material_->update
                    (
                        dEps, dt, sigma[celli], hOld, sNew, hNew
                    );
                }
                scalar matV = voidA;
                for (label k = 0; k < nVoidAxes; ++k)
                {
                    matV += voidDm[k];
                }
                if (matV != voidA + voidMatV)
                {
                    alphaL_[celli] =
                        min(max(matV*Vo[celli]/VL, scalar(0)), scalar(1));
                }
            }

            sigma[celli] = sNew - shift*symmTensor::I;
            symmTensor e(constitutiveModel::fromVoigt(hNew.epsP));
            e.xy() *= 0.5;
            e.yz() *= 0.5;
            e.xz() *= 0.5;
            epsP[celli] = e;
            epsPEq[celli] = hNew.epsPEq;
            if (pc_.valid())
            {
                pc_()[celli] = hNew.pc;
            }
        }

        // Rest relaxation (decision 0114): horizontal block of a resting
        // full cell towards K0 times the vertical stress
        if
        (
            restTime_ > 0
         && alpha_[celli] >= alphaFull
         && mag(U_[celli]) < restSpeed_
        )
        {
            const vector n(g_/mag(g_));
            symmTensor& s = sigma[celli];
            const scalar sv = n & (s & n);
            if (sv < 0)
            {
                const symmTensor Ph(symmTensor::I - sqr(n));
                const symmTensor sh(symm(Ph & s & Ph));
                s += min(dt/restTime_, scalar(1))*(restK0_*sv*Ph - sh);
            }
        }

        // 6. Lagrangian density
        rho[celli] = m[celli]/VL;

        // WP F.8f (decision 0133): packing-dependent frictional pressure
        // (Princeton, Breard et al. 2019): the effective mean stress is at
        // most Pc(solid fraction of the filled part); above it the stress
        // is scaled towards the apex (the grains cannot carry it at this
        // packing and the cell compacts)
        const scalar aP = freeSurface_ ? alphaL_[celli] : scalar(1);
        if (princetonFr_ > 0 && aP > 0 && m[celli] > 0)
        {
            const scalar es = m[celli]/(VL*max(aP, SMALL)*princetonRhoS_);
            scalar Pc = 0;
            if (es > princetonMin_)
            {
                const scalar gap = max(princetonMax_ - es, 1e-6);
                Pc = princetonFr_*Foam::pow(es - princetonMin_, princetonR_)
                    /Foam::pow(gap, princetonS_);
            }
            symmTensor& s = sigma[celli];
            const scalar P = -tr(s)/3;
            if (P > Pc)
            {
                s *= (P > 0 ? Pc/P : 0);
            }
        }
    }

    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    sigma_.correctBoundaryConditions();
    epsP_.correctBoundaryConditions();
    epsPEq_.correctBoundaryConditions();
    if (pc_.valid())
    {
        pc_().correctBoundaryConditions();
    }
}


// ************************************************************************* //
