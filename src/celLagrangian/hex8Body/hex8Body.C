/*---------------------------------------------------------------------------*\
    hex8Body: explicit hexahedral Lagrangian body (WP 5.1, decision 0070)
\*---------------------------------------------------------------------------*/

#include "hex8Body.H"
#include "tensor.H"
#include "transform.H"

#include <utility>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hex8Body, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- Faces of a hex, outward for a positive element (OpenFOAM hex model);
//  node k of a face sits at (xi, eta) = (-1,-1), (1,-1), (1,1), (-1,1)
static const label hexFaces[6][4] =
{
    {0, 4, 7, 3},
    {1, 2, 6, 5},
    {0, 1, 5, 4},
    {3, 7, 6, 2},
    {0, 3, 2, 1},
    {4, 5, 6, 7}
};

//- Flanagan-Belytschko hourglass base vectors (xi eta, eta zeta, zeta xi,
//  xi eta zeta) in OpenFOAM hex node order
static const scalar hgBase[4][8] =
{
    { 1, -1,  1, -1,  1, -1,  1, -1},
    { 1,  1, -1, -1, -1, -1,  1,  1},
    { 1, -1, -1,  1, -1,  1,  1, -1},
    {-1,  1, -1,  1,  1, -1,  1, -1}
};

static const scalar faceXi[4] = {-1, 1, 1, -1};
static const scalar faceEta[4] = {-1, -1, 1, 1};

inline tensor hughesWinget(const tensor& dW)
{
    return inv(I - 0.5*dW) & (I + 0.5*dW);
}

} // End anonymous namespace


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::hex8Body::gradients
(
    const hexElement& e,
    const pointField& x,
    FixedList<vector, 8>& B
) const
{
    B = Zero;
    scalar V = 0;
    const scalar gp = 1/Foam::sqrt(3.0);
    for (label f = 0; f < 6; ++f)
    {
        point p[4];
        for (label k = 0; k < 4; ++k)
        {
            p[k] = x[e[hexFaces[f][k]]];
        }
        for (label a = 0; a < 2; ++a)
        {
            for (label b = 0; b < 2; ++b)
            {
                const scalar xi = (a ? gp : -gp);
                const scalar eta = (b ? gp : -gp);
                vector dxi(Zero), deta(Zero);
                point xg(Zero);
                scalar N[4];
                for (label k = 0; k < 4; ++k)
                {
                    N[k] = 0.25*(1 + faceXi[k]*xi)*(1 + faceEta[k]*eta);
                    dxi += 0.25*faceXi[k]*(1 + faceEta[k]*eta)*p[k];
                    deta += 0.25*faceEta[k]*(1 + faceXi[k]*xi)*p[k];
                    xg += N[k]*p[k];
                }
                // Unit Gauss weights
                const vector nA(dxi ^ deta);
                V += (xg & nA)/3;
                for (label k = 0; k < 4; ++k)
                {
                    B[hexFaces[f][k]] += N[k]*nA;
                }
            }
        }
    }
    for (label i = 0; i < 8; ++i)
    {
        B[i] /= V;
    }
    return V;
}


void Foam::hex8Body::internalForces(vectorField& f) const
{
    f.setSize(x_.size());
    f = Zero;
    forAll(elems_, ei)
    {
        const hexElement& e = elems_[ei];
        FixedList<vector, 8> B;
        const scalar V = gradients(e, x_, B);
        for (label i = 0; i < 8; ++i)
        {
            f[e[i]] -= V*(sigma_[ei] & B[i]);
        }

        // Hourglass forces: gamma_aI = h_aI - (sum_J h_aJ x_J) . B_I
        const FixedList<vector, 4>& Q = Qhg_[ei];
        for (label a = 0; a < 4; ++a)
        {
            vector hx(Zero);
            for (label j = 0; j < 8; ++j)
            {
                hx += hgBase[a][j]*x_[e[j]];
            }
            for (label i = 0; i < 8; ++i)
            {
                const scalar gamma = hgBase[a][i] - (hx & B[i]);
                f[e[i]] -= gamma*Q[a];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hex8Body::hex8Body
(
    const word& name,
    const pointField& points,
    const List<hexElement>& elems,
    const dictionary& dict
)
:
    name_(name),
    model_(constitutiveModel::New(dict.subDict("material"))),
    rho_(dict.get<scalar>("density")),
    hourglass_(dict.getOrDefault<scalar>("hourglass", 0.05)),
    planeStrain_(dict.getOrDefault<bool>("planeStrain", false)),
    elems_(elems),
    x_(points),
    v_(points.size(), Zero),
    m_(points.size(), 0),
    fixMask_(points.size(), Zero),
    fixVel_(points.size(), Zero),
    fExt_(points.size(), Zero),
    g_(Zero),
    sigma_(elems.size(), Zero),
    hist_(elems.size()),
    Qhg_(elems.size(), FixedList<vector, 4>(Zero)),
    V0_(elems.size(), 0),
    workInt_(0),
    workExt_(0)
{
    forAll(elems_, ei)
    {
        FixedList<vector, 8> B;
        scalar V = gradients(elems_[ei], x_, B);
        if (V < 0)
        {
            // Mirror: swap the bottom and top faces
            hexElement& e = elems_[ei];
            for (label k = 0; k < 4; ++k)
            {
                std::swap(e[k], e[k + 4]);
            }
            V = gradients(e, x_, B);
        }
        if (V <= 0)
        {
            FatalErrorInFunction
                << "Body " << name_ << ": element " << ei
                << " has no positive volume" << exit(FatalError);
        }
        V0_[ei] = V;
        for (label i = 0; i < 8; ++i)
        {
            m_[elems_[ei][i]] += rho_*V/8;
        }
    }
    if (planeStrain_)
    {
        forAll(fixMask_, i)
        {
            fixMask_[i].z() = 1;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hex8Body::setVelocity(const vectorField& v)
{
    v_ = v;
    forAll(v_, i)
    {
        for (direction d = 0; d < 3; ++d)
        {
            if (fixMask_[i][d] > 0)
            {
                v_[i][d] = fixVel_[i][d];
            }
        }
    }
}


void Foam::hex8Body::fix(const label i, const vector& mask, const vector& vel)
{
    for (direction d = 0; d < 3; ++d)
    {
        if (mask[d] > 0)
        {
            fixMask_[i][d] = 1;
            fixVel_[i][d] = vel[d];
            v_[i][d] = vel[d];
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::hex8Body::volumes() const
{
    auto tV = tmp<scalarField>::New(elems_.size());
    forAll(elems_, ei)
    {
        FixedList<vector, 8> B;
        tV.ref()[ei] = gradients(elems_[ei], x_, B);
    }
    return tV;
}


Foam::tmp<Foam::vectorField> Foam::hex8Body::forces() const
{
    auto tf = tmp<vectorField>::New();
    internalForces(tf.ref());
    return tf;
}


Foam::scalar Foam::hex8Body::stableDeltaT() const
{
    const scalar c = Foam::sqrt(model_->waveModulus()/rho_);
    scalar dt = GREAT;
    forAll(elems_, ei)
    {
        const hexElement& e = elems_[ei];
        FixedList<vector, 8> B;
        const scalar V = gradients(e, x_, B);
        scalar Amax = 0;
        for (label f = 0; f < 6; ++f)
        {
            const point& a = x_[e[hexFaces[f][0]]];
            const point& b = x_[e[hexFaces[f][1]]];
            const point& cc = x_[e[hexFaces[f][2]]];
            const point& d = x_[e[hexFaces[f][3]]];
            Amax = max(Amax, 0.5*mag((cc - a) ^ (d - b)));
        }
        dt = min(dt, V/Amax/c);
    }
    return dt;
}


void Foam::hex8Body::step(const scalar dt)
{
    // Forces at x_n from sigma_n
    vectorField f;
    internalForces(f);

    // Velocity to n + 1/2; work of the external forces and constraints
    forAll(v_, i)
    {
        const vector fTot(f[i] + fExt_[i] + m_[i]*g_);
        vector vNew(v_[i] + dt*fTot/m_[i]);
        for (direction d = 0; d < 3; ++d)
        {
            if (fixMask_[i][d] > 0)
            {
                vNew[d] = fixVel_[i][d];
            }
        }
        const vector vMid(0.5*(v_[i] + vNew));
        // Reaction of the constraints: m dv/dt - fTot
        const vector fReact(m_[i]*(vNew - v_[i])/dt - f[i]);
        workExt_ += dt*(fReact & vMid);
        workInt_ -= dt*(f[i] & vMid);
        v_[i] = vNew;
    }

    // Mid-step configuration, then positions to n + 1
    const pointField xMid(x_ + 0.5*dt*v_);
    x_ += dt*v_;

    const scalar M = model_->waveModulus();
    const bool rate = model_->rateDependent();
    forAll(elems_, ei)
    {
        const hexElement& e = elems_[ei];
        FixedList<vector, 8> B;
        const scalar V = gradients(e, xMid, B);

        tensor L(Zero);
        for (label i = 0; i < 8; ++i)
        {
            L += v_[e[i]]*B[i];
        }
        // L_ij = dv_i/dx_j (v (x) B)
        const symmTensor deps(dt*symm(L));
        const tensor dW(dt*skew(L));

        const symmTensor sRot(transform(hughesWinget(dW), sigma_[ei]));
        symmTensor sNew;
        constitutiveModel::history hNew;
        if (rate)
        {
            model_->update(deps, dt, sRot, hist_[ei], sNew, hNew);
        }
        else
        {
            model_->update(deps, sRot, hist_[ei], sNew, hNew);
        }
        sigma_[ei] = sNew;
        hist_[ei] = hNew;

        // Hourglass: integrated hourglass velocity, stiffness form
        scalar sumBB = 0;
        for (label i = 0; i < 8; ++i)
        {
            sumBB += magSqr(B[i]);
        }
        const scalar k = hourglass_*M*V*sumBB/8;
        for (label a = 0; a < 4; ++a)
        {
            vector hx(Zero);
            for (label j = 0; j < 8; ++j)
            {
                hx += hgBase[a][j]*xMid[e[j]];
            }
            vector q(Zero);
            for (label i = 0; i < 8; ++i)
            {
                q += (hgBase[a][i] - (hx & B[i]))*v_[e[i]];
            }
            Qhg_[ei][a] += dt*k*q;
        }
    }
}


Foam::scalar Foam::hex8Body::hourglassEnergy() const
{
    const scalar M = model_->waveModulus();
    scalar E = 0;
    forAll(elems_, ei)
    {
        FixedList<vector, 8> B;
        const scalar V = gradients(elems_[ei], x_, B);
        scalar sumBB = 0;
        for (label i = 0; i < 8; ++i)
        {
            sumBB += magSqr(B[i]);
        }
        const scalar k = hourglass_*M*V*sumBB/8;
        for (label a = 0; a < 4; ++a)
        {
            E += magSqr(Qhg_[ei][a])/(2*k);
        }
    }
    return E;
}


Foam::scalar Foam::hex8Body::kineticEnergy() const
{
    return 0.5*sum(m_*magSqr(v_));
}


// ************************************************************************* //
