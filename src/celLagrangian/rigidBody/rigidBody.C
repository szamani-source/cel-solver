/*---------------------------------------------------------------------------*\
    rigidBody: STL surface with prescribed rigid motion (WP 4.1)
\*---------------------------------------------------------------------------*/

#include "rigidBody.H"
#include "Time.H"
#include "triangle.H"
#include "boundBox.H"
#include "HashTable.H"

#include <algorithm>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBody, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

//- Rotation by the angle |W| about the axis W/|W| (Rodrigues)
tensor rodrigues(const vector& W)
{
    const scalar theta = mag(W);
    if (theta < VSMALL)
    {
        return tensor::I;
    }
    const vector k(W/theta);
    const tensor K
    (
        0, -k.z(), k.y(),
        k.z(), 0, -k.x(),
        -k.y(), k.x(), 0
    );
    return tensor::I + Foam::sin(theta)*K + (1 - Foam::cos(theta))*(K & K);
}


//- Area-weighted centroid of a triangulated surface
point centroid(const triSurface& s)
{
    const pointField& p = s.points();
    vector sum(Zero);
    scalar area = 0;
    for (const labelledTri& f : s)
    {
        const scalar a = f.mag(p);
        sum += a*f.centre(p);
        area += a;
    }
    return area > VSMALL ? sum/area : point(Zero);
}


//- Function1 entry or zero when absent
autoPtr<Function1<vector>> functionOrZero
(
    const word& key,
    const dictionary& dict
)
{
    if (dict.found(key))
    {
        return Function1<vector>::New(key, dict);
    }
    dictionary d;
    d.add(key, vector::zero);
    return Function1<vector>::New(key, d);
}

//- Structured box of hex8 elements from a deformable entry
void boxMesh
(
    const dictionary& d,
    pointField& p,
    List<hex8Body::hexElement>& e
)
{
    const point lo(d.get<point>("min"));
    const point hi(d.get<point>("max"));
    const Vector<label> n(d.get<Vector<label>>("cells"));
    const label nx = n.x(), ny = n.y(), nz = n.z();
    auto id = [&](label i, label j, label k)
    {
        return i + (nx + 1)*(j + (ny + 1)*k);
    };
    p.setSize((nx + 1)*(ny + 1)*(nz + 1));
    for (label k = 0; k <= nz; ++k)
        for (label j = 0; j <= ny; ++j)
            for (label i = 0; i <= nx; ++i)
            {
                p[id(i, j, k)] = point
                (
                    lo.x() + (hi.x() - lo.x())*i/nx,
                    lo.y() + (hi.y() - lo.y())*j/ny,
                    lo.z() + (hi.z() - lo.z())*k/nz
                );
            }
    e.setSize(nx*ny*nz);
    label m = 0;
    for (label k = 0; k < nz; ++k)
        for (label j = 0; j < ny; ++j)
            for (label i = 0; i < nx; ++i)
            {
                hex8Body::hexElement& c = e[m++];
                c[0] = id(i, j, k);         c[1] = id(i + 1, j, k);
                c[2] = id(i + 1, j + 1, k); c[3] = id(i, j + 1, k);
                c[4] = id(i, j, k + 1);     c[5] = id(i + 1, j, k + 1);
                c[6] = id(i + 1, j + 1, k + 1); c[7] = id(i, j + 1, k + 1);
            }
}


//- Boundary faces of hex elements as an outward triangulated surface;
//  surfNode: the element node of each surface point
triSurface boundarySurface
(
    const pointField& p,
    const List<hex8Body::hexElement>& e,
    labelList& surfNode
)
{
    static const label hf[6][4] =
    {
        {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
        {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}
    };
    HashTable<label, FixedList<label, 4>> count;
    auto key = [](FixedList<label, 4> k)
    {
        std::sort(k.begin(), k.end());
        return k;
    };
    for (const auto& c : e)
    {
        for (label f = 0; f < 6; ++f)
        {
            FixedList<label, 4> q;
            for (label k = 0; k < 4; ++k) q[k] = c[hf[f][k]];
            count(key(q)) += 1;
        }
    }
    labelList local(p.size(), -1);
    DynamicList<label> nodes;
    DynamicList<labelledTri> tris;
    for (const auto& c : e)
    {
        for (label f = 0; f < 6; ++f)
        {
            FixedList<label, 4> q;
            for (label k = 0; k < 4; ++k) q[k] = c[hf[f][k]];
            if (count[key(q)] != 1) continue;
            for (label k = 0; k < 4; ++k)
            {
                if (local[q[k]] < 0)
                {
                    local[q[k]] = nodes.size();
                    nodes.append(q[k]);
                }
            }
            tris.append(labelledTri(local[q[0]], local[q[1]], local[q[2]], 0));
            tris.append(labelledTri(local[q[0]], local[q[2]], local[q[3]], 0));
        }
    }
    surfNode = nodes;
    return triSurface(tris, pointField(p, surfNode));
}


//- The surface: the STL, or the boundary of a deformable box
triSurface initialSurface(const Time& runTime, const dictionary& dict)
{
    if (dict.found("deformable"))
    {
        pointField p;
        List<hex8Body::hexElement> e;
        boxMesh(dict.subDict("deformable"), p, e);
        labelList sn;
        return boundarySurface(p, e, sn);
    }
    return triSurface
    (
        runTime.globalPath()/runTime.constant()/"triSurface"
       /dict.get<fileName>("surface").expand()
    );
}


//- Nearest point of the triangle (a, b, c) to x and its barycentric
//  weights
scalar nearestOnTri
(
    const point& a, const point& b, const point& c, const point& x,
    FixedList<scalar, 3>& w
)
{
    const point y(triPointRef(a, b, c).nearestPoint(x).point());
    const vector e0(b - a), e1(c - a), d(y - a);
    const scalar d00 = e0 & e0, d01 = e0 & e1, d11 = e1 & e1;
    const scalar d20 = d & e0, d21 = d & e1;
    const scalar den = d00*d11 - d01*d01;
    const scalar v = den > VSMALL ? (d11*d20 - d01*d21)/den : 0;
    const scalar u = den > VSMALL ? (d00*d21 - d01*d20)/den : 0;
    w[1] = v; w[2] = u; w[0] = 1 - u - v;
    return magSqr(y - x);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBody::rigidBody
(
    const Time& runTime,
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    surfaceFile_(dict.getOrDefault<fileName>("surface", fileName::null)),
    surface0_(initialSurface(runTime, dict)),
    c0_(dict.getOrDefault<point>("centre", centroid(surface0_))),
    velocity_(functionOrZero("velocity", dict)),
    omega_(functionOrZero("omega", dict)),
    friction_(dict.getOrDefault<scalar>("friction", 0)),
    adhesion_(dict.getOrDefault<scalar>("adhesion", 0)),
    moving_(dict.found("velocity") || dict.found("omega")),
    mass_(dict.getOrDefault<scalar>("mass", 0)),
    inertia0_(dict.getOrDefault<symmTensor>("inertia", symmTensor::zero)),
    v_(Zero),
    w_(Zero),
    c_(c0_),
    R_(tensor::I),
    time_(0),
    tHex_(0)
{
    if (surface0_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Body " << name_ << ": surface " << surfaceFile_
            << " has no triangles" << exit(FatalIOError);
    }
    if (friction_ < 0 || friction_ >= 90)
    {
        FatalIOErrorInFunction(dict)
            << "Body " << name_ << ": friction " << friction_
            << " deg must be in [0, 90)" << exit(FatalIOError);
    }
    if (adhesion_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Body " << name_ << ": adhesion " << adhesion_
            << " Pa must not be negative" << exit(FatalIOError);
    }
    if (dict.found("deformable"))
    {
        // Deformable body (WP 5.3, decision 0072): fixed frame, the surface
        // moves with the elements
        const dictionary& d = dict.subDict("deformable");
        pointField p;
        List<hex8Body::hexElement> e;
        boxMesh(d, p, e);
        boundarySurface(p, e, surfNode_);
        hex_.reset(new hex8Body(name_, p, e, d));
        nodeForce_.setSize(p.size(), Zero);
        surfU_.setSize(surfNode_.size(), Zero);
        moving_ = true;
        if (d.found("drive"))
        {
            const boundBox bb(d.get<boundBox>("drive"));
            DynamicList<label> dn;
            forAll(p, i)
            {
                if (bb.contains(p[i])) dn.append(i);
            }
            driveNodes_ = dn;
            const vector v0(velocity_->value(0));
            for (const label i : driveNodes_)
            {
                hex_->fix(i, vector::one, v0);
            }
        }
        surfU_ = vectorField(hex_->velocity(), surfNode_);
        if (mass_ > 0 || dict.found("inertia"))
        {
            FatalIOErrorInFunction(dict)
                << "Body " << name_ << ": a deformable body takes its mass "
                << "from its density" << exit(FatalIOError);
        }
    }
    if (mass_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Body " << name_ << ": mass " << mass_ << " must be positive"
            << exit(FatalIOError);
    }
    if (mass_ > 0)
    {
        // Free body (WP 4.3d): velocity and omega are initial values
        moving_ = true;
        v_ = velocity_->value(0);
        w_ = omega_->value(0);
        if (runTime.value() > 0)
        {
            FatalIOErrorInFunction(dict)
                << "Body " << name_ << ": a free body (mass) cannot restart"
                << " from t = " << runTime.value() << " yet"
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::rigidBody::velocity(const scalar t) const
{
    return free() ? v_ : velocity_->value(t);
}


Foam::vector Foam::rigidBody::omega(const scalar t) const
{
    return free() ? w_ : omega_->value(t);
}


void Foam::rigidBody::accelerate
(
    const vector& F,
    const vector& M,
    const scalar dt
)
{
    if (!free())
    {
        return;
    }
    v_ += dt*F/mass_;
    if (mag(inertia0_) > VSMALL)
    {
        // Inertia in the current pose; Euler's equations in the lab frame
        const tensor J(R_ & tensor(inertia0_) & R_.T());
        w_ += dt*(inv(J) & (M - (w_ ^ (J & w_))));
    }
}


Foam::vector Foam::rigidBody::velocityAt(const point& x) const
{
    if (hex_)
    {
        // Nearest surface triangle, its points' velocities over the step
        const pointField& sp = surface0_.points();
        const vectorField& v = surfU_;
        scalar best = GREAT;
        vector u(Zero);
        for (const labelledTri& t : surface0_)
        {
            FixedList<scalar, 3> w;
            const scalar d2 = nearestOnTri(sp[t[0]], sp[t[1]], sp[t[2]], x, w);
            if (d2 < best)
            {
                best = d2;
                u = w[0]*v[t[0]] + w[1]*v[t[1]] + w[2]*v[t[2]];
            }
        }
        return u;
    }
    return velocity(time_) + (omega(time_) ^ (x - c_));
}


Foam::tmp<Foam::pointField> Foam::rigidBody::points() const
{
    const pointField& p0 = surface0_.points();
    auto tp = tmp<pointField>::New(p0.size());
    pointField& p = tp.ref();
    forAll(p0, i)
    {
        p[i] = transform(p0[i]);
    }
    return tp;
}


Foam::triSurface Foam::rigidBody::surface() const
{
    triSurface s(surface0_);
    s.movePoints(points());
    return s;
}


void Foam::rigidBody::addSurfaceForce(const point& x, const vector& f)
{
    const pointField& sp = surface0_.points();
    scalar best = GREAT;
    label tBest = -1;
    FixedList<scalar, 3> wBest(Zero);
    forAll(surface0_, ti)
    {
        const labelledTri& t = surface0_[ti];
        FixedList<scalar, 3> w;
        const scalar d2 = nearestOnTri(sp[t[0]], sp[t[1]], sp[t[2]], x, w);
        if (d2 < best)
        {
            best = d2;
            tBest = ti;
            wBest = w;
        }
    }
    if (tBest >= 0)
    {
        const labelledTri& t = surface0_[tBest];
        for (label k = 0; k < 3; ++k)
        {
            nodeForce_[surfNode_[t[k]]] += wBest[k]*f;
        }
    }
}


void Foam::rigidBody::advance
(
    const scalar t,
    const scalar dtNext,
    const vector& g
)
{
    if (!hex_)
    {
        return;
    }

    // The surface the material sees moves with the velocity it saw
    if (t > time_)
    {
        pointField sp(surface0_.points());
        sp += (t - time_)*surfU_;
        surface0_.movePoints(sp);
        time_ = t;
    }

    // Elements one step ahead, with this step's forces
    hex8Body& h = hex_();
    const scalar tTo = t + dtNext;
    if (tTo > tHex_)
    {
        h.setGravity(g);
        h.externalForce() = nodeForce_;
        const scalar span = tTo - tHex_;
        const label n =
            max(label(std::ceil(span/(0.5*h.stableDeltaT()))), label(1));
        const scalar dt = span/n;
        for (label k = 0; k < n; ++k)
        {
            const vector vd(velocity_->value(tHex_ + (k + 1)*dt));
            for (const label i : driveNodes_)
            {
                h.fix(i, vector::one, vd);
            }
            h.step(dt);
        }
        tHex_ = tTo;
    }
    nodeForce_ = Zero;
    if (dtNext > 0)
    {
        surfU_ =
            (pointField(h.points(), surfNode_) - surface0_.points())/dtNext;
    }
}


void Foam::rigidBody::moveTo(const scalar t)
{
    if (t == time_ || hex_)
    {
        // Deformable bodies move in advance()
        return;
    }
    if (free())
    {
        // Free body (WP 4.3d): the current velocity over the step
        c_ += (t - time_)*v_;
        R_ = rodrigues((t - time_)*w_) & R_;
    }
    else
    {
        c_ += velocity_->integrate(time_, t);
        R_ = rodrigues(omega_->integrate(time_, t)) & R_;
    }
    time_ = t;
}


// ************************************************************************* //
