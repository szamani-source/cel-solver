/*---------------------------------------------------------------------------*\
    rigidBody: STL surface with prescribed rigid motion (WP 4.1)
\*---------------------------------------------------------------------------*/

#include "rigidBody.H"
#include "Time.H"

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
    surfaceFile_(dict.get<fileName>("surface")),
    surface0_
    (
        runTime.globalPath()/runTime.constant()/"triSurface"
       /surfaceFile_.expand()
    ),
    c0_(dict.getOrDefault<point>("centre", centroid(surface0_))),
    velocity_(functionOrZero("velocity", dict)),
    omega_(functionOrZero("omega", dict)),
    friction_(dict.getOrDefault<scalar>("friction", 0)),
    moving_(dict.found("velocity") || dict.found("omega")),
    c_(c0_),
    R_(tensor::I),
    time_(0)
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
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::rigidBody::velocity(const scalar t) const
{
    return velocity_->value(t);
}


Foam::vector Foam::rigidBody::omega(const scalar t) const
{
    return omega_->value(t);
}


Foam::vector Foam::rigidBody::velocityAt(const point& x) const
{
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


void Foam::rigidBody::moveTo(const scalar t)
{
    if (t == time_)
    {
        return;
    }
    c_ += velocity_->integrate(time_, t);
    R_ = rodrigues(omega_->integrate(time_, t)) & R_;
    time_ = t;
}


// ************************************************************************* //
