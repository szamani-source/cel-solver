/*---------------------------------------------------------------------------*\
Class
    Foam::rigidBody

Description
    A rigid body whose surface is read from an STL (or any format triSurface
    reads) and whose motion is prescribed (WP 4.1, decision 0041).

    Kinematics: the body carries a reference point c (its "centre", default
    the STL's area-weighted centroid) and a rotation R from its pose at time
    zero, so a point X0 of the surface as read sits at

        x(t) = c(t) + R(t) (X0 - c0)

    and moves with the velocity v(t) + omega(t) ^ (x - c(t)). v and omega are
    Function1<vector>s of time (constant, table, sine, ...; default zero, a
    fixed wall). Over a step from t to t + dt:

        c += integral of v over the step                   (exact)
        R  = rotation by the angle |W| about W/|W|, then R  (exact when the
             axis of omega is fixed; W = integral of omega over the step)

    The STL is the pose at time zero. A run starting later is moved there in
    one step from zero by the same formulas.

    Dictionary entry (one per body, in constant/rigidBodies):

        insert
        {
            surface     insert.stl;           // in constant/triSurface
            centre      (0.2 0.5 0.05);       // optional, default centroid
            velocity    constant (0 0 0);     // optional, m/s
            omega       constant (0 0 0);     // optional, rad/s
            friction    20;                   // optional, deg, default 0:
                                              // wall friction angle of the
                                              // material on this body
                                              // (WP 4.3a, decision 0047)
        }

    Bodies act on the material as cut-cell walls (WP 4.3a, decision 0047);
    moving bodies push it and the cut cells follow them (WP 4.3b, decision
    0053, freeSurface yes). The motion is prescribed: the force of the
    material on the body is reported, not applied to it.

    Free bodies (WP 4.3d, decision 0058): with a mass the motion follows
    from the force of the material and gravity (the solid's g):

        block
        {
            surface     block.stl;
            mass        5.4;                  // kg, of the part the mesh
                                              // holds (2D: one cell deep)
            inertia     (0.01 0 0 0.01 0 0.01); // optional, kg m^2, about
                                              // the centre at time zero
                                              // (symmTensor); without it
                                              // the body does not rotate
            velocity    (0 0 0);              // optional, initial, m/s
            omega       (0 0 0);              // optional, initial, rad/s
        }

    Velocities are staggered in time (central difference): a step moves the
    body with the velocity the material saw in its Lagrange step, then
    changes the velocity by the step's impulse (accelerate). A free body
    cannot restart from a later time yet.

SourceFiles
    rigidBody.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_rigidBody_H
#define Foam_rigidBody_H

#include "triSurface.H"
#include "Function1.H"
#include "tensor.H"

namespace Foam
{

class Time;

class rigidBody
{
    // Private Data

        //- Name of the body (its dictionary keyword)
        word name_;

        //- Surface file name (relative to constant/triSurface)
        fileName surfaceFile_;

        //- Surface at time zero, as read
        triSurface surface0_;

        //- Reference point at time zero
        point c0_;

        //- Prescribed velocity of the reference point
        autoPtr<Function1<vector>> velocity_;

        //- Prescribed angular velocity
        autoPtr<Function1<vector>> omega_;

        //- Wall friction angle [deg]
        scalar friction_;

        //- Velocity or omega given, or free (a moving body)
        bool moving_;

        //- Mass [kg]; zero: prescribed motion
        scalar mass_;

        //- Inertia about the centre at time zero [kg m^2]; zero: no
        //  rotation
        symmTensor inertia0_;

        //- Free body: current velocity and angular velocity
        vector v_;
        vector w_;

        //- Current reference point
        point c_;

        //- Current rotation from the pose at time zero
        tensor R_;

        //- Time of the current pose
        scalar time_;


public:

    //- Runtime type information
    ClassName("rigidBody");


    // Constructors

        //- Construct from the body's dictionary; the surface is read from
        //  <case>/constant/triSurface (the undecomposed case in parallel;
        //  every processor holds the whole surface). Pose at time zero.
        rigidBody(const Time& runTime, const word& name, const dictionary& dict);


    // Member Functions

        //- Name
        const word& name() const noexcept { return name_; }

        //- Surface as read (pose at time zero)
        const triSurface& surface0() const noexcept { return surface0_; }

        //- Reference point at time zero
        const point& centre0() const noexcept { return c0_; }

        //- Current reference point
        const point& centre() const noexcept { return c_; }

        //- Current rotation from the pose at time zero
        const tensor& rotation() const noexcept { return R_; }

        //- Time of the current pose
        scalar time() const noexcept { return time_; }

        //- Wall friction angle [deg] of the material on this body
        scalar friction() const noexcept { return friction_; }

        //- Does the body move (prescribed motion, or free)?
        bool moving() const noexcept { return moving_; }

        //- Is the motion free (follows from the forces, WP 4.3d)?
        bool free() const noexcept { return mass_ > 0; }

        //- Mass [kg] (free bodies)
        scalar mass() const noexcept { return mass_; }

        //- Change a free body's velocity by the force F and moment M (about
        //  the current centre) over dt. Does nothing for prescribed bodies.
        void accelerate(const vector& F, const vector& M, const scalar dt);

        //- Velocity of the reference point at time t (free: current)
        vector velocity(const scalar t) const;

        //- Angular velocity at time t (free: current)
        vector omega(const scalar t) const;

        //- Current position of a point given in the time-zero pose
        point transform(const point& X0) const
        {
            return c_ + (R_ & (X0 - c0_));
        }

        //- Velocity of the body at the current position x (current time)
        vector velocityAt(const point& x) const;

        //- Current surface points
        tmp<pointField> points() const;

        //- Current surface (a copy of the time-zero surface, moved)
        triSurface surface() const;

        //- Move the pose from time() to time t
        void moveTo(const scalar t);
};


} // End namespace Foam

#endif
