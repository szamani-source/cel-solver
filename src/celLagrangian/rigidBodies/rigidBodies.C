/*---------------------------------------------------------------------------*\
    rigidBodies: the rigid bodies of a case (WP 4.1)
\*---------------------------------------------------------------------------*/

#include "rigidBodies.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodies::rigidBodies(const Time& runTime)
:
    PtrList<rigidBody>(),
    runTime_(runTime),
    files_()
{
    IOobject io
    (
        "rigidBodies",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return;
    }

    const IOdictionary dict(io);
    const dictionary& bodies = dict.subDict("bodies");

    label n = 0;
    for (const entry& e : bodies)
    {
        if (e.isDict())
        {
            ++n;
        }
    }
    setSize(n);

    n = 0;
    for (const entry& e : bodies)
    {
        if (e.isDict())
        {
            set(n++, new rigidBody(runTime, e.keyword(), e.dict()));
        }
    }

    move();

    const label oldPrecision = Info().precision(16);
    for (const rigidBody& b : *this)
    {
        Info<< "Rigid body " << b.name() << ": "
            << b.surface0().size() << " triangles, "
            << b.surface0().nPoints() << " points, centre "
            << b.centre0() << endl;
    }
    Info().precision(oldPrecision);

    if (Pstream::master() && size())
    {
        const fileName dir(runTime.globalPath()/"postProcessing"/"rigidBodies");
        mkDir(dir);
        files_.setSize(size());
        forAll(*this, i)
        {
            files_.set(i, new OFstream(dir/(operator[](i).name() + ".dat")));
            OFstream& os = files_[i];
            os.precision(17);
            os  << "# t  centre(3)  R(9, row by row)  v(3)  omega(3)" << nl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rigidBodies::move()
{
    const scalar t = runTime_.value();
    for (rigidBody& b : *this)
    {
        b.moveTo(t);
    }
}


void Foam::rigidBodies::write()
{
    if (!Pstream::master() || empty())
    {
        return;
    }

    const scalar t = runTime_.value();
    const fileName dir
    (
        runTime_.globalPath()/"postProcessing"/"rigidBodies"/runTime_.timeName()
    );
    mkDir(dir);

    forAll(*this, i)
    {
        const rigidBody& b = operator[](i);
        OFstream& os = files_[i];
        const point& c = b.centre();
        const tensor& R = b.rotation();
        const vector v(b.velocity(t));
        const vector w(b.omega(t));

        os  << t << ' ' << c.x() << ' ' << c.y() << ' ' << c.z();
        for (direction k = 0; k < tensor::nComponents; ++k)
        {
            os  << ' ' << R[k];
        }
        os  << ' ' << v.x() << ' ' << v.y() << ' ' << v.z()
            << ' ' << w.x() << ' ' << w.y() << ' ' << w.z() << nl;
        os.flush();

        b.surface().write(dir/(b.name() + ".stl"));
    }
}


// ************************************************************************* //
