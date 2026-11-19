/*---------------------------------------------------------------------------*\
    rigidBodies: the rigid bodies of a case (WP 4.1)
\*---------------------------------------------------------------------------*/

#include "rigidBodies.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OSspecific.H"
#include "IFstream.H"

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

    // Restart (WP 6.2, decision 0081): the state written at the start time
    const fileName stateFile(runTime.timePath()/"uniform"/"rigidBodies");
    const bool restart = isFile(stateFile);
    if (restart)
    {
        IFstream is(stateFile);
        const dictionary state(is);
        for (rigidBody& b : *this)
        {
            b.readState(state.subDict(b.name()));
        }
        Info<< "Rigid bodies: state read from " << runTime.timeName()
            << "/uniform/rigidBodies" << endl;
    }
    else if (runTime.value() > 0)
    {
        for (const rigidBody& b : *this)
        {
            if (b.free() || b.deformable())
            {
                FatalErrorInFunction
                    << "Body " << b.name() << " (free or deformable) cannot "
                    << "start at t = " << runTime.value() << " without its "
                    << "state: " << stateFile << " not found" << nl
                    << "Restart from a time that celFoam wrote (WP 6.2)"
                    << exit(FatalError);
            }
        }
    }

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
            // A restarted run appends to the history of the first run
            files_.set
            (
                i,
                new OFstream
                (
                    dir/(operator[](i).name() + ".dat"),
                    IOstreamOption(),
                    restart ? IOstreamOption::APPEND_APP : IOstreamOption::NO_APPEND
                )
            );
            OFstream& os = files_[i];
            os.precision(17);
            if (!restart)
            {
                os  << "# t  centre(3)  R(9, row by row)  v(3)  omega(3)" << nl;
            }
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


void Foam::rigidBodies::writeState() const
{
    if (empty())
    {
        return;
    }
    // Every processor writes its own copy (the bodies are whole on each)
    const fileName dir(runTime_.timePath()/"uniform");
    mkDir(dir);
    dictionary state;
    for (const rigidBody& b : *this)
    {
        dictionary bd;
        b.writeState(bd);
        state.set(b.name(), bd);
    }
    OFstream os(dir/"rigidBodies");
    os.precision(17);
    state.write(os, false);
}


void Foam::rigidBodies::write()
{
    writeState();

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
