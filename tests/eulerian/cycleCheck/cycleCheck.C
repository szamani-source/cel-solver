/*---------------------------------------------------------------------------*\
Application
    celCycleCheck

Description
    WP 2.4b test driver for the full explicit cycle of eulerianSolid
    (lagrangeStep then remap, as celFoam runs them). Run by
    tests/test_eulerian_cycle.sh through tests/eulerian/cycle_checks.py.

    The initial state (U, sigma, optionally rho) is read from the 0
    directory, written by the Python script. The driver runs a fixed number
    of cycles with a fixed time step and records, every step, the totals and
    the range of each field in a wave frame, so the script can check
    conservation, speeds and extrema without writing every time step.

    system/cycleCheckDict:

    \verbatim
        deltaT      1e-5;           // time step [s]
        nSteps      400;            // cycles
        holdU       no;             // yes: U reset to its initial field
                                    // after every cycle (prescribed motion,
                                    // with solveMomentum no in controlDict)
        frameNormal (1 0 0);        // wave normal n (in plane); t = z ^ n
        streaming   (0 0 0);        // subtracted from U in the frame values
        prescribe   none;           // none | rotation: U = omega ^ (x - centre)
        omega       (0 0 1);        // in the cells and on fixedValue patches
        centre      (0.5 0.5 0);    // (overrides 0/U)
        writeSteps  (200 400);      // write cells_<step>.dat at these steps
    \endverbatim

    history.dat (master), one line per step:
        step t mass Px Py Pz sum(m|U|) sum(m epsPEq)
        then min and max of vn vt snn stt szz snt tau mean
        (tau = (snn - stt)/2, mean = (snn + stt)/2; velocities relative to
        the streaming velocity) and the wall time of the Lagrange step and
        of the remap for that step.

    cells_<step>.dat (cells_<step>.dat.procN in parallel), one line per cell:
        x y z Ux Uy Uz s_xx s_xy s_xz s_yy s_yz s_zz epsPEq rho
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "OFstream.H"
#include "clockTime.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("WP 2.4b: Lagrange step plus remap checks");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    eulerianSolid solid(mesh);

    IOdictionary dict
    (
        IOobject
        (
            "cycleCheckDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const scalar dt = dict.get<scalar>("deltaT");
    const label nSteps = dict.get<label>("nSteps");
    const bool holdU = dict.getOrDefault<Switch>("holdU", false);
    vector nHat(dict.getOrDefault<vector>("frameNormal", vector(1, 0, 0)));
    nHat /= mag(nHat);
    const vector tHat(vector(0, 0, 1) ^ nHat);
    const vector streaming(dict.getOrDefault<vector>("streaming", Zero));
    const labelHashSet writeSteps
    (
        dict.getOrDefault<labelList>("writeSteps", labelList())
    );

    const scalarField& V = mesh.V().field();
    const vectorField& C = mesh.C().primitiveField();

    if (holdU && solid.solveMomentum())
    {
        FatalErrorInFunction
            << "holdU yes needs solveMomentum no in controlDict"
            << exit(FatalError);
    }

    // Optional prescribed rigid rotation v = omega ^ (x - centre), set in
    // the cells and on fixedValue patches (their face centres)
    const word prescribe(dict.getOrDefault<word>("prescribe", "none"));
    if (prescribe == "rotation")
    {
        const vector omega(dict.get<vector>("omega"));
        const vector centre(dict.get<vector>("centre"));
        volVectorField& U = solid.URef();
        U.primitiveFieldRef() = omega ^ (C - centre);
        forAll(U.boundaryField(), patchi)
        {
            fvPatchVectorField& pU = U.boundaryFieldRef()[patchi];
            if (pU.fixesValue())
            {
                pU == (omega ^ (mesh.C().boundaryField()[patchi] - centre));
            }
        }
    }
    else if (prescribe != "none")
    {
        FatalIOErrorInFunction(dict)
            << "prescribe none or rotation" << exit(FatalIOError);
    }

    solid.URef().correctBoundaryConditions();
    solid.sigmaRef().correctBoundaryConditions();
    const volVectorField U0("U0", solid.U());

    Info<< "deltaT = " << dt << " (stable " << solid.stableDeltaT() << "), "
        << nSteps << " steps" << nl << endl;

    autoPtr<OFstream> hist;
    if (Pstream::master())
    {
        hist.reset(new OFstream("history.dat"));
        hist().precision(17);
        hist()
            << "# step t mass Px Py Pz sum_m|U| sum_m_epsPEq"
            << " vn(min max) vt snn stt szz snt tau mean tLagrange tRemap"
            << nl;
    }

    auto record = [&](const label step, const scalar tL, const scalar tR)
    {
        const scalarField m(solid.rho().primitiveField()*V);
        const vectorField& U = solid.U().primitiveField();
        const symmTensorField& s = solid.sigma().primitiveField();

        const scalar mass = gSum(m);
        const vector P(gSum(m*U));
        const scalar Pscale = gSum(m*mag(U));
        const scalar ep = gSum(m*solid.epsPEq().primitiveField());

        const label nq = 8;
        List<scalarField> q(nq, scalarField(U.size()));
        forAll(U, celli)
        {
            const vector u(U[celli] - streaming);
            const vector tn(s[celli] & nHat);
            const scalar snn = tn & nHat;
            const scalar stt = tHat & (s[celli] & tHat);
            q[0][celli] = u & nHat;
            q[1][celli] = u & tHat;
            q[2][celli] = snn;
            q[3][celli] = stt;
            q[4][celli] = s[celli].zz();
            q[5][celli] = tn & tHat;
            q[6][celli] = 0.5*(snn - stt);
            q[7][celli] = 0.5*(snn + stt);
        }
        scalarList lo(nq), hi(nq);
        for (label k = 0; k < nq; ++k)
        {
            lo[k] = gMin(q[k]);
            hi[k] = gMax(q[k]);
        }

        if (Pstream::master())
        {
            OFstream& os = hist();
            os  << step << ' ' << step*dt << ' ' << mass
                << ' ' << P.x() << ' ' << P.y() << ' ' << P.z()
                << ' ' << Pscale << ' ' << ep;
            for (label k = 0; k < nq; ++k)
            {
                os << ' ' << lo[k] << ' ' << hi[k];
            }
            os  << ' ' << tL << ' ' << tR << nl;
        }
    };

    auto writeCells = [&](const label step)
    {
        const word base("cells_" + Foam::name(step) + ".dat");
        const fileName out
        (
            Pstream::parRun()
          ? fileName(base + ".proc" + Foam::name(Pstream::myProcNo()))
          : fileName(base)
        );
        OFstream os(runTime.globalPath()/out);
        os.precision(17);
        os  << "# x y z Ux Uy Uz s_xx s_xy s_xz s_yy s_yz s_zz epsPEq rho"
            << nl;
        forAll(C, celli)
        {
            const vector& u = solid.U()[celli];
            const symmTensor& s = solid.sigma()[celli];
            os  << C[celli].x() << ' ' << C[celli].y() << ' ' << C[celli].z()
                << ' ' << u.x() << ' ' << u.y() << ' ' << u.z();
            for (direction d = 0; d < 6; ++d)
            {
                os << ' ' << s[d];
            }
            os  << ' ' << solid.epsPEq()[celli]
                << ' ' << solid.rho()[celli] << nl;
        }
    };

    record(0, 0, 0);
    if (writeSteps.found(0))
    {
        writeCells(0);
    }

    scalar tLagrange = 0, tRemap = 0;
    for (label step = 1; step <= nSteps; ++step)
    {
        const Foam::clockTime c1;
        solid.lagrangeStep(dt);
        const scalar tL = c1.elapsedTime();

        const Foam::clockTime c2;
        solid.remap(dt);
        const scalar tR = c2.elapsedTime();

        tLagrange += tL;
        tRemap += tR;

        if (holdU)
        {
            solid.URef() == U0;
        }

        record(step, tL, tR);
        if (writeSteps.found(step))
        {
            writeCells(step);
        }
    }

    const scalar cellSteps =
        scalar(max(nSteps, 1))*returnReduce(mesh.nCells(), sumOp<label>())
       /Pstream::nProcs();
    Info<< "Lagrange step per cell per step = " << tLagrange/cellSteps
        << " s" << nl
        << "Remap per cell per step = " << tRemap/cellSteps << " s" << nl
        << "Cycle per cell per step = " << (tLagrange + tRemap)/cellSteps
        << " s (wall time x processors)" << nl;

    Info<< "End" << nl << endl;
    return 0;
}

// ************************************************************************* //
