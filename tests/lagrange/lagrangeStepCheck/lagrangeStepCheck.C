/*---------------------------------------------------------------------------*\
Application
    celLagrangeStepCheck

Description
    WP 2.3 test driver for eulerianSolid::lagrangeStep (run by
    tests/test_lagrange_step.sh). Sets an initial state by hand, runs Lagrange
    steps only (no remap) and writes the results for tests/lagrange/check.py.

    After each step the driver checks the contract of decision 0014
    (rho V_L = m with V_L = V + dt sum(phiL out)), then puts rho back to its
    value before the step: with no remap, this keeps the cell mass fixed
    (the "null remap" of linear acoustics), so total momentum sum(m U) can
    be compared from step to step.

    system/lagrangeCheckDict:

    \verbatim
        setup       wave;       // wave | rotation | shear

        // wave: plane wave in a box of side L with n cells per side
        mode        P;          // P | S
        k           (1 1 0);    // integer wave vector: phase 2 pi (k & x)/L
        polarisation (0 0 1);   // S only: made normal to k
        profile     sine;       // sine | step | pulse
        amplitude   1e-6;       // strain amplitude
        L           1;
        n           50;
        pulse       (0.5 0.7);  // pulse: (k^ & x)/L range with the wave
        crossings   1;          // run time = crossings L / c

        // rotation: v = omega ^ (x - centre), U held fixed (solveMomentum no)
        omega       (0 0 1);
        centre      (0.5 0.5 0);
        turns       1;
        steps       2000;       // per turn
        preStress   (xx xy xz yy yz zz);

        // shear: v = (gammaDot (y - centre.y), 0, 0), U held fixed
        gammaDot    1;
        gammaEnd    1;
        steps       2000;
        preStress   (...);
    \endverbatim

    Output: history.dat (master; per step: step, t, volume-mean sigma (6),
    largest deviation of sigma from the mean, max |U|, total momentum (3),
    sum m|U|, largest |rho V_L - m|/m), and cells.dat (cells.dat.procN in
    parallel): centre, U, sigma, epsPEq of every cell at the end.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"
#include "OFstream.H"
#include "mathematicalConstants.H"
#include "clockTime.H"

#include <cmath>

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("WP 2.3: Lagrange step checks (no remap)");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    eulerianSolid solid(mesh);

    IOdictionary dict
    (
        IOobject
        (
            "lagrangeCheckDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word setup(dict.get<word>("setup"));
    const scalarField& V = mesh.V().field();
    const vectorField& C = mesh.C().primitiveField();

    const scalar rho0 = gMin(solid.rho().primitiveField());
    const scalar M = solid.material().waveModulus();
    const scalar G = solid.material().shearModulus();
    const scalar K = M - 4.0*G/3.0;
    const scalar lambda = K - 2.0*G/3.0;

    volVectorField& U = solid.URef();
    volSymmTensorField& sigma = solid.sigmaRef();

    scalar dt = solid.stableDeltaT();
    label nSteps = 0;

    if (setup == "wave")
    {
        const word mode(dict.get<word>("mode"));
        const vector kInt(dict.get<vector>("k"));
        const word profile(dict.get<word>("profile"));
        const scalar A = dict.get<scalar>("amplitude");
        const scalar L = dict.get<scalar>("L");
        const label n = dict.get<label>("n");
        const scalar h = L/n;
        const vector khat(kInt/mag(kInt));

        vector p(khat);
        scalar c = Foam::sqrt(M/rho0);
        if (mode == "S")
        {
            p = dict.get<vector>("polarisation");
            p -= (p & khat)*khat;
            p /= mag(p);
            c = Foam::sqrt(G/rho0);
        }
        else if (mode != "P")
        {
            FatalIOErrorInFunction(dict) << "mode P or S" << exit(FatalIOError);
        }

        forAll(C, celli)
        {
            scalar e = 0;
            if (profile == "sine")
            {
                e = A*Foam::sin
                (
                    constant::mathematical::twoPi*(kInt & C[celli])/L
                );
            }
            else if (profile == "step")
            {
                // From integer cell indices, so every cell of a diagonal is
                // on the same side (decision 0013)
                label s = 0;
                for (direction d = 0; d < 3; ++d)
                {
                    s += label(kInt[d])*label(std::floor(C[celli][d]/h));
                }
                s = ((s % n) + n) % n;
                e = (2*s < n) ? A : 0;
            }
            else if (profile == "pulse")
            {
                const scalarList r(dict.get<scalarList>("pulse"));
                const scalar x = (khat & C[celli])/L;
                e = (x >= r[0] && x < r[1]) ? A : 0;
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "profile sine, step or pulse" << exit(FatalIOError);
            }

            // Right-going simple wave: u = f(s - c t) p, strain e = f'
            U[celli] = -c*e*p;
            if (mode == "P")
            {
                sigma[celli] = lambda*e*symmTensor::I + 2*G*e*sqr(khat);
            }
            else
            {
                sigma[celli] = G*e*twoSymm(khat*p);
            }
        }

        const scalar tEnd = dict.get<scalar>("crossings")*L/c;
        nSteps = label(std::ceil(tEnd/dt));
        dt = tEnd/nSteps;
        Info<< "Wave " << mode << " k " << kInt << " " << profile
            << ", c = " << c << ", " << nSteps << " steps" << nl;
    }
    else if (setup == "rotation" || setup == "shear")
    {
        if (solid.solveMomentum())
        {
            FatalErrorInFunction
                << "setup " << setup << " needs solveMomentum no"
                << exit(FatalError);
        }
        const vector centre(dict.get<vector>("centre"));
        const symmTensor s0(dict.get<symmTensor>("preStress"));

        auto velocity = [&](const vector& x)
        {
            if (setup == "rotation")
            {
                return vector(dict.get<vector>("omega") ^ (x - centre));
            }
            return vector
            (
                dict.get<scalar>("gammaDot")*(x.y() - centre.y()), 0, 0
            );
        };

        forAll(C, celli)
        {
            U[celli] = velocity(C[celli]);
            sigma[celli] = s0;
        }
        forAll(U.boundaryField(), patchi)
        {
            fvPatchVectorField& pU = U.boundaryFieldRef()[patchi];
            if (pU.fixesValue())
            {
                const vectorField& Cf = mesh.C().boundaryField()[patchi];
                vectorField v(Cf.size());
                forAll(Cf, i)
                {
                    v[i] = velocity(Cf[i]);
                }
                pU == v;
            }
        }

        scalar tEnd = 0;
        if (setup == "rotation")
        {
            tEnd = dict.get<scalar>("turns")*constant::mathematical::twoPi
                /mag(dict.get<vector>("omega"));
            nSteps = label
            (
                std::round(dict.get<scalar>("turns")*dict.get<label>("steps"))
            );
        }
        else
        {
            tEnd = dict.get<scalar>("gammaEnd")/dict.get<scalar>("gammaDot");
            nSteps = dict.get<label>("steps");
        }
        dt = tEnd/nSteps;
        Info<< "Prescribed " << setup << ", " << nSteps << " steps" << nl;
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "setup wave, rotation or shear" << exit(FatalIOError);
    }

    U.correctBoundaryConditions();
    sigma.correctBoundaryConditions();

    Info<< "deltaT = " << dt << " (stable " << solid.stableDeltaT() << ")"
        << nl << endl;

    autoPtr<OFstream> hist;
    if (Pstream::master())
    {
        hist.reset(new OFstream("history.dat"));
        hist().precision(17);
        hist()
            << "# step t s_xx s_xy s_xz s_yy s_yz s_zz spread maxU "
            << "Px Py Pz sum_m|U| rhoErr" << nl;
    }

    const scalar Vtot = gSum(V);

    auto record = [&](const label step, const scalar rhoErr)
    {
        const scalarField m(solid.rho().primitiveField()*V);
        const symmTensor sMean
        (
            gSum(solid.sigma().primitiveField()*V)/Vtot
        );
        const scalar spread =
            gMax(mag(solid.sigma().primitiveField() - sMean));
        const vector P(gSum(m*solid.U().primitiveField()));
        const scalar Pscale = gSum(m*mag(solid.U().primitiveField()));
        const scalar maxU = gMax(mag(solid.U().primitiveField()));
        if (Pstream::master())
        {
            hist()
                << step << ' ' << step*dt;
            for (direction d = 0; d < 6; ++d)
            {
                hist() << ' ' << sMean[d];
            }
            hist()
                << ' ' << spread << ' ' << maxU
                << ' ' << P.x() << ' ' << P.y() << ' ' << P.z()
                << ' ' << Pscale << ' ' << rhoErr << nl;
        }
    };

    record(0, 0);

    const Foam::clockTime startClock;
    for (label step = 1; step <= nSteps; ++step)
    {
        const scalarField rhoOld(solid.rho().primitiveField());
        const scalarField m(rhoOld*V);

        solid.lagrangeStep(dt);

        // Contract (decision 0014): rho V_L = m
        const scalarField VL
        (
            V + dt*V*fvc::surfaceIntegrate(solid.phiL())().primitiveField()
        );
        const scalar rhoErr =
            gMax(mag(solid.rho().primitiveField()*VL - m)/m);

        // Null remap: keep the cell mass
        solid.rhoRef().primitiveFieldRef() = rhoOld;
        solid.rhoRef().correctBoundaryConditions();

        record(step, rhoErr);
    }
    const scalar elapsed = startClock.elapsedTime();

    Info<< "Run time per cell per step = "
        << elapsed/max(nSteps, 1)/returnReduce(mesh.nCells(), sumOp<label>())
          *Pstream::nProcs()
        << " s (wall time x processors)" << nl;

    {
        const fileName out
        (
            Pstream::parRun()
          ? fileName("cells.dat.proc" + Foam::name(Pstream::myProcNo()))
          : fileName("cells.dat")
        );
        OFstream os(runTime.globalPath()/out);
        os.precision(17);
        os  << "# x y z Ux Uy Uz s_xx s_xy s_xz s_yy s_yz s_zz epsPEq" << nl;
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
            os << ' ' << solid.epsPEq()[celli] << nl;
        }
    }

    Info<< "End" << nl << endl;
    return 0;
}

// ************************************************************************* //
