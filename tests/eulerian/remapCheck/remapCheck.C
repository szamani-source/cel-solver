/*---------------------------------------------------------------------------*\
Application
    celRemapCheck

Description
    WP 2.4: drives eulerianSolid::remap alone with a prescribed motion (no
    Lagrange step). Run by tests/test_remap.sh through
    tests/eulerian/remap_test.py, which writes the mesh, the initial fields
    (0/rho, U, sigma, epsP, epsPEq) and system/remapCheckDict, and checks
    what this program prints and writes.

    Each step: phiL from the motion, the Lagrangian density consistent with
    it (rho_L = rho V/V_L, as the Lagrange step would leave it), remap(dt).

    system/remapCheckDict:
    \verbatim
        motion      translate;  // translate | rotate | swirl | compress
        U           (1 1 0);    // translate: velocity [m/s]
        omega       1;          // rotate: angular velocity [1/s] about z
        centre      (0.5 0.5 0);// rotate
        rIn         0.4;        // rotate: rigid rotation inside rIn,
        rOut        0.48;       //   cosine taper to rest at rOut
        A           0.1;        // swirl: psi = A sin^2(pi x/L) sin^2(pi y/L)
                                // compress: u = (A Foam::sin(2 pi x/L), 0, 0)
        L           1;          // swirl, compress: period [m]
        reverse     false;      // flip the sign of phiL every step
        endTime     1;          // [s]
        courant     0.25;       // material Courant number umax*dt/h
        nSteps      0;          // if > 0: fixed number of steps instead
    \endverbatim

    Prints (one line each, parsed by remap_test.py):
        steps N dt X
        maxCourantSum X        largest dt*sum(outflow)/V over cells and steps
        conservation X         largest relative change over all steps of
                               mass, momentum and sum(m q) for every field
                               component (scaled by sum(m |q|))
        range <field> <cmpt> min0 max0 minRun maxRun
    and writes the fields at endTime with 17 digits.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "eulerianSolid.H"

using namespace Foam;

// Stream function of the rotation: u_theta = omega r for r < rIn, tapered to
// 0 at rOut. psi(r) = -int_0^r u_theta ds (u = (psi_y, -psi_x)).
struct rotation
{
    scalar omega, rIn, rOut;

    scalar uTheta(const scalar r) const
    {
        if (r <= rIn) return omega*r;
        if (r >= rOut) return 0;
        return omega*r*0.5*(1 + Foam::cos(constant::mathematical::pi*(r - rIn)/(rOut - rIn)));
    }

    scalar psi(const scalar r) const
    {
        if (r <= rIn) return -0.5*omega*r*r;
        const scalar r1 = min(r, rOut);
        // Simpson on [rIn, r1]
        const label n = 400;
        const scalar h = (r1 - rIn)/n;
        scalar s = uTheta(rIn) + uTheta(r1);
        for (label i = 1; i < n; ++i)
        {
            s += (i % 2 ? 4 : 2)*uTheta(rIn + i*h);
        }
        return -0.5*omega*rIn*rIn - s*h/3;
    }
};


int main(int argc, char *argv[])
{
    argList::addNote("WP 2.4: remap with prescribed motion");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "remapCheckDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    eulerianSolid solid(mesh);

    const word motion(dict.get<word>("motion"));
    const scalar endTime = dict.get<scalar>("endTime");
    const bool reverse = dict.getOrDefault<bool>("reverse", false);
    const scalar L = dict.getOrDefault<scalar>("L", 1);
    const scalar A = dict.getOrDefault<scalar>("A", 0);
    const vector Utrans = dict.getOrDefault<vector>("U", Zero);
    const point centre = dict.getOrDefault<point>("centre", Zero);
    rotation rot
    {
        dict.getOrDefault<scalar>("omega", 0),
        dict.getOrDefault<scalar>("rIn", 0),
        dict.getOrDefault<scalar>("rOut", 0)
    };
    const scalar pi = constant::mathematical::pi;

    // Cell size h: from the volume and the number of geometric directions
    const label nD = mesh.nGeometricD();
    scalar depth = 1;
    if (nD == 2)
    {
        const boundBox& bb = mesh.bounds();
        depth = bb.span().z();
    }
    const scalar h =
        Foam::pow(gMin(mesh.V().field())/depth, 1.0/scalar(nD));

    // Stream function for the 2D motions (psi at a point)
    auto psi = [&](const point& p) -> scalar
    {
        if (motion == "rotate")
        {
            return rot.psi(mag(vector(p.x() - centre.x(), p.y() - centre.y(), 0)));
        }
        // swirl
        return A*sqr(Foam::sin(pi*p.x()/L))*sqr(Foam::sin(pi*p.y()/L));
    };

    // Flux through one face (global face index), positive along its Sf
    const pointField& pts = mesh.points();
    auto faceFlux = [&](const label facei, const vector& Sf, const point& Cf)
        -> scalar
    {
        if (motion == "translate")
        {
            return Utrans & Sf;
        }
        if (motion == "compress")
        {
            return A*Foam::sin(2*pi*Cf.x()/L)*Sf.x();
        }
        // 2D stream function: edge a -> b of the face in the lowest z plane
        const face& f = mesh.faces()[facei];
        scalar zMin = GREAT;
        for (const label pti : f) zMin = min(zMin, pts[pti].z());
        DynamicList<point> ab;
        scalar zMax = -GREAT;
        for (const label pti : f)
        {
            zMax = max(zMax, pts[pti].z());
            if (mag(pts[pti].z() - zMin) < 1e-9*(1 + mag(zMin)))
            {
                ab.append(pts[pti]);
            }
        }
        if (ab.size() != 2) return 0;   // a face normal to z
        const vector e(ab[1] - ab[0]);
        const scalar s = Sf & vector(e.y(), -e.x(), 0);
        const scalar dpsi = (zMax - zMin)*(psi(ab[1]) - psi(ab[0]));
        return s > 0 ? dpsi : -dpsi;
    };

    // Max material speed for the time step
    scalar uMax = 0;
    if (motion == "translate") uMax = cmptMax(cmptMag(Utrans));
    else if (motion == "rotate") uMax = rot.omega*rot.rIn;
    else if (motion == "swirl") uMax = A*pi/L;
    else if (motion == "compress") uMax = A;
    else
    {
        FatalIOErrorInFunction(dict) << "Unknown motion " << motion
            << exit(FatalIOError);
    }

    label nSteps = dict.getOrDefault<label>("nSteps", 0);
    if (nSteps <= 0)
    {
        const scalar co = dict.get<scalar>("courant");
        nSteps = label(std::ceil(endTime*uMax/(co*h) - 1e-9));
    }
    const scalar dt = endTime/nSteps;

    // Fill phiL once (the motions are steady); the sign flips if reverse
    surfaceScalarField& phi = solid.phiLRef();
    {
        const vectorField& Sf = mesh.Sf().primitiveField();
        const vectorField& Cf = mesh.Cf().primitiveField();
        scalarField& phiI = phi.primitiveFieldRef();
        forAll(phiI, facei)
        {
            phiI[facei] = faceFlux(facei, Sf[facei], Cf[facei]);
        }
        auto& phiB = phi.boundaryFieldRef();
        forAll(phiB, patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!p.coupled())
            {
                phiB[patchi] = 0;   // walls: no flux
                continue;
            }
            const label start = p.patch().start();
            forAll(phiB[patchi], i)
            {
                phiB[patchi][i] = faceFlux(start + i, p.Sf()[i], p.Cf()[i]);
            }
        }
    }

    const scalarField& V = mesh.V().field();

    // Sum of dt*phi out of each cell (internal and coupled faces)
    auto netOut = [&](scalarField& out, scalarField& outflow)
    {
        out = 0;
        outflow = 0;
        const scalarField& phiI = phi.primitiveField();
        forAll(mesh.owner(), facei)
        {
            const label P = mesh.owner()[facei];
            const label N = mesh.neighbour()[facei];
            out[P] += dt*phiI[facei];
            out[N] -= dt*phiI[facei];
            outflow[phiI[facei] > 0 ? P : N] += dt*mag(phiI[facei]);
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& p = mesh.boundary()[patchi];
            if (!p.coupled()) continue;
            const scalarField& pp = phi.boundaryField()[patchi];
            forAll(p.faceCells(), i)
            {
                out[p.faceCells()[i]] += dt*pp[i];
                if (pp[i] > 0) outflow[p.faceCells()[i]] += dt*pp[i];
            }
        }
    };

    // Conserved totals: mass and sum(m q) of every component (+ sum m|q|)
    auto totals = [&](List<scalar>& tot, List<scalar>& scale)
    {
        const scalarField m(solid.rho().primitiveField()*V);
        tot.setSize(0);
        scale.setSize(0);
        tot.append(gSum(m));
        scale.append(gSum(m));
        auto add = [&](const scalarField& q)
        {
            tot.append(gSum(m*q));
            scale.append(gSum(m*mag(q)) + VSMALL);
        };
        for (direction c = 0; c < 3; ++c) add(solid.U().primitiveField().component(c));
        for (direction c = 0; c < 6; ++c) add(solid.sigma().primitiveField().component(c));
        for (direction c = 0; c < 6; ++c) add(solid.epsP().primitiveField().component(c));
        add(solid.epsPEq().primitiveField());
    };

    // Ranges of every component, initial and over the run
    struct rangeRec { word name; label cmpt; scalar min0, max0, minR, maxR; };
    DynamicList<rangeRec> ranges;
    auto track = [&](bool init)
    {
        label k = 0;
        auto rec = [&](const word& n, label c, const scalarField& q)
        {
            const scalar lo = gMin(q), hi = gMax(q);
            if (init)
            {
                ranges.append(rangeRec{n, c, lo, hi, lo, hi});
            }
            else
            {
                ranges[k].minR = min(ranges[k].minR, lo);
                ranges[k].maxR = max(ranges[k].maxR, hi);
            }
            ++k;
        };
        rec("rho", 0, solid.rho().primitiveField());
        for (direction c = 0; c < 3; ++c) rec("U", c, solid.U().primitiveField().component(c));
        for (direction c = 0; c < 6; ++c) rec("sigma", c, solid.sigma().primitiveField().component(c));
        for (direction c = 0; c < 6; ++c) rec("epsP", c, solid.epsP().primitiveField().component(c));
        rec("epsPEq", 0, solid.epsPEq().primitiveField());
    };

    List<scalar> tot0, scale0;
    totals(tot0, scale0);
    track(true);

    scalarField out(mesh.nCells()), outflow(mesh.nCells());
    scalar maxCo = 0;
    scalar maxCons = 0;

    for (label step = 0; step < nSteps; ++step)
    {
        if (reverse && step > 0)
        {
            phi *= -1.0;
        }
        netOut(out, outflow);
        maxCo = max(maxCo, gMax(outflow/V));

        // Lagrangian density: same mass in the moved volume V + out
        volScalarField& rho = solid.rhoRef();
        rho.primitiveFieldRef() *= V/(V + out);

        solid.remap(dt);

        List<scalar> tot, scale;
        totals(tot, scale);
        forAll(tot, i)
        {
            maxCons = max(maxCons, mag(tot[i] - tot0[i])/scale0[i]);
        }
        track(false);
    }

    const label oldPrecision = Info().precision(17);
    Info<< "steps " << nSteps << " dt " << dt << nl
        << "cellSize " << h << nl
        << "maxCourantSum " << maxCo << nl
        << "conservation " << maxCons << nl;
    for (const rangeRec& r : ranges)
    {
        Info<< "range " << r.name << ' ' << r.cmpt << ' ' << r.min0 << ' '
            << r.max0 << ' ' << r.minR << ' ' << r.maxR << nl;
    }
    Info().precision(oldPrecision);

    runTime.setTime(endTime, 1);
    solid.rhoRef().write();
    solid.URef().write();
    solid.sigmaRef().write();
    solid.epsPRef().write();
    solid.epsPEqRef().write();

    Info<< "End" << endl;
    return 0;
}


// ************************************************************************* //
