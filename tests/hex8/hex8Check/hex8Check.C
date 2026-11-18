/*---------------------------------------------------------------------------*\
Application
    celHex8Check

Description
    WP 5.1 done-when checks of the explicit hex8 body (decision 0070):

    [1] Patch test (MacNeal and Harder 1985, seven distorted hexahedra in a
        unit cube): every node moves with a linear velocity field (strain
        1e-3, all components); after one step the stress is the same in
        every element (to 1e-12 of its size) and close to the exact
        2000 / 400 Pa (E 1e6, nu 0.25; within 1e-3, the finite-strain
        difference), and the force on each inner node is zero (to 1e-10 of
        V sigma / L).
    [2] One-dimensional strain wave: steel bar of 50 elements, lateral
        motion held, one end driven at 1 m/s: behind the front the stress is
        -rho c v (within 2 %), the front at c t (within 2 elements).
    [3] Energy and hourglass: a free elastic cube (4 x 4 x 4 elements, one
        corner node pushed out of line so hourglass modes start), 2000
        steps: the total energy stays within 1 % of the
        start (kinetic, elastic strain energy from the stresses, and the
        hourglass energy), momentum kept to 1e-12.

    Prints each check with PASS / FAIL; exit code 1 on any failure.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IStringStream.H"
#include "hex8Body.H"

using namespace Foam;

namespace
{

dictionary dictOf(const string& s)
{
    IStringStream is(s);
    return dictionary(is);
}


//- Structured block of nx x ny x nz elements of size h
void block
(
    const label nx, const label ny, const label nz, const scalar h,
    pointField& p, List<hex8Body::hexElement>& e
)
{
    auto id = [&](label i, label j, label k)
    {
        return i + (nx + 1)*(j + (ny + 1)*k);
    };
    p.setSize((nx + 1)*(ny + 1)*(nz + 1));
    for (label k = 0; k <= nz; ++k)
        for (label j = 0; j <= ny; ++j)
            for (label i = 0; i <= nx; ++i)
                p[id(i, j, k)] = point(i*h, j*h, k*h);
    e.setSize(nx*ny*nz);
    label n = 0;
    for (label k = 0; k < nz; ++k)
        for (label j = 0; j < ny; ++j)
            for (label i = 0; i < nx; ++i)
            {
                hex8Body::hexElement& c = e[n++];
                c[0] = id(i, j, k);         c[1] = id(i + 1, j, k);
                c[2] = id(i + 1, j + 1, k); c[3] = id(i, j + 1, k);
                c[4] = id(i, j, k + 1);     c[5] = id(i + 1, j, k + 1);
                c[6] = id(i + 1, j + 1, k + 1); c[7] = id(i, j + 1, k + 1);
            }
}


bool report(const word& name, const bool ok, const string& detail)
{
    Info<< (ok ? "[PASS] " : "[FAIL] ") << name.c_str() << nl
        << "    " << detail.c_str() << nl;
    return ok;
}


bool patchTest()
{
    pointField p(16);
    p[0] = point(0, 0, 0); p[1] = point(1, 0, 0);
    p[2] = point(1, 1, 0); p[3] = point(0, 1, 0);
    p[4] = point(0, 0, 1); p[5] = point(1, 0, 1);
    p[6] = point(1, 1, 1); p[7] = point(0, 1, 1);
    p[8] = point(0.249, 0.342, 0.192);  p[9] = point(0.826, 0.288, 0.288);
    p[10] = point(0.850, 0.649, 0.263); p[11] = point(0.273, 0.750, 0.230);
    p[12] = point(0.320, 0.186, 0.643); p[13] = point(0.677, 0.305, 0.683);
    p[14] = point(0.788, 0.693, 0.644); p[15] = point(0.165, 0.745, 0.702);

    const label quads[7][8] =
    {
        {8, 9, 10, 11, 12, 13, 14, 15},
        {0, 1, 2, 3, 8, 9, 10, 11},
        {4, 5, 6, 7, 12, 13, 14, 15},
        {0, 1, 5, 4, 8, 9, 13, 12},
        {3, 2, 6, 7, 11, 10, 14, 15},
        {0, 3, 7, 4, 8, 11, 15, 12},
        {1, 2, 6, 5, 9, 10, 14, 13}
    };
    List<hex8Body::hexElement> e(7);
    forAll(e, i)
    {
        for (label k = 0; k < 8; ++k) e[i][k] = quads[i][k];
    }
    hex8Body body
    (
        "patch", p, e,
        dictOf("density 1; material { type linearElastic; E 1e6; nu 0.25; }")
    );

    const symmTensor eps(1e-3, 0.5e-3, 0.5e-3, 1e-3, 0.5e-3, 1e-3);
    vectorField v(p.size());
    forAll(p, i)
    {
        v[i] = eps & p[i];
        if (i < 8)
        {
            body.fix(i, vector::one, v[i]);
        }
    }
    body.setVelocity(v);
    body.step(1.0);

    const symmTensorField& s = body.sigma();
    scalar dev = 0;
    forAll(s, i)
    {
        dev = max(dev, cmptMax(cmptMag(s[i] - s[0])));
    }
    dev /= 2000;
    const symmTensor exact(2000, 400, 400, 2000, 400, 2000);
    const scalar err = cmptMax(cmptMag(s[0] - exact))/2000;
    const vectorField f(body.forces());
    scalar fin = 0;
    for (label i = 8; i < 16; ++i)
    {
        fin = max(fin, mag(f[i]));
    }
    fin /= 2000*0.1;     // V sigma / L of an element
    const bool ok = dev <= 1e-12 && err <= 1e-3 && fin <= 1e-10;
    return report
    (
        "Patch test (MacNeal-Harder, 7 distorted hexahedra)", ok,
        "stress spread " + Foam::name(dev) + " (<= 1e-12), against exact "
      + Foam::name(err) + " (<= 1e-3), inner node force "
      + Foam::name(fin) + " (<= 1e-10)"
    );
}


bool waveTest()
{
    const label n = 50;
    const scalar h = 0.01, v0 = 1, rho = 7800, E = 2e11, nu = 0.3;
    pointField p;
    List<hex8Body::hexElement> e;
    block(n, 1, 1, h, p, e);
    hex8Body body
    (
        "bar", p, e,
        dictOf("density 7800; material { type linearElastic; E 2e11; nu 0.3; }")
    );
    forAll(p, i)
    {
        if (p[i].x() < 0.5*h)
        {
            body.fix(i, vector::one, vector(v0, 0, 0));
        }
        else
        {
            body.fix(i, vector(0, 1, 1), Zero);
        }
    }
    const scalar M = E*(1 - nu)/((1 + nu)*(1 - 2*nu));
    const scalar c = Foam::sqrt(M/rho);
    const scalar dt = 0.5*body.stableDeltaT();
    const scalar tEnd = 0.5*n*h/c;
    const label nSteps = label(tEnd/dt);
    for (label k = 0; k < nSteps; ++k)
    {
        body.step(dt);
    }
    const scalar t = nSteps*dt;
    const scalar sRef = -rho*c*v0;
    scalar mean = 0;
    const label nBehind = n/4;
    for (label i = 2; i < nBehind; ++i)
    {
        mean += body.sigma()[i].xx();
    }
    mean /= (nBehind - 2);
    // Front: first element (from the driven end) below half the stress
    label front = n;
    for (label i = 0; i < n; ++i)
    {
        if (body.sigma()[i].xx() > 0.5*sRef)
        {
            front = i;
            break;
        }
    }
    const scalar xf = front*h, xRef = c*t;
    const bool ok = mag(mean/sRef - 1) <= 0.02 && mag(xf - xRef) <= 2*h;
    return report
    (
        "One-dimensional strain wave in a steel bar", ok,
        "stress behind the front " + Foam::name(mean/1e6) + " MPa against "
      + Foam::name(sRef/1e6) + " MPa (2 %), front at " + Foam::name(xf)
      + " m against c t = " + Foam::name(xRef) + " m (2 elements)"
    );
}


bool energyTest()
{
    pointField p;
    List<hex8Body::hexElement> e;
    block(4, 4, 4, 0.025, p, e);
    p[0] += vector(-0.004, -0.003, -0.002);     // distort a corner
    hex8Body body
    (
        "cube", p, e,
        dictOf("density 2700; material { type linearElastic; E 7e10; nu 0.33; }")
    );
    vectorField v(p.size());
    forAll(p, i)
    {
        // Alternating (hourglass-like) pattern plus a rotation
        const label s = (label(p[i].x()/0.025 + 0.5) + label(p[i].y()/0.025 + 0.5)
          + label(p[i].z()/0.025 + 0.5)) % 2;
        v[i] = vector(s ? 0.1 : -0.1, 0.05, 0) + (vector(0, 0, 2) ^ p[i]);
    }
    body.setVelocity(v);
    const scalar KE0 = body.kineticEnergy();
    const vector P0(sum(body.mass()*body.velocity()));
    const scalar dt = 0.5*body.stableDeltaT();
    const scalar E = 7e10, nu = 0.33;
    auto stored = [&]()
    {
        // Elastic strain energy from the stresses, plus the hourglass part
        const scalarField V(body.volumes());
        scalar U = 0;
        forAll(V, ei)
        {
            const symmTensor& s = body.sigma()[ei];
            U += 0.5*V[ei]*((1 + nu)*(s && s) - nu*sqr(tr(s)))/E;
        }
        return U + body.hourglassEnergy();
    };
    // Central difference keeps 0.5 v(n-1/2) M v(n+1/2) + U(x_n) (the
    // kinetic energy from the two half-step velocities around x_n)
    scalar worst = 0, hgMax = 0;
    body.step(dt);
    scalar U = stored();
    vectorField vPrev(body.velocity());
    scalar E0 = -1;
    for (label k = 1; k < 2000; ++k)
    {
        body.step(dt);
        const scalar KE = 0.5*sum(body.mass()*(vPrev & body.velocity()));
        const scalar Etot = KE + U;
        if (E0 < 0) E0 = Etot;
        worst = max(worst, mag(Etot - E0)/E0);
        hgMax = max(hgMax, body.hourglassEnergy()/E0);
        U = stored();
        vPrev = body.velocity();
    }
    const vector P(sum(body.mass()*body.velocity()));
    const scalar dP = mag(P - P0)/max(mag(P0), VSMALL);
    const bool ok = worst <= 0.01 && dP <= 1e-12;
    return report
    (
        "Energy and hourglass, free elastic cube, 2000 steps", ok,
        "kinetic + strain + hourglass energy: largest change "
      + Foam::name(worst) + " (<= 0.01), hourglass energy up to "
      + Foam::name(hgMax) + " of the start, momentum "
      + Foam::name(dP) + " (<= 1e-12), kinetic energy at the end "
      + Foam::name(body.kineticEnergy()/KE0) + " of the start"
    );
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noBanner();
    argList::noCheckProcessorDirectories();
    argList args(argc, argv, false, false, false);

    bool ok = patchTest();
    ok = waveTest() && ok;
    ok = energyTest() && ok;
    Info<< nl << "WP 5.1 checks: " << (ok ? "ALL PASS" : "FAILED") << nl;
    return ok ? 0 : 1;
}
