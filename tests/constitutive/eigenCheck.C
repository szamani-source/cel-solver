/*---------------------------------------------------------------------------*\
Application
    eigenCheck

Description
    WP 2.5b unit test of the kernel's symmetric 3x3 eigen-solver
    (symmetricEigen, kernel/celEigen.H), which Mohr-Coulomb uses for its
    principal stresses.

    Matrices are built as R diag(lambda) R^T with a known random rotation R,
    so the exact eigenvalues are known. Cases: distinct, hydrostatic (three
    equal), triaxial (two equal, both orders), nearly repeated (gap 1e-9),
    zero, pure shear, diagonal, wide magnitude range, and 10000 random ones.
    Each case requires, relative to max |lambda|:
    - eigenvalues sorted descending and within 1e-13 of the exact ones,
    - V orthonormal to 1e-13,
    - V diag(w) V^T = A to 1e-13 (unique even when eigenvectors are not),
    - A v_k = w_k v_k to 1e-13.

    Build and run: tests/test_constitutive_cpp.sh
\*---------------------------------------------------------------------------*/

#include "celEigen.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

using namespace cel::constitutive;

namespace
{

std::mt19937_64 rng(25);

//- Random rotation from a random unit quaternion
void randomRotation(double R[3][3])
{
    std::normal_distribution<double> n(0, 1);
    double q[4];
    double norm = 0;
    for (double& v : q)
    {
        v = n(rng);
        norm += v*v;
    }
    norm = std::sqrt(norm);
    for (double& v : q) v /= norm;
    const double w = q[0], x = q[1], y = q[2], z = q[3];
    R[0][0] = 1 - 2*(y*y + z*z); R[0][1] = 2*(x*y - z*w); R[0][2] = 2*(x*z + y*w);
    R[1][0] = 2*(x*y + z*w); R[1][1] = 1 - 2*(x*x + z*z); R[1][2] = 2*(y*z - x*w);
    R[2][0] = 2*(x*z - y*w); R[2][1] = 2*(y*z + x*w); R[2][2] = 1 - 2*(x*x + y*y);
}


//- One case: A = R diag(lam) R^T (R = identity if rotate is false).
//  Returns the worst relative error found.
double check(const double lam[3], bool rotate)
{
    double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    if (rotate)
    {
        randomRotation(R);
    }
    double A[3][3];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            A[i][j] = 0;
            for (int k = 0; k < 3; ++k)
            {
                A[i][j] += R[i][k]*lam[k]*R[j][k];
            }
        }
    }

    double w[3], V[3][3];
    symmetricEigen(A, w, V);

    double exact[3] = {lam[0], lam[1], lam[2]};
    std::sort(exact, exact + 3, [](double a, double b) { return a > b; });
    const double scale =
        std::max({std::abs(lam[0]), std::abs(lam[1]), std::abs(lam[2]), 1e-300});

    double err = 0;
    if (!(w[0] >= w[1] && w[1] >= w[2]))
    {
        err = 1;  // not sorted
    }
    for (int k = 0; k < 3; ++k)
    {
        err = std::max(err, std::abs(w[k] - exact[k])/scale);
    }
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            double vv = 0, rec = 0, av = 0;
            for (int i = 0; i < 3; ++i)
            {
                vv += V[i][a]*V[i][b];
                rec += V[a][i]*w[i]*V[b][i];
                av += A[a][i]*V[i][b];
            }
            err = std::max(err, std::abs(vv - (a == b ? 1.0 : 0.0)));
            err = std::max(err, std::abs(rec - A[a][b])/scale);
            err = std::max(err, std::abs(av - w[b]*V[a][b])/scale);
        }
    }
    return err;
}

} // End anonymous namespace


int main()
{
    const double tol = 1e-13;
    int failures = 0;

    struct Named
    {
        const char* name;
        double lam[3];
        bool rotate;
    };
    const Named cases[] =
    {
        {"distinct", {3e5, -1e5, -2e5}, true},
        {"hydrostatic (three equal)", {-1e5, -1e5, -1e5}, true},
        {"triaxial compression (s1 = s2)", {-1e5, -1e5, -3e5}, true},
        {"triaxial extension (s2 = s3)", {-1e5, -3e5, -3e5}, true},
        {"nearly repeated (gap 1e-9)", {-1e5, -1e5*(1 + 1e-9), -3e5}, true},
        {"zero matrix", {0, 0, 0}, true},
        {"pure shear", {1e5, 0, -1e5}, true},
        {"already diagonal, unsorted", {-2, 5, 1}, false},
        {"wide magnitude range", {1e9, 1e-3, -1e2}, true},
        {"tiny values", {1e-200, -3e-200, 2e-200}, true},
    };

    for (const Named& c : cases)
    {
        const double err = check(c.lam, c.rotate);
        const bool ok = err < tol;
        std::printf("%s  %-34s max rel error %.2e\n", ok ? "PASS" : "FAIL", c.name, err);
        failures += !ok;
    }

    // Random spectra, including random repeated pairs
    std::uniform_real_distribution<double> u(-1, 1);
    double worst = 0;
    for (int n = 0; n < 10000; ++n)
    {
        double lam[3] = {u(rng)*1e6, u(rng)*1e6, u(rng)*1e6};
        if (n % 4 == 1) lam[1] = lam[0];
        if (n % 4 == 2) lam[2] = lam[1] = lam[0];
        worst = std::max(worst, check(lam, true));
    }
    const bool ok = worst < tol;
    std::printf("%s  %-34s max rel error %.2e\n", ok ? "PASS" : "FAIL",
                "10000 random (with repeats)", worst);
    failures += !ok;

    std::printf("\n%s (tolerance %.0e)\n",
                failures ? "FAILED" : "All eigen checks passed", tol);
    return failures ? 1 : 0;
}
