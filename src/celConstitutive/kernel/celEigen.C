#include "celEigen.H"

#include <cmath>
#include <utility>

namespace cel
{
namespace constitutive
{

void symmetricEigen(const double A[3][3], double w[3], double V[3][3])
{
    double a[3][3];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            a[i][j] = (j >= i) ? A[i][j] : A[j][i];
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Cyclic Jacobi. Convergence is quadratic; 50 sweeps is far more than the
    // 5 or 6 a 3x3 matrix needs, and the loop stops once the off-diagonal
    // part is zero or negligible against the diagonal.
    for (int sweep = 0; sweep < 50; ++sweep)
    {
        const double off =
            std::abs(a[0][1]) + std::abs(a[0][2]) + std::abs(a[1][2]);
        const double diag =
            std::abs(a[0][0]) + std::abs(a[1][1]) + std::abs(a[2][2]);
        if (off <= 1e-22*diag || off < 1e-300)
        {
            break;
        }

        for (int p = 0; p < 2; ++p)
        {
            for (int q = p + 1; q < 3; ++q)
            {
                const double apq = a[p][q];
                if (apq == 0)
                {
                    continue;
                }

                // Rotation angle (Numerical Recipes, jacobi)
                const double theta = (a[q][q] - a[p][p])/(2.0*apq);
                double t = 1.0/(std::abs(theta) + std::sqrt(theta*theta + 1.0));
                if (theta < 0)
                {
                    t = -t;
                }
                if (!std::isfinite(t))
                {
                    t = 0.5/theta;  // theta*theta overflowed: t ~ 1/(2 theta)
                }
                const double c = 1.0/std::sqrt(t*t + 1.0);
                const double s = t*c;

                // a <- J^T a J with J the rotation in the (p, q) plane
                for (int k = 0; k < 3; ++k)
                {
                    const double akp = a[k][p];
                    const double akq = a[k][q];
                    a[k][p] = c*akp - s*akq;
                    a[k][q] = s*akp + c*akq;
                }
                for (int k = 0; k < 3; ++k)
                {
                    const double apk = a[p][k];
                    const double aqk = a[q][k];
                    a[p][k] = c*apk - s*aqk;
                    a[q][k] = s*apk + c*aqk;
                }
                a[p][q] = a[q][p] = 0;

                for (int k = 0; k < 3; ++k)
                {
                    const double vkp = V[k][p];
                    const double vkq = V[k][q];
                    V[k][p] = c*vkp - s*vkq;
                    V[k][q] = s*vkp + c*vkq;
                }
            }
        }
    }

    for (int i = 0; i < 3; ++i)
    {
        w[i] = a[i][i];
    }

    // Sort descending, permuting the eigenvector columns with the values
    for (int i = 0; i < 2; ++i)
    {
        int m = i;
        for (int j = i + 1; j < 3; ++j)
        {
            if (w[j] > w[m])
            {
                m = j;
            }
        }
        if (m != i)
        {
            std::swap(w[i], w[m]);
            for (int k = 0; k < 3; ++k)
            {
                std::swap(V[k][i], V[k][m]);
            }
        }
    }
}


void voigtToMatrix(const double s[6], double A[3][3])
{
    A[0][0] = s[0];
    A[1][1] = s[1];
    A[2][2] = s[2];
    A[0][1] = A[1][0] = s[3];
    A[1][2] = A[2][1] = s[4];
    A[2][0] = A[0][2] = s[5];
}


void spectralToVoigt(const double w[3], const double V[3][3], double s[6])
{
    auto entry = [&](int i, int j)
    {
        return V[i][0]*w[0]*V[j][0] + V[i][1]*w[1]*V[j][1] + V[i][2]*w[2]*V[j][2];
    };
    s[0] = entry(0, 0);
    s[1] = entry(1, 1);
    s[2] = entry(2, 2);
    s[3] = entry(0, 1);
    s[4] = entry(1, 2);
    s[5] = entry(2, 0);
}

} // End namespace constitutive
} // End namespace cel
