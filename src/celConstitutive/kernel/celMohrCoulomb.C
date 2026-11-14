#include "celMohrCoulomb.H"
#include "celEigen.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cel
{
namespace constitutive
{

namespace
{

constexpr double pi = 3.14159265358979323846;

//- A plane in sorted principal space: yield normal nf, flow normal ng
struct Plane
{
    double nf[3];
    double ng[3];
};

double dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//- Principal elastic matrix times v: (K - 2G/3) tr(v) + 2G v
void principalD(const Elastic& el, const double v[3], double out[3])
{
    const double lambda = el.K - 2.0*el.G/3.0;
    const double tr = v[0] + v[1] + v[2];
    for (int i = 0; i < 3; ++i)
    {
        out[i] = lambda*tr + 2.0*el.G*v[i];
    }
}

//- Return onto the active planes (one or two): solve
//  A dg = b, A_ij = nf_i . Dp ng_j, b_i = nf_i . w - k; s = w - Dp sum dg_j ng_j
void returnTo
(
    const Elastic& el,
    const Plane* planes,
    int n,
    const double w[3],
    double k,
    double s[3],
    double dg[2]
)
{
    double Dng[2][3];
    for (int j = 0; j < n; ++j)
    {
        principalD(el, planes[j].ng, Dng[j]);
    }
    double b[2];
    for (int i = 0; i < n; ++i)
    {
        b[i] = dot3(planes[i].nf, w) - k;
    }

    if (n == 1)
    {
        dg[0] = b[0]/dot3(planes[0].nf, Dng[0]);
        dg[1] = 0;
    }
    else
    {
        const double a00 = dot3(planes[0].nf, Dng[0]);
        const double a01 = dot3(planes[0].nf, Dng[1]);
        const double a10 = dot3(planes[1].nf, Dng[0]);
        const double a11 = dot3(planes[1].nf, Dng[1]);
        const double det = a00*a11 - a01*a10;
        dg[0] = (b[0]*a11 - a01*b[1])/det;
        dg[1] = (a00*b[1] - a10*b[0])/det;
    }

    for (int i = 0; i < 3; ++i)
    {
        s[i] = w[i];
        for (int j = 0; j < n; ++j)
        {
            s[i] -= dg[j]*Dng[j][i];
        }
    }
}

bool sorted(const double s[3], double tol)
{
    return s[0] >= s[1] - tol && s[1] >= s[2] - tol;
}

} // End anonymous namespace


MohrCoulomb::MohrCoulomb(const Parameters& p)
:
    el_(parameter(p, "E"), parameter(p, "nu")),
    c_(parameter(p, "c")),
    sinPhi_(std::sin(parameter(p, "phi")*pi/180.0)),
    cosPhi_(std::cos(parameter(p, "phi")*pi/180.0)),
    sinPsi_(std::sin(parameter(p, "psi")*pi/180.0))
{
    const double phi = parameter(p, "phi");
    const double psi = parameter(p, "psi");
    if (c_ < 0)
    {
        throw std::invalid_argument("mohrCoulomb: c must not be negative");
    }
    if (phi < 0 || phi >= 90)
    {
        throw std::invalid_argument("mohrCoulomb: phi must be in [0, 90) degrees");
    }
    if (psi < 0 || psi > phi)
    {
        throw std::invalid_argument("mohrCoulomb: psi must be in [0, phi]");
    }
    if (phi == 0 && c_ == 0)
    {
        throw std::invalid_argument("mohrCoulomb: phi = 0 needs c > 0");
    }
}


void MohrCoulomb::update
(
    const Voigt deps,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    updateReturn(deps, sigOld, hOld, sigNew, hNew);
}


MohrCoulomb::Return MohrCoulomb::updateReturn
(
    const Voigt deps,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    Voigt trial;
    el_.addStress(sigOld, deps, trial);

    hNew = hOld;
    hNew.apex = false;

    double A[3][3], w[3], V[3][3];
    voigtToMatrix(trial, A);
    symmetricEigen(A, w, V);

    const double sphi = sinPhi_;
    const double spsi = sinPsi_;
    const double k = 2.0*c_*cosPhi_;
    const double scale =
        std::max({std::abs(w[0]), std::abs(w[1]), std::abs(w[2]), c_, 1.0});
    const double fTrial = (w[0] - w[2]) + (w[0] + w[2])*sphi - k;

    if (fTrial <= yieldTol*scale)
    {
        std::copy(trial, trial + 6, sigNew);
        return Return::elastic;
    }

    const Plane main
    {
        {1 + sphi, 0, -1 + sphi},
        {1 + spsi, 0, -1 + spsi}
    };
    const Plane edge12
    {
        {1 + sphi, -1 + sphi, 0},
        {1 + spsi, -1 + spsi, 0}
    };
    const Plane edge23
    {
        {0, 1 + sphi, -1 + sphi},
        {0, 1 + spsi, -1 + spsi}
    };

    const double tol = 1e-12*scale;
    Return ret = Return::plane;
    bool apex = false;
    double s[3], dg[2];

    returnTo(el_, &main, 1, w, k, s, dg);
    if (!sorted(s, tol))
    {
        Plane planes[2] = {main, edge12};
        ret = Return::compressionEdge;
        if (s[1] > s[0])
        {
            // s1 dropped below s2: triaxial compression corner s1 = s2,
            // active planes (s1, s3) and (s2, s3)
            planes[1] = edge23;
        }
        else
        {
            // s3 rose above s2: triaxial extension corner s2 = s3
            ret = Return::extensionEdge;
        }
        returnTo(el_, planes, 2, w, k, s, dg);
        const bool ok = dg[0] >= -tol && dg[1] >= -tol && sorted(s, tol);
        if (!ok)
        {
            if (sphi == 0)
            {
                throw std::runtime_error("mohrCoulomb: Tresca corner return failed");
            }
            apex = true;
        }
    }
    if (!apex && sphi > 0 && 0.5*(s[0] + s[2]) > c_*cosPhi_/sphi + tol)
    {
        apex = true;
    }
    if (apex)
    {
        s[0] = s[1] = s[2] = c_*cosPhi_/sphi;
        ret = Return::apex;
    }

    spectralToVoigt(s, V, sigNew);

    // History from the plastic corrector: deps_p = C (trial - sigNew)
    plasticCorrector(el_, trial, sigNew, hOld, apex, hNew);

    return ret;
}


double MohrCoulomb::yield(const Voigt sigma) const
{
    double A[3][3], w[3], V[3][3];
    voigtToMatrix(sigma, A);
    symmetricEigen(A, w, V);
    return (w[0] - w[2]) + (w[0] + w[2])*sinPhi_ - 2.0*c_*cosPhi_;
}

} // End namespace constitutive
} // End namespace cel
