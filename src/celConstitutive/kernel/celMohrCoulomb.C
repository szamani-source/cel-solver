#include "celMohrCoulomb.H"
#include "celEigen.H"

#include <algorithm>
#include <cmath>
#include <limits>
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


//- Solve the n x n system A x = b (n <= 3) by Gaussian elimination with
//  partial pivoting. Returns false if A is singular to the test used by
//  models.py: |det A| <= 1e-12 prod |A_ii|.
bool solveSmall(int n, double A[3][3], double b[3], double x[3])
{
    double diag = 1;
    for (int i = 0; i < n; ++i)
    {
        diag *= std::abs(A[i][i]);
    }

    double det = 1;
    for (int k = 0; k < n; ++k)
    {
        int piv = k;
        for (int i = k + 1; i < n; ++i)
        {
            if (std::abs(A[i][k]) > std::abs(A[piv][k]))
            {
                piv = i;
            }
        }
        if (piv != k)
        {
            for (int j = 0; j < n; ++j)
            {
                std::swap(A[k][j], A[piv][j]);
            }
            std::swap(b[k], b[piv]);
            det = -det;
        }
        det *= A[k][k];
        if (A[k][k] == 0)
        {
            return false;
        }
        for (int i = k + 1; i < n; ++i)
        {
            const double f = A[i][k]/A[k][k];
            for (int j = k; j < n; ++j)
            {
                A[i][j] -= f*A[k][j];
            }
            b[i] -= f*b[k];
        }
    }
    if (std::abs(det) <= 1e-12*diag)
    {
        return false;
    }
    for (int i = n - 1; i >= 0; --i)
    {
        double r = b[i];
        for (int j = i + 1; j < n; ++j)
        {
            r -= A[i][j]*x[j];
        }
        x[i] = r/A[i][i];
    }
    return true;
}

} // End anonymous namespace


MohrCoulomb::MohrCoulomb(const Parameters& p)
:
    el_(parameter(p, "E"), parameter(p, "nu")),
    c_(parameter(p, "c")),
    sinPhi_(std::sin(parameter(p, "phi")*pi/180.0)),
    cosPhi_(std::cos(parameter(p, "phi")*pi/180.0)),
    sinPsi_(std::sin(parameter(p, "psi")*pi/180.0)),
    tension_(p.count("sigmaT") > 0),
    sigmaT_(parameter(p, "sigmaT", 0))
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
    if (tension_)
    {
        const double apex =
            phi == 0
          ? std::numeric_limits<double>::infinity()
          : c_*cosPhi_/sinPhi_;
        if (!(sigmaT_ >= 0 && sigmaT_ <= apex*(1 + 1e-12)))
        {
            throw std::invalid_argument
            (
                "mohrCoulomb: sigmaT must be in [0, c cot(phi)] = [0, "
              + std::to_string(apex) + "] Pa"
            );
        }
    }

    law_ = FrictionLaw::New(p, phi);
    nl_ = Nonlocal::New(p, phi);
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


void MohrCoulomb::updateRate
(
    const Voigt deps,
    double dt,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    if (!law_)
    {
        update(deps, sigOld, hOld, sigNew, hNew);
        return;
    }
    Voigt trial;
    el_.addStress(sigOld, deps, trial);
    muIReturn
    (
        *law_, el_, trial, dt,
        [&](double mu, Voigt s, History& h)
        {
            if (mu == law_->muS() && !law_->dilatant())
            {
                update(deps, sigOld, hOld, s, h);
                return;
            }
            // phi_m in degrees and back, as mu_i.py, so both sides round alike
            const double phi = std::atan(mu)*180.0/pi*pi/180.0;
            Friction fr{std::sin(phi), std::cos(phi), sigmaT_, sinPsi_};
            if (law_->dilatant())
            {
                const double t = law_->tanPsi(hOld, mu);
                fr.sinPsi = t/std::sqrt(1 + t*t);
                // Contraction is bounded so that the plane return's
                // a_f.D.a_g = 4G + 4 sin(phi) sin(psi) (lambda + G) stays
                // at least 2G (positive plastic multiplier)
                const double lambda = el_.K - 2.0*el_.G/3.0;
                fr.sinPsi =
                    std::max
                    (
                        fr.sinPsi,
                        -0.5*el_.G/((lambda + el_.G)*fr.sinPhi)
                    );
            }
            if (tension_)
            {
                fr.sigmaT = std::min(sigmaT_, c_*fr.cosPhi/fr.sinPhi);
            }
            returnAt(fr, deps, sigOld, hOld, s, h);
        },
        sigNew, hNew
    );
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
    return returnAt
    (
        Friction{sinPhi_, cosPhi_, sigmaT_, sinPsi_}, deps, sigOld, hOld, sigNew, hNew
    );
}


MohrCoulomb::Return MohrCoulomb::returnAt
(
    const Friction& fr,
    const Voigt deps,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    const double sigT = fr.sigmaT;
    const double cphi = fr.cosPhi;

    Voigt trial;
    el_.addStress(sigOld, deps, trial);

    hNew = hOld;
    hNew.apex = false;

    double A[3][3], w[3], V[3][3];
    voigtToMatrix(trial, A);
    symmetricEigen(A, w, V);

    const double sphi = fr.sinPhi;
    const double spsi = fr.sinPsi;
    const double k = 2.0*c_*cphi;
    const double scale =
        std::max({std::abs(w[0]), std::abs(w[1]), std::abs(w[2]), c_, 1.0});
    const double fTrial = (w[0] - w[2]) + (w[0] + w[2])*sphi - k;
    const double tTrial =
        tension_ ? w[0] - sigT : -std::numeric_limits<double>::infinity();

    if (fTrial <= yieldTol*scale && tTrial <= yieldTol*scale)
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
    bool mcDone = false;
    double s[3], dg[2];

    if (fTrial > yieldTol*scale)
    {
        mcDone = true;
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
                    throw std::runtime_error
                    (
                        "mohrCoulomb: Tresca corner return failed"
                    );
                }
                apex = true;
            }
        }
        if (!apex && sphi > 0 && 0.5*(s[0] + s[2]) > c_*cphi/sphi + tol)
        {
            apex = true;
        }
        if (apex)
        {
            s[0] = s[1] = s[2] = c_*cphi/sphi;
            ret = Return::apex;
        }
    }

    if (tension_ && (!mcDone || s[0] > sigT + tol))
    {
        ret = tensionReturn(fr, w, scale, s);
        apex = (ret == Return::tensionApex);
    }

    spectralToVoigt(s, V, sigNew);

    // History from the plastic corrector: deps_p = C (trial - sigNew)
    plasticCorrector(el_, trial, sigNew, hOld, apex, hNew);

    return ret;
}


MohrCoulomb::Return MohrCoulomb::tensionReturn
(
    const Friction& fr,
    const double w[3],
    double scale,
    double s[3]
) const
{
    const double sigT = fr.sigmaT;
    const double sphi = fr.sinPhi;
    const double spsi = fr.sinPsi;
    const double k = 2.0*c_*fr.cosPhi;

    // Planes in the order of _MC_T_SETS in models.py: main, s1 = s2 edge
    // (edge_23), s2 = s3 edge (edge_12), then the Rankine planes T1, T2, T3
    const Plane planes[6] =
    {
        {{1 + sphi, 0, -1 + sphi}, {1 + spsi, 0, -1 + spsi}},
        {{0, 1 + sphi, -1 + sphi}, {0, 1 + spsi, -1 + spsi}},
        {{1 + sphi, -1 + sphi, 0}, {1 + spsi, -1 + spsi, 0}},
        {{1, 0, 0}, {1, 0, 0}},
        {{0, 1, 0}, {0, 1, 0}},
        {{0, 0, 1}, {0, 0, 1}}
    };
    const double rhs[6] = {k, k, k, sigT, sigT, sigT};

    double Dng[6][3];
    for (int j = 0; j < 6; ++j)
    {
        principalD(el_, planes[j].ng, Dng[j]);
    }

    const double M = el_.K + 4.0*el_.G/3.0;
    const double tol = 1e-10*std::max(scale, sigT);

    // Active sets: one, two, then three planes, in lexicographic order
    int sets[41][3];
    int nSet[41];
    int m = 0;
    for (int i = 0; i < 6; ++i)
    {
        sets[m][0] = i; nSet[m++] = 1;
    }
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            sets[m][0] = i; sets[m][1] = j; nSet[m++] = 2;
        }
    }
    for (int i = 0; i < 6; ++i)
    {
        for (int j = i + 1; j < 6; ++j)
        {
            for (int l = j + 1; l < 6; ++l)
            {
                sets[m][0] = i; sets[m][1] = j; sets[m][2] = l; nSet[m++] = 3;
            }
        }
    }

    for (int a = 0; a < m; ++a)
    {
        const int n = nSet[a];
        const int* act = sets[a];

        double A[3][3], b[3], dg[3];
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                A[i][j] = dot3(planes[act[i]].nf, Dng[act[j]]);
            }
            b[i] = dot3(planes[act[i]].nf, w) - rhs[act[i]];
        }
        if (!solveSmall(n, A, b, dg))
        {
            continue;
        }

        bool ok = true;
        for (int i = 0; i < n; ++i)
        {
            ok = ok && dg[i]*M >= -tol;
        }
        if (!ok)
        {
            continue;
        }

        double t[3];
        for (int i = 0; i < 3; ++i)
        {
            t[i] = w[i];
            for (int j = 0; j < n; ++j)
            {
                t[i] -= dg[j]*Dng[act[j]][i];
            }
        }
        if (!sorted(t, tol))
        {
            continue;
        }
        const double fMC = (t[0] - t[2]) + (t[0] + t[2])*sphi - k;
        if (fMC > tol || t[0] - sigT > tol)
        {
            continue;
        }

        std::copy(t, t + 3, s);

        int nMC = 0, nT = 0;
        bool compEdge = false, extEdge = false;
        for (int i = 0; i < n; ++i)
        {
            (act[i] < 3 ? nMC : nT) += 1;
            compEdge = compEdge || act[i] == 1;
            extEdge = extEdge || act[i] == 2;
        }
        if (nMC && nT) return Return::tensionCorner;
        if (nT == 1) return Return::tensionPlane;
        if (nT == 2) return Return::tensionEdge;
        if (nT == 3) return Return::tensionApex;
        if (nMC == 3) return Return::apex;
        if (compEdge) return Return::compressionEdge;
        if (extEdge) return Return::extensionEdge;
        return Return::plane;
    }

    throw std::runtime_error
    (
        "mohrCoulomb: tension cut-off return found no active set"
    );
}


double MohrCoulomb::tensionYield(const Voigt sigma) const
{
    if (!tension_)
    {
        return -std::numeric_limits<double>::infinity();
    }
    double A[3][3], w[3], V[3][3];
    voigtToMatrix(sigma, A);
    symmetricEigen(A, w, V);
    return w[0] - sigmaT_;
}


double MohrCoulomb::yield(const Voigt sigma) const
{
    double A[3][3], w[3], V[3][3];
    voigtToMatrix(sigma, A);
    symmetricEigen(A, w, V);
    return (w[0] - w[2]) + (w[0] + w[2])*sinPhi_ - 2.0*c_*cosPhi_;
}

void MohrCoulomb::updateFluidity
(
    const Voigt deps,
    double dt,
    double g,
    const Voigt sigOld,
    const History& hOld,
    Voigt sigNew,
    History& hNew
) const
{
    Voigt trial;
    el_.addStress(sigOld, deps, trial);
    if (!nl_ || !nl_->update(el_, sigOld, trial, dt, g, hOld, sigNew, hNew))
    {
        updateRate(deps, dt, sigOld, hOld, sigNew, hNew);
    }
}


} // End namespace constitutive
} // End namespace cel
