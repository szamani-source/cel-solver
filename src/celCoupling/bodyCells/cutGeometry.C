/*---------------------------------------------------------------------------*\
    cutGeometry: exact edge, face and cell coverage by a closed triangulated
    surface (WP 4.2, decision 0042)
\*---------------------------------------------------------------------------*/

#include "cutGeometry.H"
#include "SortableList.H"

namespace Foam
{
namespace cutGeometry
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

//- Parameter s of the crossing of segment a + s d with triangle (p0, p1, p2),
//  false when none (Moller-Trumbore, edges and corners of the triangle
//  included with a small relative tolerance)
bool segTri
(
    const point& a,
    const vector& d,
    const point& p0,
    const point& p1,
    const point& p2,
    scalar& s
)
{
    const vector e1(p1 - p0);
    const vector e2(p2 - p0);
    const vector pv(d ^ e2);
    const scalar det = e1 & pv;
    const scalar scale = mag(e1)*mag(e2)*mag(d);
    if (mag(det) <= 1e-14*scale)
    {
        return false;               // parallel (or degenerate triangle)
    }
    const scalar inv = 1.0/det;
    const vector tv(a - p0);
    const scalar u = (tv & pv)*inv;
    const scalar tol = 1e-12;
    if (u < -tol || u > 1 + tol)
    {
        return false;
    }
    const vector qv(tv ^ e1);
    const scalar v = (d & qv)*inv;
    if (v < -tol || u + v > 1 + tol)
    {
        return false;
    }
    s = (e2 & qv)*inv;
    return s >= -tol && s <= 1 + tol;
}


//- Clip a polygon to the half-space (x - p).n <= 0
void clipPoly(DynamicList<point>& poly, const vector& n, const point& p)
{
    DynamicList<point> out(poly.size() + 2);
    const label m = poly.size();
    for (label i = 0; i < m; ++i)
    {
        const point& a = poly[i];
        const point& b = poly[(i + 1) % m];
        const scalar fa = (a - p) & n;
        const scalar fb = (b - p) & n;
        if (fa <= 0)
        {
            out.append(a);
        }
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0))
        {
            out.append(a + (fa/(fa - fb))*(b - a));
        }
    }
    poly.transfer(out);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

scalar edgeCovered
(
    const point& a,
    const point& b,
    const triSurface& s,
    const labelUList& tris
)
{
    const pointField& P = s.points();
    const vector d(b - a);

    DynamicList<scalar> sv;
    DynamicList<label> tv;         // +1 entering the body, -1 leaving
    for (const label ti : tris)
    {
        const labelledTri& f = s[ti];
        scalar t;
        if (segTri(a, d, P[f[0]], P[f[1]], P[f[2]], t))
        {
            const vector n((P[f[1]] - P[f[0]]) ^ (P[f[2]] - P[f[0]]));
            sv.append(min(max(t, scalar(0)), scalar(1)));
            tv.append((n & d) < 0 ? 1 : -1);
        }
    }
    if (sv.empty())
    {
        return -1;
    }

    // Sort, then merge equal crossings of the same kind (a segment through a
    // shared triangle edge or vertex meets every triangle there)
    SortableList<scalar> ss(sv);
    const labelList& order = ss.indices();
    DynamicList<scalar> s2;
    DynamicList<label> t2;
    forAll(order, k)
    {
        const scalar t = sv[order[k]];
        const label ty = tv[order[k]];
        if (s2.size() && mag(t - s2.last()) < 1e-10 && ty == t2.last())
        {
            continue;
        }
        s2.append(t);
        t2.append(ty);
    }

    // Inside at a if the first crossing leaves the body
    bool inside = (t2[0] == -1);
    scalar covered = 0;
    scalar last = 0;
    forAll(s2, k)
    {
        if (inside)
        {
            covered += s2[k] - last;
        }
        inside = (t2[k] == 1);
        last = s2[k];
    }
    if (inside)
    {
        covered += 1 - last;
    }
    return min(max(covered, scalar(0)), scalar(1));
}


scalar faceCovered
(
    const UList<point>& pts,
    const vector& e,
    const UList<scalar>& edgeFrac,
    const triSurface& s,
    const labelUList& tris
)
{
    const pointField& P = s.points();
    const label np = pts.size();

    point c(Zero);
    for (const point& p : pts)
    {
        c += p;
    }
    c /= np;

    // Twice the covered area, from Green's theorem
    scalar twoA = 0;

    // Covered parts of the edges, all along the edge direction
    for (label i = 0; i < np; ++i)
    {
        const vector d(pts[(i + 1) % np] - pts[i]);
        twoA += edgeFrac[i]*(((pts[i] - c) ^ d) & e);
    }

    // Cut lines of the triangles with the face plane, clipped to the face
    label nSeg = 0;
    for (const label ti : tris)
    {
        const labelledTri& f = s[ti];
        scalar h[3];
        for (label k = 0; k < 3; ++k)
        {
            h[k] = (P[f[k]] - pts[0]) & e;
        }
        point q[2];
        label nq = 0;
        for (label k = 0; k < 3 && nq < 2; ++k)
        {
            const label k1 = (k + 1) % 3;
            if ((h[k] < 0 && h[k1] > 0) || (h[k] > 0 && h[k1] < 0))
            {
                q[nq++] = P[f[k]] + (h[k]/(h[k] - h[k1]))*(P[f[k1]] - P[f[k]]);
            }
            else if (h[k] == 0 && h[k1] != 0)
            {
                q[nq++] = P[f[k]];
            }
        }
        if (nq < 2)
        {
            continue;
        }

        // Orientation: body on the left, tangent e ^ m (m: in-plane part of
        // the triangle's outward normal)
        const vector n((P[f[1]] - P[f[0]]) ^ (P[f[2]] - P[f[0]]));
        const vector m(n - (n & e)*e);
        const vector t(e ^ m);
        point p0 = q[0], p1 = q[1];
        if (((p1 - p0) & t) < 0)
        {
            std::swap(p0, p1);
        }

        // Clip to the face: keep (x - pts[i]).(d_i ^ e) <= 0 for every edge
        scalar s0 = 0, s1 = 1;
        const vector dq(p1 - p0);
        for (label i = 0; i < np && s0 < s1; ++i)
        {
            const vector on((pts[(i + 1) % np] - pts[i]) ^ e);
            const scalar f0 = (p0 - pts[i]) & on;
            const scalar df = dq & on;
            if (mag(df) < VSMALL)
            {
                if (f0 > 0) s1 = s0 - 1;   // wholly outside
                continue;
            }
            const scalar sc = -f0/df;
            if (df > 0) s1 = min(s1, sc); else s0 = max(s0, sc);
        }
        if (s1 <= s0)
        {
            continue;
        }
        const point a(p0 + s0*dq);
        const point b(p0 + s1*dq);
        twoA += ((a - c) ^ (b - c)) & e;
        ++nSeg;
    }

    if (nSeg == 0)
    {
        return -1;
    }
    return 0.5*twoA;
}


scalar cellCovered
(
    const UList<vector>& n,
    const UList<point>& p,
    const UList<scalar>& cA,
    const point& c,
    const triSurface& s,
    const labelUList& tris
)
{
    const pointField& P = s.points();
    scalar threeV = 0;

    // Triangles clipped to the cell; (x - c).n is constant on a plane
    for (const label ti : tris)
    {
        const labelledTri& f = s[ti];
        DynamicList<point> poly(8);
        poly.append(P[f[0]]);
        poly.append(P[f[1]]);
        poly.append(P[f[2]]);
        forAll(n, fi)
        {
            clipPoly(poly, n[fi], p[fi]);
            if (poly.size() < 3) break;
        }
        if (poly.size() < 3)
        {
            continue;
        }
        vector A(Zero);
        for (label i = 1; i + 1 < poly.size(); ++i)
        {
            A += 0.5*((poly[i] - poly[0]) ^ (poly[i + 1] - poly[0]));
        }
        threeV += (poly[0] - c) & A;
    }

    // Covered parts of the faces
    forAll(n, fi)
    {
        threeV += ((p[fi] - c) & n[fi])*cA[fi];
    }

    return threeV/3;
}

} // End namespace cutGeometry
} // End namespace Foam

// ************************************************************************* //
