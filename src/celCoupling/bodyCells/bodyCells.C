/*---------------------------------------------------------------------------*\
    bodyCells: exact cell and face coverage by rigid bodies (WP 4.2)
\*---------------------------------------------------------------------------*/

#include "bodyCells.H"
#include "cutGeometry.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"
#include "fvc.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bodyCells::bodyCells(const fvMesh& mesh, const rigidBodies& bodies)
:
    mesh_(mesh),
    bodies_(bodies),
    search_(bodies.size()),
    cellBb_(mesh.nCells()),
    alpha_
    (
        IOobject("alphaBody", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh, dimensionedScalar(dimless, Zero)
    ),
    open_
    (
        IOobject("bodyFaceOpen", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh, dimensionedScalar(dimless, 1)
    ),
    wallArea_
    (
        IOobject("wallArea", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh, dimensionedVector(dimArea, Zero)
    ),
    wallVelocity_
    (
        IOobject("wallVelocity", mesh.time().timeName(), mesh,
            IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh, dimensionedVector(dimVelocity, Zero)
    ),
    mixNbr_(mesh.nCells(), -1),
    cellBody_(mesh.nCells(), -1),
    nThin_(0),
    nSearched_(0),
    nCut_(0)
{
    forAll(bodies_, i)
    {
        search_.set(i, new triSurfaceSearch(bodies_[i].surface0()));
    }
    const pointField& pts = mesh_.points();
    scalar hSum = 0;
    forAll(cellBb_, celli)
    {
        cellBb_[celli] = treeBoundBox(pts, mesh_.cellPoints()[celli]);
        hSum += cmptMin(cellBb_[celli].span());
    }
    const scalar h =
        returnReduce(hSum, sumOp<scalar>())
       /max(returnReduce(mesh_.nCells(), sumOp<label>()), label(1));

    // Coarse facets on curved parts: neighbouring triangles longer than 4
    // cells along their shared edge meeting at 10 to 45 deg
    forAll(bodies_, i)
    {
        const triSurface& s = bodies_[i].surface0();
        const vectorField nf(s.faceNormals());
        const edgeList& edges = s.edges();
        const labelListList& ef = s.edgeFaces();
        const pointField& lp = s.localPoints();
        label nCoarse = 0;
        scalar longest = 0;
        forAll(edges, e)
        {
            if (ef[e].size() != 2) continue;
            const scalar len = edges[e].mag(lp);
            const scalar c = nf[ef[e][0]] & nf[ef[e][1]];
            const scalar ang =
                180/constant::mathematical::pi
               *Foam::acos(min(max(c, scalar(-1)), scalar(1)));
            if (len > 4*h && ang > 10 && ang < 45)
            {
                ++nCoarse;
                longest = max(longest, len);
            }
        }
        if (nCoarse)
        {
            WarningInFunction
                << "Body " << bodies_[i].name() << ": " << nCoarse
                << " facet edges up to " << longest/h << " cells long join "
                << "facets at 10 to 45 deg: a curved surface faceted at the "
                << "mesh scale; refine the STL" << endl;
        }
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::bodyCells::addBody
(
    const label bodyi,
    scalarField& cellV,
    scalarField& faceA,
    scalarField& bestCover,
    labelList& bestBody
) const
{
    const rigidBody& body = bodies_[bodyi];
    const triSurface& surf = body.surface0();
    const triSurfaceSearch& srch = search_[bodyi];
    const indexedOctree<treeDataTriSurface>& tree = srch.tree();

    const tensor& R = body.rotation();
    const point& c = body.centre();
    const point& c0 = body.centre0();
    auto toBody = [&](const point& x) { return c0 + (R.T() & (x - c)); };

    // Body bounding box in the world frame (corners of the frame-0 box)
    const boundBox bb0(surf.points(), false);
    boundBox bbW;
    const pointField corners(treeBoundBox(bb0).points());
    for (const point& q : corners)
    {
        bbW.add(body.transform(q));
    }
    const scalar tolW = 1e-9*mag(bbW.span());
    bbW.inflate(1e-9);

    const pointField& mp = mesh_.points();
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();
    const labelUList& own = mesh_.faceOwner();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField magSf(mag(mesh_.faceAreas()));

    // Candidate cells and their faces
    DynamicList<label> cand;
    forAll(cellBb_, celli)
    {
        if (cellBb_[celli].overlaps(bbW))
        {
            cand.append(celli);
        }
    }
    const_cast<bodyCells&>(*this).nSearched_ += cand.size();

    bitSet faceDone(mesh_.nFaces());
    Map<scalar> faceCov;           // covered area per candidate face

    auto boxOf = [&](const labelUList& pl)
    {
        boundBox b;
        for (const label pi : pl)
        {
            b.add(toBody(mp[pi]));
        }
        treeBoundBox tb(b);
        tb.inflate(1e-6);
        return tb;
    };

    // 1. Faces: covered fraction of each edge, then covered area
    DynamicList<label> pendFace;        // faces whose edges need inside tests
    DynamicList<point> pendPoint;
    DynamicList<label> pendEdge;        // edge index within the face
    HashTable<scalarField, label> edgeFrac;
    for (const label celli : cand)
    {
        for (const label facei : cells[celli])
        {
            if (faceDone.test(facei)) continue;
            faceDone.set(facei);
            const face& f = faces[facei];
            const labelList tris(tree.findBox(boxOf(f)));
            scalarField ef(f.size(), 0);
            forAll(f, k)
            {
                const point a(toBody(mp[f[k]]));
                const point b(toBody(mp[f.nextLabel(k)]));
                const scalar r = cutGeometry::edgeCovered(a, b, surf, tris);
                if (r < 0)
                {
                    pendFace.append(facei);
                    pendEdge.append(k);
                    pendPoint.append(0.5*(a + b));
                }
                ef[k] = max(r, scalar(0));
            }
            edgeFrac.insert(facei, ef);
        }
    }
    if (pendPoint.size())
    {
        const boolList in(srch.calcInside(pointField(pendPoint)));
        forAll(in, i)
        {
            if (in[i]) edgeFrac[pendFace[i]][pendEdge[i]] = 1;
        }
    }
    forAllConstIters(edgeFrac, iter)
    {
        const label facei = iter.key();
        const scalarField& ef = iter.val();
        const face& f = faces[facei];
        pointField fp(f.size());
        forAll(f, k)
        {
            fp[k] = toBody(mp[f[k]]);
        }
        const vector e(R.T() & (Sf[facei]/magSf[facei]));
        const labelList tris(tree.findBox(boxOf(f)));
        scalar A = cutGeometry::faceCovered(fp, e, ef, surf, tris);
        if (A < 0)
        {
            A = average(ef)*magSf[facei];   // all 0 or all 1 when not cut
        }
        A = min(max(A, scalar(0)), magSf[facei]);
        faceCov.insert(facei, A);
        faceA[facei] += A;
    }

    // 2. Cells: covered volume
    DynamicList<label> pendCell;
    DynamicList<point> pendC;
    const scalarField& V = mesh_.V();
    const vectorField& Cc = mesh_.C();
    for (const label celli : cand)
    {
        const cell& cf = cells[celli];
        const labelList tris(tree.findBox(boxOf(cf.labels(faces))));
        scalar cov;
        if (tris.empty())
        {
            scalar sumA = 0, sumMag = 0;
            for (const label facei : cf)
            {
                sumA += faceCov[facei];
                sumMag += magSf[facei];
            }
            if (sumA > 0.5*sumMag) cov = V[celli];
            else if (sumA > 0) cov = (sumA/sumMag)*V[celli];
            else
            {
                cov = 0;   // may still be inside a body (no face cut): test
                pendCell.append(celli);
                pendC.append(toBody(Cc[celli]));
            }
        }
        else
        {
            List<vector> n(cf.size());
            List<point> p(cf.size());
            List<scalar> cA(cf.size());
            forAll(cf, k)
            {
                const label facei = cf[k];
                const scalar sgn = (own[facei] == celli) ? 1 : -1;
                n[k] = R.T() & (sgn*Sf[facei]/magSf[facei]);
                p[k] = toBody(mesh_.faceCentres()[facei]);
                cA[k] = faceCov[facei];
            }
            cov = cutGeometry::cellCovered
            (
                n, p, cA, toBody(Cc[celli]), surf, tris
            );
            cov = min(max(cov, scalar(0)), V[celli]);
        }
        cellV[celli] += cov;
        if (cov > bestCover[celli])
        {
            bestCover[celli] = cov;
            bestBody[celli] = bodyi;
        }
    }
    if (pendC.size())
    {
        // Cells whose faces are all open but which hold the whole body
        // inside are impossible for a closed body larger than a cell; a cell
        // with no faces cut and its centre inside the body is covered.
        const boolList in(srch.calcInside(pointField(pendC)));
        forAll(in, i)
        {
            if (in[i])
            {
                cellV[pendCell[i]] += V[pendCell[i]];
                bestCover[pendCell[i]] = V[pendCell[i]];
                bestBody[pendCell[i]] = bodyi;
            }
        }
    }
    (void)tolW;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bodyCells::update()
{
    nSearched_ = 0;
    nCut_ = 0;

    // Deformable bodies (WP 5.3): the surface moved, rebuild its search
    forAll(bodies_, i)
    {
        if (bodies_[i].deformable())
        {
            search_.set(i, new triSurfaceSearch(bodies_[i].surface0()));
        }
    }

    scalarField cellV(mesh_.nCells(), 0);
    scalarField faceA(mesh_.nFaces(), 0);
    scalarField bestCover(mesh_.nCells(), 0);
    labelList bestBody(mesh_.nCells(), -1);

    forAll(bodies_, bodyi)
    {
        addBody(bodyi, cellV, faceA, bestCover, bestBody);
    }

    const scalarField& V = mesh_.V();
    const scalarField magSf(mag(mesh_.faceAreas()));
    const vectorField& Sf = mesh_.faceAreas();
    const labelUList& own = mesh_.faceOwner();
    const labelUList& nei = mesh_.faceNeighbour();

    scalarField& a = alpha_.primitiveFieldRef();
    a = min(cellV/V, scalar(1));

    scalarField openAll(mesh_.nFaces());
    forAll(openAll, facei)
    {
        openAll[facei] = max(1 - faceA[facei]/magSf[facei], scalar(0));
    }
    open_.primitiveFieldRef() =
        SubField<scalar>(openAll, mesh_.nInternalFaces());
    forAll(open_.boundaryField(), patchi)
    {
        fvsPatchScalarField& pf = open_.boundaryFieldRef()[patchi];
        const label start = mesh_.boundaryMesh()[patchi].start();
        forAll(pf, i)
        {
            pf[i] = openAll[start + i];
        }
    }

    // Wall face: minus the open face areas times outward normals
    vectorField& W = wallArea_.primitiveFieldRef();
    W = Zero;
    forAll(own, facei)
    {
        const vector s(openAll[facei]*Sf[facei]);
        W[own[facei]] -= s;
        if (facei < mesh_.nInternalFaces())
        {
            W[nei[facei]] += s;
        }
    }
    vectorField& U = wallVelocity_.primitiveFieldRef();
    U = Zero;
    const vectorField& Cc = mesh_.C();
    mixNbr_ = -1;
    const cellList& cells = mesh_.cells();
    forAll(W, celli)
    {
        const bool cut = a[celli] > 1e-12 && a[celli] < 1 - 1e-12;
        if (!cut)
        {
            W[celli] = Zero;     // round-off of a closed cell
        }
        const scalar mW = mag(W[celli]);
        if (bestBody[celli] >= 0)
        {
            const rigidBody& body = bodies_[bestBody[celli]];
            if (cut && body.deformable() && mW >= VSMALL)
            {
                // Mean over the wall face (decision 0073)
                U[celli] = Zero;
                for
                (
                    const point& p
                  : rigidBody::faceSamples
                    (
                        Cc[celli], cellBb_[celli].span(), W[celli]/mW
                    )
                )
                {
                    U[celli] += 0.25*body.velocityAt(p);
                }
            }
            else
            {
                U[celli] = body.velocityAt(Cc[celli]);
            }
        }
        if (!cut) continue;
        ++nCut_;
        if (mW < VSMALL) continue;
        const vector nb(W[celli]/mW);
        scalar best = 0;
        for (const label facei : cells[celli])
        {
            const scalar sgn = (own[facei] == celli) ? 1 : -1;
            const scalar score = -openAll[facei]*((sgn*Sf[facei]) & nb);
            if (score > best)
            {
                best = score;
                mixNbr_[celli] =
                    facei < mesh_.nInternalFaces()
                  ? (own[facei] == celli ? nei[facei] : own[facei])
                  : -1;
            }
        }
    }
    cellBody_ = bestBody;

    // Thin parts: a cut cell with no face neighbour deeper in a body (on
    // this processor): the body ends within the cell on every side
    nThin_ = 0;
    forAll(a, celli)
    {
        if (a[celli] <= 1e-12 || a[celli] >= 1 - 1e-12) continue;
        bool backed = false;
        for (const label facei : cells[celli])
        {
            if (facei >= mesh_.nInternalFaces()) continue;
            const label other =
                own[facei] == celli ? nei[facei] : own[facei];
            if (a[other] > a[celli] + 1e-9)
            {
                backed = true;
                break;
            }
        }
        if (!backed) ++nThin_;
    }
    if (nThin() > 0)
    {
        WarningInFunction
            << nThin() << " cut cells have no neighbour deeper in a body: a body "
            << "is thinner than about two cells there (decision 0042); "
            << "material may pass through. Thicken it to at least two cells"
            << endl;
    }

    alpha_.correctBoundaryConditions();
    wallArea_.correctBoundaryConditions();
    wallVelocity_.correctBoundaryConditions();
}


Foam::label Foam::bodyCells::nThin() const
{
    return returnReduce(nThin_, sumOp<label>());
}


void Foam::bodyCells::write() const
{
    const word tn(mesh_.time().timeName());
    auto writeAt = [&](const auto& fld)
    {
        typename std::remove_cv_t<std::remove_reference_t<decltype(fld)>>
            out(IOobject(fld.name(), tn, mesh_, IOobject::NO_READ,
                IOobject::NO_WRITE, IOobject::NO_REGISTER), fld);
        out.write();
    };
    writeAt(alpha_);
    writeAt(open_);
    writeAt(wallArea_);
    writeAt(wallVelocity_);

    volScalarField mix
    (
        IOobject("mixNeighbour", tn, mesh_, IOobject::NO_READ,
            IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh_, dimensionedScalar(dimless, -1)
    );
    forAll(mixNbr_, celli)
    {
        mix[celli] = mixNbr_[celli];
    }
    mix.write();
}


// ************************************************************************* //
