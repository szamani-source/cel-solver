/*---------------------------------------------------------------------------*\
Class
    Foam::bodyCells

Description
    Which Eulerian cells and faces the rigid bodies cover, exactly, at their
    current pose (WP 4.2, decision 0042). Per cell: body fraction
    (alphaBody), the wall face of a cut cell as area times unit normal into
    the body (wallArea = minus the sum of open face areas times outward
    normals, so a uniform stress gives no force), the body velocity there
    (wallVelocity, at the cell centre) and the neighbour a small cut cell is
    mixed with (mixNeighbour: across the face with the largest open area
    facing away from the body; -1 when that face is a boundary). Per face:
    open fraction (bodyFaceOpen, 1 where no body).

    Each body's surface is searched in its own frame (the STL as read, one
    search tree for the whole run); the mesh's cells near the body are
    mapped into that frame. Only cells whose bounding box meets the body's
    are looked at. Bodies must be closed surfaces with outward normals and
    must not overlap each other. The mesh must be made of convex cells with
    planar faces (Cartesian hexes in practice).

    Warnings (decision 0042 addendum, WP 4.3a): at construction, facets
    longer than 4 cells meeting at 10 to 45 deg (a curved surface faceted at
    the mesh scale); at each update, cut cells with no face neighbour more
    covered than themselves (a body thinner than about two cells: material
    may pass).

SourceFiles
    bodyCells.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_bodyCells_H
#define Foam_bodyCells_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "rigidBodies.H"
#include "triSurfaceSearch.H"

namespace Foam
{

class bodyCells
{
    // Private Data

        const fvMesh& mesh_;
        const rigidBodies& bodies_;

        //- Search tree per body, on the surface as read (body frame)
        PtrList<triSurfaceSearch> search_;

        //- Cell bounding boxes (world frame; the mesh does not move)
        List<treeBoundBox> cellBb_;

        volScalarField alpha_;
        surfaceScalarField open_;
        volVectorField wallArea_;
        volVectorField wallVelocity_;
        labelList mixNbr_;

        //- Body covering most of each cell (-1: none)
        labelList cellBody_;

        //- Cut cells with no face neighbour more covered than themselves
        //  in the last update (a body thinner than about two cells there)
        label nThin_;

        //- Number of cells looked at in the last update, and cut cells
        label nSearched_;
        label nCut_;


    // Private Member Functions

        //- Add one body's covered cell volume and face areas
        void addBody
        (
            const label bodyi,
            scalarField& cellV,
            scalarField& faceA,
            scalarField& bestCover,
            labelList& bestBody
        ) const;


public:

    // Constructors

        bodyCells(const fvMesh& mesh, const rigidBodies& bodies);


    // Member Functions

        //- Recompute everything at the bodies' current poses
        void update();

        //- Write alphaBody, bodyFaceOpen, wallArea, wallVelocity, mixNeighbour
        void write() const;

        const volScalarField& alphaBody() const noexcept { return alpha_; }
        const surfaceScalarField& faceOpen() const noexcept { return open_; }
        const volVectorField& wallArea() const noexcept { return wallArea_; }
        const volVectorField& wallVelocity() const noexcept
        {
            return wallVelocity_;
        }
        const labelList& mixNeighbour() const noexcept { return mixNbr_; }
        const labelList& cellBody() const noexcept { return cellBody_; }

        //- Cut cells without a more covered face neighbour (all processors)
        label nThin() const;

        label nSearched() const noexcept { return nSearched_; }
        label nCut() const noexcept { return nCut_; }
};


} // End namespace Foam

#endif
