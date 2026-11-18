/*---------------------------------------------------------------------------*\
Class
    Foam::rigidBodies

Description
    The rigid bodies of a case (WP 4.1, decision 0041), read from
    constant/rigidBodies:

        bodies
        {
            <name> { surface <file>; centre ...; velocity ...; omega ...; }
            ...
        }

    See rigidBody.H for the entries. A case without constant/rigidBodies
    has no bodies and nothing changes.

    Output (master processor only):
      postProcessing/rigidBodies/<name>.dat   one line per write time:
          t  centre (3)  rotation tensor (9, row by row)  v (3)  omega (3)
      postProcessing/rigidBodies/<time>/<name>.stl   the moved surface

SourceFiles
    rigidBodies.C
\*---------------------------------------------------------------------------*/

#ifndef Foam_rigidBodies_H
#define Foam_rigidBodies_H

#include "rigidBody.H"
#include "PtrList.H"
#include "OFstream.H"

namespace Foam
{

class rigidBodies
:
    public PtrList<rigidBody>
{
    // Private Data

        const Time& runTime_;

        //- One history file per body (master only)
        PtrList<OFstream> files_;


public:

    // Constructors

        //- Read constant/rigidBodies if present; bodies are moved to the
        //  current time
        explicit rigidBodies(const Time& runTime);


    // Member Functions

        //- Move every body to the current time
        void move();

        //- Write the pose of every body and its moved surface
        void write();
};


} // End namespace Foam

#endif
