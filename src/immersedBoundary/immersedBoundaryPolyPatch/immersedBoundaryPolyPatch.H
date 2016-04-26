/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::immersedBoundaryPolyPatch

Description
    Immersed boundary patch

SourceFiles
    immersedBoundaryPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryPolyPatch_H
#define immersedBoundaryPolyPatch_H

#include "polyPatch.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryPolyPatch
:
    public polyPatch
{
    // Private data

        //- Triangular surface representing immersed boundary.
        //  Name of tri surface will be identical to the name of the patch
        triSurfaceMesh surface_;

        //- Internal or external flow calculation
        Switch internalFlow_;

public:

    //- Runtime type information
    TypeName("immersedBoundary");


    // Constructors

        //- Construct from components
        immersedBoundaryPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct from dictionary
        immersedBoundaryPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~immersedBoundaryPolyPatch();


    // Member Functions

        // Access

            //- Return immersed boundary surface mesh
            const triSurface& surface() const
            {
                return surface_;
            }

            //- Return true if solving for flow inside the immersed boundary
            bool internalFlow() const
            {
                return internalFlow_;
            }

            //- Return triSurface search object
            const triSurfaceSearch& triSurfSearch() const
            {
                return surface_;
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
