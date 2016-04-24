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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryPolyPatch.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, immersedBoundaryPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        immersedBoundaryPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryPolyPatch::makeTriSurfSearch() const
{
    if (debug)
    {
        Info<< "void immersedBoundaryPolyPatch::makeTriSurfSearch() const : "
            << "creating triSurface search algorithm"
            << endl;
    }

    if (triSurfSearchPtr_)
    {
        FatalErrorIn(__PRETTY_FUNCTION__)
            << "triSurface search algorithm already exist"
            << abort(FatalError);
    }

    triSurfSearchPtr_ = new triSurfaceSearch(ibMesh_.surface());
}


void Foam::immersedBoundaryPolyPatch::clearOut()
{
    deleteDemandDrivenData(triSurfSearchPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    ibMesh_
    (
        IOobject
        (
            name  + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(false),
    triSurfSearchPtr_(NULL)
{}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    ibMesh_
    (
        IOobject
        (
            name  + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(dict.lookup("internalFlow")),
    triSurfSearchPtr_(NULL)
{
    if (size() > 0)
    {
        FatalIOErrorIn(__PRETTY_FUNCTION__, dict)
            << "Faces detected in the immersedBoundaryPolyPatch.  "
            << "This is not allowed: please make sure that the patch size "
            << "equals zero."
            << abort(FatalIOError);
    }
}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(pp.internalFlow()),
    triSurfSearchPtr_(NULL)
{}


Foam::immersedBoundaryPolyPatch::immersedBoundaryPolyPatch
(
    const immersedBoundaryPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    ibMesh_
    (
        IOobject
        (
            pp.name() + ".ftr",
            bm.mesh().time().constant(), // instance
            "triSurface",                // local
            bm.mesh(),                   // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    internalFlow_(pp.internalFlow()),
    triSurfSearchPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::immersedBoundaryPolyPatch::~immersedBoundaryPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::triSurfaceSearch&
Foam::immersedBoundaryPolyPatch::triSurfSearch() const
{
    if (!triSurfSearchPtr_)
    {
        makeTriSurfSearch();
    }

    return *triSurfSearchPtr_;
}


void Foam::immersedBoundaryPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("internalFlow") << internalFlow_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
