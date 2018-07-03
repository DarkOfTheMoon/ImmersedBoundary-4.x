/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "immersedBoundaryForces.H"
#include "immersedBoundaryFvPatch.H"
#include "immersedBoundaryFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(immersedBoundaryForces, 0);

    addToRunTimeSelectionTable
    (
        functionObject, immersedBoundaryForces, dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::immersedBoundaryForces::immersedBoundaryForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces
    (
        name,
        obr,
        dict
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::immersedBoundaryForces::~immersedBoundaryForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::immersedBoundaryForces::calcForcesMoment()
{
    force_[0] = vector::zero;
    force_[1] = vector::zero;
    force_[2] = vector::zero;

    moment_[0] = vector::zero;
    moment_[1] = vector::zero;
    moment_[2] = vector::zero;


    if (directForceDensity_)
    {
        const volVectorField& fD = obr_.lookupObject<volVectorField>(fDName_);

        const fvMesh& mesh = fD.mesh();

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            // Check and cast into immersed boundary type
            if
            (
                isA<immersedBoundaryFvPatchVectorField>
                (
                    fD.boundaryField()[patchI]
                )
            )
            {
                // Found immersed boundary patch and field.
                // Cast into immersed boundary type
                const immersedBoundaryFvPatch& ibPatch =
                    refCast<const immersedBoundaryFvPatch>
                    (
                        mesh.boundary()[patchI]
                    );

                const immersedBoundaryFvPatchVectorField fDpatch =
                    refCast<const immersedBoundaryFvPatchVectorField>
                    (
                        fD.boundaryField()[patchI]
                    );

                // Get face area vectors for triangles
                const vectorField& Sfb = ibPatch.triSf();
                scalarField sA = mag(Sfb);

                // Calculate distance for triangles
                vectorField Md = ibPatch.triCf() - coordSys_.origin();

                // Normal force =
                // surfaceNormal*(surfaceUnitNormal & forceDensity)
                // The first operation will be done on ibPoints, the data will
                // then be distributed onto the ib surface
                // for surface integration
                vectorField fN =
                    Sfb*
                    ibPatch.toTriFaces
                    (
                        ibPatch.ibNormals() & fDpatch.ibValue()
                    );

                force_[0] += sum(fN);
                moment_[0] += sum(Md ^ fN);

                // Tangential force (total force minus normal fN)
                vectorField fT = sA*fDpatch.triValue() - fN;

                force_[1] += sum(fT);
                moment_[1] += sum(Md ^ fT);
            }
        }
    }
    else
    {
        const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        const fvMesh& mesh = U.mesh();

        volSymmTensorField stress = devRhoReff();

        // Scale pRef by density for incompressible simulations
        scalar pRef = pRef_/rho(p);

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchI = iter.key();

            // Check and cast into immersed boundary type
            if
            (
                isA<immersedBoundaryFvPatchVectorField>
                (
                    U.boundaryField()[patchI]
                )
            )
            {
                // Found immersed boundary patch and field.
                // Cast into immersed boundary type
                const immersedBoundaryFvPatch& ibPatch =
                    refCast<const immersedBoundaryFvPatch>
                    (
                        mesh.boundary()[patchI]
                    );

                const immersedBoundaryFvPatchSymmTensorField stressPatch =
                    refCast<const immersedBoundaryFvPatchSymmTensorField>
                    (
                        stress.boundaryField()[patchI]
                    );

                const immersedBoundaryFvPatchScalarField pPatch =
                    refCast<const immersedBoundaryFvPatchScalarField>
                    (
                        p.boundaryField()[patchI]
                    );

                // Get face area vectors for triangles
                const vectorField& Sfb = ibPatch.triSf();
                scalarField sA = mag(Sfb);

                // Calculate distance for triangles
                vectorField Md = ibPatch.triCf() - coordSys_.origin();

                // Pressure force is an integral of interpolated pressure
                // on triangular faces
                vectorField pf = Sfb*(pPatch.triValue() - pRef);

                force_[0] += rho(p)*sum(pf);
                moment_[0] += rho(p)*sum(Md ^ pf);

                // Shear force is dotted with a normal in the IB point
                // and integrated on triangular faces
                vectorField vf =
                    sA*
                    ibPatch.toTriFaces
                    (
                        ibPatch.ibNormals() & stress[patchI]
                    );

                force_[1] += sum(vf);
                moment_[1] += sum(Md ^ vf);
            }
        }
    }

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


// ************************************************************************* //
