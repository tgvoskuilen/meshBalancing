/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "dynamicRefineBalancedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineBalancedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineBalancedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::label Foam::dynamicRefineBalancedFvMesh::topParentID(label p)
{
    label nextP = meshCutter().history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::dynamicRefineBalancedFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::~dynamicRefineBalancedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedFvMesh::update()
{
    //Part 1 - Call normal update from dynamicRefineFvMesh
    bool hasChanged = dynamicRefineFvMesh::update();
    
    // Part 2 - Load Balancing
    bool enableBalancing = true;     //TODO: Read from dictionary
    scalar allowableImbalance = 0.1; //TODO: Read from dictionary
    if ( Pstream::parRun() && hasChanged && enableBalancing )
    {
        //First determine current level of imbalance
        label nGlobalCells = globalData().nTotalCells();
        
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/nGlobalCells;
        
        Info<< "Maximum imbalance = " << 100*maxImbalance << " %" << endl;
        
        //If imbalanced, construct grouping fields and re-balance
        if( maxImbalance > allowableImbalance )
        {
            Info<< "Re-balancing problem" << endl;
                        
            const labelIOList& cellLevel = meshCutter().cellLevel();
            Map<label> coarseIDmap(nCells());
            labelList uniqueIndex(nCells(),0);
            
            label nCoarse = 0;

            forAll(cells(), cellI)
            {
                if( cellLevel[cellI] > 0 )
                {
                    uniqueIndex[cellI] = nCells() + topParentID
                    (
                        meshCutter().history().parentIndex(cellI)
                    );
                }
                else
                {
                    uniqueIndex[cellI] = cellI;
                }
                
                if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
                {
                    ++nCoarse;
                }
            }
            
            // Convert to local, sequential indexing and calculate coarse
            // points and weights
            labelList localIndices(nCells(),0);
            pointField coarsePoints(nCoarse,vector::zero);
            scalarField coarseWeights(nCoarse,0.0);
            
            forAll(uniqueIndex, cellI)
            {
                localIndices[cellI] = coarseIDmap[uniqueIndex[cellI]];
                
                label w = (1 << (3*cellLevel[cellI]));
                coarseWeights[localIndices[cellI]] += 1.0;
                coarsePoints[localIndices[cellI]] += C()[cellI]/w;
            }
            
            //Set up decomposer                
            autoPtr<decompositionMethod> decomposer
            (
                decompositionMethod::New
                (
                    IOdictionary
                    (
                        IOobject
                        (
                            "balanceParDict",
                            time().system(),
                            *this,
                            IOobject::MUST_READ_IF_MODIFIED,
                            IOobject::NO_WRITE
                        )
                    )
                )
            );
            
            labelList finalDecomp = decomposer().decompose
            (
                *this, 
                localIndices,
                coarsePoints,
                coarseWeights
            );
        
            scalar tolDim = globalMeshData::matchTol_ * bounds().mag();
            
            fvMeshDistribute distributor(*this, tolDim);
            
            autoPtr<mapDistributePolyMesh> map =
                  distributor.distribute(finalDecomp);
                  
            meshCutter_.distribute(map);
        }
    }

    return hasChanged;
}


// ************************************************************************* //
