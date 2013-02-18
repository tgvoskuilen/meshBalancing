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
label Foam::dynamicRefineBalancedFvMesh::topParentID(label p)
{
    nextP = meshCutter().history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParent(nextP);
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
    //Part 1 - Copied from dynamicRefineFvMesh::update
    // If dynamicRefineFvMesh put this code in a separate "doUpdate" function
    // and set update so it just called doUpdate(), then this code would not
    // have to be repeated. It could also just call doUpdate().
    
    

    // Re-read dictionary. Choosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );

    label refineInterval = readLabel(refineDict.lookup("refineInterval"));

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        changing(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorIn("dynamicRefineFvMesh::update()")
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }




    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        label maxCells = readLabel(refineDict.lookup("maxCells"));

        if (maxCells <= 0)
        {
            FatalErrorIn("dynamicRefineFvMesh::update()")
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        label maxRefinement = readLabel(refineDict.lookup("maxRefinement"));

        if (maxRefinement <= 0)
        {
            FatalErrorIn("dynamicRefineFvMesh::update()")
                << "Illegal maximum refinement level " << maxRefinement << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        const word fieldName(refineDict.lookup("field"));

        const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

        const scalar lowerRefineLevel =
            readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel =
            readScalar(refineDict.lookup("upperRefineLevel"));
        const scalar unrefineLevel =
            readScalar(refineDict.lookup("unrefineLevel"));
        const label nBufferLayers =
            readLabel(refineDict.lookup("nBufferLayers"));

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(nCells());

        if (globalData().nTotalCells() < maxCells)
        {
            // Determine candidates for refinement (looking at field only)
            selectRefineCandidates
            (
                lowerRefineLevel,
                upperRefineLevel,
                vFld,
                refineCell
            );

            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells,
                    maxRefinement,
                    refineCell
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                // Refine/update mesh and map fields
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                // Update refineCell. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, cellI)
                    {
                        label oldCellI = cellMap[cellI];

                        if (oldCellI < 0)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else if (reverseCellMap[oldCellI] != cellI)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else
                        {
                            newRefineCell.set(cellI, refineCell.get(oldCellI));
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }

                // Extend with a buffer layer to prevent neighbouring points
                // being unrefined.
                for (label i = 0; i < nBufferLayers; i++)
                {
                    extendMarkedCells(refineCell);
                }

                hasChanged = true;
            }
        }


        {
            // Select unrefineable points that are not marked in refineCell
            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    minCellField(vFld)
                )
            );

            label nSplitPoints = returnReduce
            (
                pointsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitPoints > 0)
            {
                // Refine/update mesh
                unrefine(pointsToUnrefine);

                hasChanged = true;
            }
        }


        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occassionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<refinementHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;
              
    }

    changing(hasChanged);
    
    
    // LOAD BALANCING SECTION
    // TODO: Read 'balancingInterval' 'enableBalancing' and 'maxImbalance'
    //       from dictionary
    if ( time().timeIndex() > 0 && Pstream::parRun() && hasChanged )
    {
        //First determine current level of imbalance
        label nGlobalCells = globalData().nTotalCells();
        
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/nGlobalCells;
        
        bool imbalanced = (maxImbalance > 0.1);
        
        Info<< "Maximum imbalance = " << 100*maxImbalance << " %" << endl;
        
        //If imbalanced, construct grouping fields
        if( imbalanced )
        {
            Info<< "Re-balancing problem" << endl;
            
  //VERSION 1 - TO BE REPLACED
            //calc weights for distribution that doesn't break refinement
            scalar maxV = gMax(V());
            scalar dx = Foam::pow(maxV,1.0/3.0); //Assume dx=dy=dz
            scalar dy = dx; //TODO
            scalar dz = dx; //TODO
            
            point minCC = bounds().min() + 0.5*vector(dx,dy,dz);
            Tensor<scalar> T(0.0);
            T.xx() = 1.0/dx;
            T.yy() = 1.0/dy;
            T.zz() = 1.0/dz;
            
            //Digitize mesh cell centroid coordinates to integer-scaled values
            pointField cellRelCoord = (T & cellCentres()) - minCC;

            // Round to nearest integer values
            forAll(cellRelCoord, cellI)
            {
                point& c = cellRelCoord[cellI];

                c.x() = std::floor(c.x() + 0.5);
                c.y() = std::floor(c.y() + 0.5);
                c.z() = std::floor(c.z() + 0.5);
            }
            
            // Determine span of local mesh
            //label nx = label( mesh.bounds().span().x()/dx + 0.5 );
            label ny = label( bounds().span().y()/dy + 0.5 );
            label nz = label( bounds().span().z()/dz + 0.5 );

            //Assign indices to cells based on cellRelCoord
            Map<label> coarseIDmap(nCells());
            labelList globalIndices(nCells(),0);
            label localID = 0;
            
            forAll(globalIndices, cellI)
            {
                point& c = cellRelCoord[cellI];
                
                label i = label(c.x()+0.5);
                label j = label(c.y()+0.5);
                label k = label(c.z()+0.5);
                
                globalIndices[cellI] = j + ny*k + ny*nz*i;
                if( coarseIDmap.insert(globalIndices[cellI], localID) )
                {
                    ++localID;
                }
            }
            
            label nCoarse = localID;
            
            // Convert to local, sequential indexing
            labelList localIndices(nCells(),0);
            forAll(globalIndices, cellI)
            {
                localIndices[cellI] = coarseIDmap[globalIndices[cellI]];
            }
            
            Pout << "Proc has " << nCoarse << " blocks" << endl;
            
            
            
            //localIndices = fineToCoarse
            // now make coarsePoints and coarseWeights
            pointField coarsePoints(nCoarse,vector::zero);
            scalarField coarseWeights(nCoarse,0.0);
            
            forAll(localIndices, cellI)
            {
                point& c = cellRelCoord[cellI];
                
                label i = label(c.x()+0.5);
                label j = label(c.y()+0.5);
                label k = label(c.z()+0.5);
                
                label& li = localIndices[cellI];
                
                coarsePoints[li] = vector
                (
                    (i+0.5)*dx,
                    (j+0.5)*dy,
                    (k+0.5)*dz
                ) + minCC; //Not sure if I need to add this back or not
                coarseWeights[li] += 1.0;
            }
            
            
            /*
            
   //VERSION 2 - not tested yet, much simpler though
            
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
            
            Pout << "Proc has " << nCoarse << " blocks" << endl;
            */
            
            
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
