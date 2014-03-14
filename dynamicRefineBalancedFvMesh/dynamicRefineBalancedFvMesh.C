/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
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
#include "fvCFD.H"

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


Foam::Pair<scalar> Foam::dynamicRefineBalancedFvMesh::readRefinementPoints()
{
    dictionary refineDict
}
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
    
    return Pair<scalar>
    (
        readScalar(refineDict.lookup("unrefineLevel")), 
        readScalar(refineDict.lookup("lowerRefineLevel"))
    );
}


void Foam::dynamicRefineBalancedFvMesh::updateRefinementField()
{
    Info<< "Calculating internal refinement field" << endl;
    
    volScalarField& intRefFld = *internalRefinementFieldPtr_;
    
    // Set the internal refinement field to zero to start with
    intRefFld = dimensionedScalar("zero",dimless,0.0);
    
    // Get the cell level field from dynamicRefineFvMesh
    const labelList& cellLevel = meshCutter().cellLevel();
    
    // Read the points at which refinement and unrefinement occur from the
    // dynamicMeshDict entries
    Pair<scalar> refinePoints = readRefinementPoints();
    
    // First gradients
    List<word> gradFieldNames = gradFields_.toc();
    
    Field<scalar> cubeRtV = Foam::pow(this->V(),1.0/3.0);
    Field<scalar> refFld(nCells(),0.0);
    
    forAll(gradFieldNames, i)
    {
        word fldName = gradFieldNames[i];
        scalar wgt = gradFields_[fldName].first();
        label maxLevel = static_cast<label>(gradFields_[fldName].second()+0.5);
        
        const volScalarField& fld = this->lookupObject<volScalarField>(fldName);
        
        refFld = wgt * mag(fvc::grad(fld)) * cubeRtV;
        
        // Limit the value of refFld based on its max level
        forAll(refFld, cellI)
        {
            if( cellLevel[cellI] >= maxLevel )
            {
                refFld[cellI] = min
                (
                    refFld[cellI],
                    0.5*(refinePoints.first() + refinePoints.second())
                );
            }
        }
 
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    // Then curls
    List<word> curlFieldNames = curlFields_.toc();
    
    
    forAll(curlFieldNames, i)
    {
        word fldName = curlFieldNames[i];
        scalar wgt = curlFields_[fldName].first();
        label maxLevel = static_cast<label>(curlFields_[fldName].second()+0.5);
        
        const volVectorField& fld = this->lookupObject<volVectorField>(fldName);
        
        refFld = wgt * mag(fvc::curl(fld)) * cubeRtV;
        
        // Limit the value of refFld based on its max level
        forAll(refFld, cellI)
        {
            if( cellLevel[cellI] >= maxLevel )
            {
                refFld[cellI] = min
                (
                    refFld[cellI],
                    0.5*(refinePoints.first() + refinePoints.second())
                );
            }
        }
        
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    // The set refinement physical regions (force the mesh to stay refined
    // near key features)
    forAll(refinedRegions_, regionI)
    {
        const entry& region = refinedRegions_[regionI];
        
        autoPtr<topoSetSource> source = 
            topoSetSource::New(region.keyword(), *this, region.dict());
    
        refFld *= 0.0;
        
        cellSet selectedCellSet
        (
            *this,
            "cellSet",
            nCells()/10+1 //estimate
        );
        
        source->applyToSet
        (
            topoSetSource::NEW,
            selectedCellSet
        );
        
        const labelList cells = selectedCellSet.toc();
        
        label minLevel = readLabel(region.dict().lookup("minLevel"));
        
        forAll(cells, i)
        {
            const label& cellI = cells[i];
            
            if( cellLevel[cellI] < minLevel ) //force refinement
            {
                refFld[cellI] = refinePoints.second() + 1.0;
            }
            else if( cellLevel[cellI] == minLevel ) // keep from coarsening
            {
                refFld[cellI] = 0.5*(refinePoints.first() + refinePoints.second());
            }
            // else do nothing
        }
        
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    // Enforce refinement at walls
    if( minWallLevel_ > 0 )
    {
        refFld *= 0.0;
        
        forAll(this->boundary(), patchI)
        {
            if( isType<wallFvPatch>(this->boundary()[patchI]) )
            {
                fvPatchScalarField& rfPf = intRefFld.boundaryField()[patchI];
                
                const labelList& pFaceCells = this->boundary()[patchI].faceCells();
                
                forAll(rfPf, pFaceI)
                {
                    label pfCellI = pFaceCells[pFaceI];
                    
                    if( cellLevel[pfCellI] < minWallLevel_ )
                    { //force it to refine
                        refFld[pfCellI] = refinePoints.second() + 1.0;
                    }
                    else if( cellLevel[pfCellI] == minWallLevel_ )
                    { //keep it from coarsening
                        refFld[pfCellI] = 0.5*(refinePoints.first() + refinePoints.second());
                    }
                }
            }
        }
        
        intRefFld.internalField() = max
        (
            intRefFld.internalField(),
            refFld
        );
    }
    
    
    
    intRefFld.correctBoundaryConditions();
    
    Info<<"Min,max refinement field = " << Foam::min(intRefFld).value() << ", "
        << Foam::max(intRefFld).value() << endl;
        
}

void Foam::dynamicRefineBalancedFvMesh::readRefinementDict()
{
    dictionary dynamicMeshDict
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
        )
    );
    
    if( dynamicMeshDict.isDict("refinementControls") )
    {
        dictionary refineControlDict = 
            dynamicMeshDict.subDict("refinementControls");
        
        enableRefinementControl_ = 
            Switch(refineControlDict.lookup("enableRefinementControl"));
        
        if( enableRefinementControl_ )
        {
            // Overwrite field name entry in dynamicRefineFvMeshCoeffs?
            // For now you just have to be smart and enter 
            // 'internalRefinementField' for the field name manually
            
            // Read HashTable of gradient-refinement scalars
            if( refineControlDict.found("gradients") )
            {
                gradFields_ = HashTable< Pair<scalar> >
                (
                    refineControlDict.lookup("gradients")
                );
            }
            
            // Read HashTable of curl-refinement vectors
            if( refineControlDict.found("curls") )
            {
                curlFields_ = HashTable< Pair<scalar> >
                (
                    refineControlDict.lookup("curls")
                );
            }
            
            // Read refinement regions
            if( refineControlDict.found("regions") )
            {
                refinedRegions_ = PtrList<entry>
                (
                    refineControlDict.lookup("regions")
                );
            }
            
            // Read minimum wall level
            if( refineControlDict.found("minWallLevel") )
            {
                minWallLevel_ = readScalar
                (
                    refineControlDict.lookup("minWallLevel")
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
            
Foam::dynamicRefineBalancedFvMesh::dynamicRefineBalancedFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io),
    internalRefinementFieldPtr_(NULL),
    gradFields_(),
    curlFields_(),
    refinedRegions_(),
    minWallLevel_(0),
    enableRefinementControl_(false)
{
    readRefinementDict();
    
    if( enableRefinementControl_ )
    {
        internalRefinementFieldPtr_ = new volScalarField
        (
            IOobject
            (
                "internalRefinementField",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::~dynamicRefineBalancedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedFvMesh::update()
{
    //Part 0 - Update internally calculated refinement field
    readRefinementDict();
    
    if( enableRefinementControl_ )
    {
        updateRefinementField();
    }
    
    //Part 1 - Call normal update from dynamicRefineFvMesh
    bool hasChanged = dynamicRefineFvMesh::update();
    
    // Part 2 - Load Balancing
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
    
    Switch enableBalancing = refineDict.lookup("enableBalancing");
    
    if ( Pstream::parRun() && hasChanged )
    {
        const scalar allowableImbalance = 
            readScalar(refineDict.lookup("allowableImbalance"));
            
        //First determine current level of imbalance - do this for all
        // parallel runs with a changing mesh, even if balancing is disabled
        label nGlobalCells = globalData().nTotalCells();
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/idealNCells;
        
        Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;
        
        //If imbalanced, construct weighted coarse graph (level 0) with node
        // weights equal to their number of subcells. This partitioning works
        // as long as the number of level 0 cells is several times greater than
        // the number of processors.
        if( maxImbalance > allowableImbalance && enableBalancing)
        {
            Info<< "Re-balancing dynamically refined mesh" << endl;
                        
            const labelIOList& cellLevel = meshCutter().cellLevel();
            Map<label> coarseIDmap(100);
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
            
            // Convert to local sequential indexing and calculate coarse
            // points and weights
            labelList localIndex(nCells(),0);
            pointField coarsePoints(nCoarse,vector::zero);
            scalarField coarseWeights(nCoarse,0.0);
            
            forAll(uniqueIndex, cellI)
            {
                localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];
                
                // If 2D refinement (quadtree) is ever implemented, this '3'
                // should be set in general as the number of refinement
                // dimensions.
                label w = (1 << (3*cellLevel[cellI]));
                
                coarseWeights[localIndex[cellI]] += 1.0;
                coarsePoints[localIndex[cellI]] += C()[cellI]/w;
            }
            
            //Set up decomposer - a separate dictionary is used here so
            // you can use a simple partitioning for decomposePar and
            // ptscotch for the rebalancing (or any chosen algorithms)
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
                localIndex,
                coarsePoints,
                coarseWeights
            );

            scalar tolDim = globalMeshData::matchTol_ * bounds().mag();
            
            
            fvMeshDistribute distributor(*this, tolDim);
            
            autoPtr<mapDistributePolyMesh> map =
                  distributor.distribute(finalDecomp);
                  
            meshCutter_.distribute(map);
            
            //Correct values on all cyclic patches
            correctBoundaries<scalar>();
            correctBoundaries<vector>();
            correctBoundaries<sphericalTensor>();
            correctBoundaries<symmTensor>();
            correctBoundaries<tensor>();
        }
    }

    return hasChanged;
}


// ************************************************************************* //
