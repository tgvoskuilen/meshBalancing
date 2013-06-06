meshBalancing
=============

This contains a single library, dynamicRefineBalancedFvMesh which is 
based on dynamicRefineFvMesh but adds mesh balancing for parallel cases
to the update() function.

Note: The redistributeParPlus function does NOT work. For some reason it does not
redistribute the cellLevel and pointLevel fields (hexRef8.distribute crashes) so
it is only redistributing the refinement history. The cellLevel volScalarField is
mapped, but pointLevel is a bit more tricky.

A workaround may be to reconstruct the cellLevel and pointLevel manually from the
distributed refinement history.

## OpenFOAM Source Changes (for 2.1.x)

  1. [ __CRASH__ ] Add guard in src/dynamicMesh/polyTopoChange/refinementHistory.C in
     refinementHistory::distribute at line 927 to catch
     when fully un-refined cells are transferred
     
        if( newVisibleCells[i] >= 0 )
        {
            visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
        }
        
  2. [ __ANNOYANCE__ ] Add `if (debug)` guard around print statements in
     `refinementHistory::countProc` in the same file as edit #1 to prevent excessive
     printing to stdout during mesh balancing operations
     
  3. [ __CRASH__ ] Add the following to src/finiteVolume/fvMesh/nearWallDist.C
     at `nearWallDist::correct()` inside the
     `mesh_.changing()` condition but before the patch sizes are set. If the
     total number of patches on a processor increases, there is a crash unless
     we increase the size of nearWallDist first.
          
        // If the number of patches on this processor increased, we need to
        // increase the number of patches in nearWallDist to match
        if( mesh_.boundary().size() != size() )
        {
            //Resize nearWallDist to match mesh patches
            resize(mesh_.boundary().size());

            //Set any newly created (unset) patches
            forAll(*this, patchI)
            {
                if( !set(patchI) )
                {
                    set
                    (
                        patchI, 
                        fvPatchField<scalar>::New
                        (
                            calculatedFvPatchScalarField::typeName,
                            mesh_.boundary()[patchI],
                            mesh_.V()
                        )
                    );
                }
            }
        }
        
  4. [ __CRASH__ ] Mapping of patches had an error for mixed patches. In
     src/dynamicMesh/fvMeshAdder/fvMeshAdderTemplates.C in function
     `fvMeshAdder::MapVolField` around line 263 there is a section where imported patch values
     are mapped to existing patches. This correctly maps the patch value, but
     omits mapping of the refValue, refGradient, and other stored fields in 
     mixed patches. The refValue is then left
     set to a section of allocated but unset memory if all the original faces 
     of that patch get moved to another processor. Remove the manual loop
     near the end of the function and use the rmap function with a reversed map
     by replacing
     
        forAll(newFld, i)
        {
            label oldFaceI = newToAdded[i];
        
            if (oldFaceI >= 0 && oldFaceI < addedFld.size())
            {
                newFld[i] = addedFld[oldFaceI];
            } 
        }
        
    with
    
        labelList addedToNew(addedFld.size(),-1);
        forAll(newFld, i)
        {
            label oldFaceI = newToAdded[i];

            if (oldFaceI >= 0 && oldFaceI < addedFld.size())
            {
                addedToNew[oldFaceI] = i;
            } 
        }
        
        newFld.rmap(addedFld, addedToNew);
     
  5.  [ __CRASH__ ] When using a chemistry solver, the DimensionedField `deltaTChem` in 
      basicChemistryModel.H is not properly mapped/distributed. This is because both
      `fvMeshDistributor::distribute` and `fvMeshAdder::add` only consider GeometricFields,
      ignoring DimensionedFields. However, `fvMesh::mapFields` looks at both Geometric
      and Dimensioned fields.

      An easy workaround is to change `deltaTChem` to a GeometricField by making the following
      edits, then recompiling the entire src directory since a lot of libraries link to this.
      
      1. In basicChemistryModel.H change `DimensionedField<scalar, volMesh> deltaTChem_;`
         to `GeometricField<scalar, fvPatchField, volMesh> deltaTChem_;` and add
         `#include GeometricField.H` and `#include fvPatchField.H` to the list of included
         headers

      2. In basicChemistryModelI.H change both instances of `return deltaTChem_;` to
         `return deltaTChem_.dimensionedInternalField();`
      
      

