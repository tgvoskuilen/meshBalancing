meshBalancing
=============

This contains a single library, dynamicRefineBalancedFvMesh which is 
based on dynamicRefineFvMesh but adds mesh balancing for parallel cases
to the update() function.

## OpenFOAM Source Changes

  1. [CRASH] Add guard in refinementHistory::distribute at line 957 to catch
     when fully un-refined cells are transferred
     
        if( newVisibleCells[i] >= 0 )
        {
            visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
        }
        
  2. [ANNOYANCE] Add `if (debug)` guard around print statements in
     `refinementHistory::countProc`
     
  3. [CRASH] Add the following to `nearWallDist::correct()` inside the
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
        
  4. [CRASH] Mapping of patches had an error for mixed patches. In
     `fvMeshAdder::MapVolField` there is a section where imported patch values
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
     
     
