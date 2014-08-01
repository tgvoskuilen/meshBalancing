meshBalancing
=============

This contains a single library, dynamicRefineBalancedFvMesh which is 
based on dynamicRefineFvMesh but adds mesh balancing for parallel cases
to the update() function.

## Usage

To use the balancing library with an existing 'DyM' solver you have to make the following
edits to the case (in addition to the changes to the source code in the section below)

  1. Add a linking statement to the new library in controlDict
  
        libs
        (
            "libdynamicFvMesh-dev.so"
        );

  2. Add the balanceParDict dictionary to the system folder (example in the library folder)
  3. Change the mesh type in dynamicMeshDict to `dynamicRefineBalancedFvMesh`
  4. Add the following two entries to the dynamicRefineFvMeshCoeffs dictionary (see the example
     in the library folder)
  
        enableBalancing true;
        allowableImbalance 0.15;

  5. You can also add a `refinementControl` entry to enable refinement based on
     field gradients, field curl, and specified regions of refinement.
  
## Notes

To use this with interDyMFoam, you have to move the call to createPrghCorrTypes
inside correctPhi.H to avoid a crash when the number of patches changes and
recompile the solver.

## OpenFOAM Source Changes (for 2.2.x)

To use this library with interDyMFoam and similar solvers, you do not need to
make all the listed edits. Fixing issues 1, 2, 3, 4, and 6 should allow most
simulations to be run without error. If you are using non-Newtonian viscosity
models you will have to fix issue 7 too.

  1. [ __CRASH__ ] Add guard in src/dynamicMesh/polyTopoChange/polyTopoChange/refinementHistory.C in
     refinementHistory::distribute at line 927 to catch
     when fully un-refined cells are transferred
     
        if( newVisibleCells[i] >= 0 )
        {
            visibleCells_[constructMap[i]] = newVisibleCells[i] + offset;
        }
        
  2. [ __ANNOYANCE__ ] Add `if (debug)` guard around print statements in
     `refinementHistory::countProc` in the same file as edit #1 to prevent excessive
     printing to stdout during mesh balancing operations
     
  3. [ __CRASH__ ] Add the following to src/finiteVolume/fvMesh/wallDist/nearWallDist.C
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
        
  4. [ __WARNINGS__ ] There was a small change in version 2.2.x in
     src/dynamicMesh/fvMeshAdder/fvMeshAdderTemplates.C which results in an
     excessive number of warnings. In the `MapVolField` function, there are
     two calls to `calcPatchMap`. The map given to the patches may contain
     unmapped cells until the last mapping. To supress the warnings, change the
     unmapped value entry from -1 to 0 (the value it used to be) on lines 139
     and 201. A similar edit is required in `MapSurfaceField` on lines 447
     and 508.
     
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
      
      The more complex workaround is to add distribution of DimensionedFields to
      the appropriate functions. See the 2.3.x branch for details on this.

  6.  [ __METHOD ERROR__ ] In the current implementation of dynamicRefineFvMesh, on which
      this is based, a refined cell can be coarsened if the minimum refinementField value
      in any of its child cells is less than the threshhold in the dictionary. This is
      nonsense and leads to oscillatory refinement. Consider the following cell, showing the
      value of `refinementField` in each cell
      
         +------+------+
         |      |      |
         |  0   | 9.9  |
         |      |      |
         +------+------+
         |      |      |
         | 9.9  | 9.9  |
         |      |      |
         +------+------+
         
     where refinement is triggered at a value of 10 and unrefinement at a value of, say, 0.5. 
     This cell, rather than being left
     alone, will be coarsened. It should use the maximum value rather than the minimum
     value. To change this, in `dynamicRefineFvMesh.H` on line 112 change `minCellField`
     to `maxCellField`. In `dynamicRefineFvMesh.C`, make the following 4 changes. Line 
     numbers here may not match your line numbers exactly, so use some common sense.
     Changes 1-3 will all be in the existing `minCellField` function and change 4 is the only
     location in the file that referenced the `minCellField` function.
     
     1. Line 651: change `minCellField` to `maxCellField`
     2. Line 653: change `GREAT` to `-GREAT`
     3. Line 661: change `min` to `max`
     4. Line 1342: change `minCellField` to `maxCellField`
     
  7. [ __CRASH__ ] When using viscosity models other than Newtonian in multiphase systems, each
     model creates a field named "nu" which conflict with each other when re-balancing the mesh.
     To fix it, for example in `BirdCarreau.C`, change `"nu"` on line 79 to `"BirdCarreauNu."+name`
     so the field has a unique name. A similar modification can be made for the other viscosity
     models if needed.
