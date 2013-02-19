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
