#!/usr/bin/python

import os
import sys
import math

toolPath = '../../../foamTools/python'

sys.path.append(toolPath)

import pyOpenFOAM
from meshBuilder import *

# Set the mesh base unit to millimeters (conversion factor to meters)
baseUnit = 0.001

# target cell size (mm)
cellSize = 2./3.

# Core dimensions
Li = 4  # mm, inlet length
Wi = 4   # mm, inlet width
Lo = 80  # mm, outlet length
Wo = 40  # mm, outlet width



# Build the mesh
sqIn = Block((-Li,-Wi/2.,-Wi/2.),(0,Wi/2.,Wi/2.),(1,1,1))

sq1 = Block((0, Wi/2.,Wi/2.),(Lo,Wo/2.,Wo/2.),(1,1,1))
sq2 = Block((0,-Wi/2.,Wi/2.),(Lo,Wi/2.,Wo/2.),(1,1,1))
sq3 = Block((0,-Wo/2.,Wi/2.),(Lo,-Wi/2.,Wo/2.),(1,1,1))

sq4 = Block((0, Wi/2.,-Wi/2.),(Lo, Wo/2.,Wi/2.),(1,1,1))
sq5 = Block((0,-Wi/2.,-Wi/2.),(Lo, Wi/2.,Wi/2.),(1,1,1))
sq6 = Block((0,-Wo/2.,-Wi/2.),(Lo,-Wi/2.,Wi/2.),(1,1,1))

sq7 = Block((0, Wi/2.,-Wo/2.),(Lo, Wo/2.,-Wi/2.),(1,1,1))
sq8 = Block((0,-Wi/2.,-Wo/2.),(Lo, Wi/2.,-Wi/2.),(1,1,1))
sq9 = Block((0,-Wo/2.,-Wo/2.),(Lo,-Wi/2.,-Wi/2.),(1,1,1))

mesh = BlockMesh([sqIn, sq1, sq2, sq3, sq4, sq5, sq6, sq7, sq8, sq9],"patch","air")
    

# tag the inlet and outlet faces based on a plane point and normal (all untagged faces are default)
mesh.tagFaces((-Li,0,0),(-1,0,0),"fuel",patchType="patch")

mesh.tagFaces((-Li/2.,Wi/2,0),(0,1,0),"walls",patchType="wall")
mesh.tagFaces((-Li/2.,-Wi/2,0),(0,-1,0),"walls",patchType="wall")
mesh.tagFaces((-Li/2.,0,Wi/2),(0,0,1),"walls",patchType="wall")
mesh.tagFaces((-Li/2.,0,-Wi/2),(0,0,-1),"walls",patchType="wall")

# set the mesh scale in cells/arb unit
mesh.scale = 1. / cellSize

# Set the mesh base unit
mesh.baseUnit = baseUnit

f = open('constant/polyMesh/blockMeshDict','w')
f.write(str(mesh))
f.close()
    

