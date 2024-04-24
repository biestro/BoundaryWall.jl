"""
Computational implementation of the Boundary Wall Method by M.G.E. da 
Luz, et al. (1997). Quantum scattering from arbitrary boundaries. Phys.
Rev. E. APS.

Solve scattering from boundaries and lattices

Copyright © 2023 - 2024 Alberto Ruiz Biestro. All rights reserved. (Will update license later on.)
"""

module BoundaryWallMethod

# core
using LinearAlgebra, StaticArrays

include("BWM_geometry.jl")
include("BWM_main.jl")
include("BWM_grids.jl")

# functions from BoundaryWall
#
export boundaryWallWave, boundaryWallVec, planeWave, gaussianWave, shapedWave, gradient

# from BWM_grids.jl
export ObliqueGrid, RectangularGrid, CenteredRectangular, SquareGrid, HexagonalGrid, HoneyLattice, TriangularGrid, buildGrid

# from BWM_geometry.jl
export createCircle, createEllipse, createConfocalBilliard, calcArcLength, calcDistances, calcMidpoints, divideCurve, divideResonatorCurve

end
