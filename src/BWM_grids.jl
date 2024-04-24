# Lattices and stuff

module BWM_grids
using StaticArrays

export ObliqueGrid, RectangularGrid, CenteredRectangular, SquareGrid, HexagonalGrid, HoneyLattice, TriangularGrid, buildGrid


abstract type BravaisLattice end;

struct ObliqueGrid <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
  b::Float64
  ϕ::Float64
end

struct RectangularGrid <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
  b::Float64
end

struct CenteredRectangular <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
  b::Float64
end

struct SquareGrid <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
end

struct HexagonalGrid <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
end

struct HoneyLattice <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
end

struct TriangularGrid <: BravaisLattice
  o::SVector{2, Float64}
  a::Float64
  b::Float64
  s::Float64 # shift
end

"""
  Builds a grid for certain lattice types centered at g.o (origin)
"""
function buildGrid(g::SquareGrid, nx::Int, ny::Int)
  return vec([SVector(g.a * i - g.o[1], g.a * j - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
end

function buildGrid(g::RectangularGrid, nx::Int, ny::Int)
  return vec([SVector(g.a * i - g.o[1], g.b * j - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
end

function buildGrid(g::HexagonalGrid, nx::Int, ny::Int)
  return vec([SVector(g.a * i - g.o[1] + g.a/2 * (j%2), (g.a*j)* sqrt(3)/2 - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
end

function buildGrid(g::CenteredRectangular, nx::Int, ny::Int)
  return vec([SVector(g.a * i - g.o[1] + g.a/2 * (j%2), (g.a*j)* g.b/2 - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
end

function buildGrid(g::TriangularGrid, _n::Int)
  _points = SVector{2, Float64}[]

  for i in 0:_n
    for j in 0:i
      push!(_points, SVector(g.o[1] + i*g.a - g.s*j, g.o[2] + g.b*j)) #sqrt(3)/2*
    end
  end 
  return _points
end

function buildGrid(g::HoneyLattice, _n::Int)
  _centers = buildGrid(TriangularGrid(g.o, g.a, g.a*sqrt(3)/2, g.a/2), _n)
  [push!(_centers, SVector(c[1]+g.a/2, c[2]+g.a/3)) for c in _centers[1:end-_n-1]]

  return _centers
end

end
