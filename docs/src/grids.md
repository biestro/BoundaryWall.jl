# Lattices & Grids

```@example grids
using StaticArrays # hide
using BoundaryWall # hide
using CairoMakie # hide
```

These methods implement grids similar to a Bravais lattice, as well as 
their repetition scheme.

## Constructors

## Methods

All `buildGrid` functions construct a grid with $n\\_1\times n\\_2$ 
elements (depending on their basis vectors).

```@docs
buildGrid
```

```@example grids
buildGrid(HexagonalGrid(SVector(0.0, 0.0), 0.5), 8, 6) |> scatter
```

One can use these as centers for other boundaries.

```@example grids
N       = 20
R       = 0.2
centers = buildGrid(HexagonalGrid(SVector(0.0, 0.0), 0.5), 8, 6)
circ    = [createEllipse(R,2R/3, LinRange(-pi,pi,N), rand()*2pi, c) for c in centers]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
scatter(x,y)
```