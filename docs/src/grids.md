# Lattices & Grids

```@example grids
using StaticArrays # hide
using BoundaryWall # hide
using CairoMakie # hide
```

These methods implement grids similar to a Bravais lattice, as well as 
their repetition scheme.

## Constructors

<!--```@docs-->
<!--RectangularGrid-->
<!--```-->

<!--```@docs-->
<!--SquareGrid-->
<!--```-->

<!--```@docs-->
<!--HexagonalGrid-->
<!--```-->

<!--```@docs-->
<!--TriangularGrid-->
<!--```-->

<!--```@docs-->
<!--HoneyLattice-->
<!--```-->

## Methods

All `buildGrid` functions construct a grid with $n\\_1\times n\\_2$ 
elements (depending on their basis vectors).

```@docs
buildGrid
```

```@example grids
buildGrid(HexagonalGrid(SVector(0.0, 0.0), 0.5), 8, 6) |> scatter
```
