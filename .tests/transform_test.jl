using GLMakie, Meshes
using CoordRefSystems: Polar

grid = CartesianGrid(10,10)

mesh = grid |> StdCoords()
Proj(Polar)

viz(mesh |> Meshes.proj(Polar), showsegments=true)