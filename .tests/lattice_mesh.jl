using HOHQMesh

using GLMakie



include("../src/GeometryUtils.jl")

begin
nx, ny  = (5,5)
origin  = (0.0, 0.0)
n_holes = nx*ny ÷ 3
grid    = GeometryUtils.SquareGrid(origin, 2.0)
R       = grid.a / 4.5
N       = 90
TH      = LinRange(-pi/2,pi/2,N)

centers = GeometryUtils.buildGrid(grid, nx, ny)
centers_3d = map(x -> [x[1], x[2], 0.0], centers)
circ_points = [GeometryUtils.createEllipse(R, R/2, TH, 0.0, c) for c in centers]
x  = vcat(getindex.(circ_points, 1)...)
y  = vcat(getindex.(circ_points, 2)...)
xm = vcat(getindex.(circ_points, 3)...)
ym = vcat(getindex.(circ_points, 4)...)
ds = vcat(getindex.(circ_points, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)

end

begin
scattering_project = newProject("ellipse", "out")
setMeshFileFormat!(scattering_project, "ABAQUS")

ellipse_1  = [LinRange(0.0,1.0,N+1) x[1:N+1]  y[1:N+1] zeros(N+1)]
ellipse_2  = [LinRange(0.0,1.0,N+1) -x[1:N+1]  -y[1:N+1] zeros(N+1)]

ellipse_spline = newSplineCurve("left",  N, ellipse_1)
ellipse_spline_2 = newSplineCurve("left",  N, ellipse_2)

setPolynomialOrder!(scattering_project,10)
addCurveToInnerBoundary!(scattering_project, ellipse_spline, "el1")
addCurveToInnerBoundary!(scattering_project, ellipse_spline_2, "el2")

xf = maximum(x[1:N]) + 2R ; x0 = minimum(x[1:N]) - 2R; 
y0 = minimum(y[1:N]) - 2R; yf = maximum(y[1:N]) + 2R
n_grid = (100,90)
bounds = [yf, x0, y0, xf]
N_grid = [n_grid[1], n_grid[2], 0]
addBackgroundGrid!(scattering_project, bounds, N_grid)

generate_mesh(scattering_project)

plotProject!(scattering_project,5)

end

begin
scattering_project = newProject("lattice", "out")
setMeshFileFormat!(scattering_project, "ABAQUS")

circ = newCircularArcCurve.("rod",   
                             centers_3d, 
                             R, 
                             0.0, 360.0) 



xf = maximum(x) + 2R ; x0 = minimum(x) - 2R; 
y0 = minimum(y) - 2R; yf = maximum(y) + 2R
n_grid = (100,90)
bounds = [yf, x0, y0, xf]
N_grid = [n_grid[1], n_grid[2], 0]
addBackgroundGrid!(scattering_project, bounds, N_grid)
[addCurveToInnerBoundary!(scattering_project, c, "circle$i") for (i,c) in enumerate(circ)]
generate_mesh(scattering_project)

plotProject!(scattering_project, 5)


end

using Meshes
_points, _connec = GeometryUtils.abaqusToMeshes("./out/lattice.inp")

MESH = Meshes.SimpleMesh(_points, _connec)
COORDS = coords.(MESH.vertices)
XDOM, YDOM = first.(COORDS), last.(COORDS)

# boundary wall
using StaticArrays
include("../src/BoundaryWall.jl")
waveVector = SVector(5.0, 0.0)


wave = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, -0.25im, ds, rij, length(ds), N, 4, Inf)


let 
  fig = Figure()
  ax = Axis(fig[1,1])
  viz!(ax, MESH, showsegments=false, color=abs2.(wave), colormap=:turbo)
  ax.aspect=DataAspect()
  fig
end