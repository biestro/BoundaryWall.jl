# Beam splitters

Using the [beam shaping](incident.md) formalism, this tutorial impinges a Gaussian
beam against a negative $\delta$ potential (note that it is only by setting a negative 
potential strength that one correctly models the optical behaviour).

```@example beam
using WGLMakie
using StaticArrays: SVector
using LinearAlgebra
using Meshes
using BoundaryWall

begin # definitions
HBAR        = 1.0
MASS        = HBAR/2
SIGMA       = (2*MASS/HBAR^2)*(1/4*im)
N           = 100
NX          = 150
NY          = 150
ϕ           = 180.0
waveVector  = 5*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
end

centers = ((LinRange(0.0, 1.0, 4)).^0.5) * 3.0 # barrier = [createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), th, SVector(x0, 0.0)) for (x0, th) in zip([0.0, 2.0], [pi/6, pi/4])]
barrier = [createEllipse(0.0, 8.0, LinRange(-pi/2,pi/2,N), pi/4, SVector(x0, 0.0)) for x0 in centers[1]]
x  = vcat(getindex.(barrier, 1)...)
y  = vcat(getindex.(barrier, 2)...)
xm = vcat(getindex.(barrier, 3)...)
ym = vcat(getindex.(barrier, 4)...)
ds = vcat(getindex.(barrier, 5)...)
rij = calcDistances(xm, ym)

# domain
x0, xf = (-8, 8)
y0, yf = (-8, 8)
xdom = LinRange(x0, xf, NX)
ydom = LinRange(y0, yf, NY)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coords.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)

banded = 5

@time wave = boundaryWallWave(waveVector, (k,r)->gaussianWave(k,r - SVector(0.0, 0.0), 10.0; abstol=1e-8), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, -5.0);
```

```@example beam
function plotBeam(_xdom, _ydom, _wave, _N) # hide
  fig = Figure(size=(600,600)) # hide
  ax = Axis(fig[1,1],xtickalign=1,ytickalign=1,xticksmirrored=1,yticksmirrored=1,xgridvisible=false, ygridvisible=false,xminorticksvisible=true,yminorticksvisible=true,xminortickalign=1,yminortickalign=1,xminorticks=IntervalsBetween(2),yminorticks=IntervalsBetween(2),) # hide
  heatmap!(ax,_xdom, _ydom, abs2.(reshape(_wave, length(_xdom), length(_ydom))), interpolate=true, colormap=:linear_kbgyw_5_98_c62_n256) # hide

  [lines!(ax, createEllipse(0.0, 7.0, LinRange(-pi/2,pi/2,_N), pi/4, SVector(x0, 0.0))[1:2]...,color=(:white,1.0)) for x0 in centers[1]] # hide
  xlims!(ax, _xdom[1],_xdom[end]) # hide
  ylims!(ax, _ydom[1],_ydom[end]) # hide
  ax.aspect=DataAspect() # hide
  return fig # hide
end # hide
f = plotBeam(xdom, ydom, wave, N)
```