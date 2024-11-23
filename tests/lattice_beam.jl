using BoundaryWall
using GLMakie
using StatsBase
using LinearAlgebra
using Meshes
using StaticArrays

N = 50 
NX, NY = 100,100
SIGMA = -0.25im
T = LinRange(-pi, pi, N)
R = 0.5
g = SquareGrid(SVector(0.0, 0.0), 3.5R)
centers = buildGrid(g, 6, 6)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)
indices_to_delete = sort(sample(eachindex(centers), replace=false, length(centers)÷2))
# deleteat!(centers, tuple(indices_to_delete...))

circ = [createCircle(R, T, c) for c in centers]
# circ = [createEllipse(0.0, R, LinRange(-pi/2,pi/2,N),0.0, c) for c in centers]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...) 
rij = calcDistances(xm, ym)

# domain
x0, xf = (-7, 7)
y0, yf = (-7, 7)
xdom = LinRange(x0, xf, NX)
ydom = LinRange(y0, yf, NY)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coords.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)

banded=2
# @time wave = boundaryWallWave(SVector(-5.0,0.0), planeWave, x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, 10.0);
@time wave = boundaryWallWave(SVector(-5.0,0.0), (k,r)->gaussianWave(k,r - SVector(-5.0, 0.0), 25.0; abstol=1e-9), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, Inf);

fig = Figure()
ax = Axis(fig[1,1])
ax.aspect=DataAspect()
xlims!(ax, x0,xf)
ylims!(ax, y0,yf)
dens = Observable(abs2.(reshape(wave, NX,NY)))
heatmap!(ax, xdom, ydom, dens, colormap=:balance, interpolate=true)
scatter!(ax, x,y, color=:black)
lines(fig[1,2], abs2.(reshape(wave, NX,NY)[end,:]), ydom)
rowsize!(fig.layout, 1,Aspect(1,1))
fig
