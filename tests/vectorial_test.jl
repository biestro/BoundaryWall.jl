using GLMakie, StaticArrays, LinearAlgebra
using Meshes
using StatsBase: sample
using Distributions: Uniform
using CircularArrays
using QuadGK
using ThreadsX
using Meshes


include("../src/GeometryUtils.jl")
include("../src/BoundaryWall.jl")

N = 400
R = 0.8
σ = -0.25im
T = LinRange(-pi, pi, N+1)
circ = [GeometryUtils.createCircle(R, T, (0.0, i)) for i in [-2.0, 0.0, 2.0]]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)

rMid = SVector.(xm, ym)
rPos = CircularArray(SVector.(x,y)[1:end-1])

waveNumber = 7.0156
th = pi/2
waveVector = SVector(cosd(th), sind(th)) * waveNumber
Mxx, Mxy, Myy, Mzz = BoundaryWall.calcGreenTensor(waveNumber, 
                             deg2rad(th), 
                             σ, 
                             rij,
                             ds, 
                             4, 
                             rMid, 
                             rPos,
                             N)
# interval (0.5, 0.5000000000001137) 
# interval (0.49999999999994316, 0.5)
# interval (0.5, 0.5000000000000284) -> error at 0.5000000000000142
# interval (0.4999999999999716, 0.5) -> error at 0.4999999999999858:
Myy[diagind(size(Myy)...)] .= 0
lines(abs2.(diag(Myy,3)))
heatmap(abs2.(Myy), colorscale=log10)


function polarizedField(k::SVector{2, Float64}, x::Float64, y::Float64, ϕ::Float64)
  return exp(im*ϕ)*exp(-im * (k[1] * x + k[2] * y))
end


Ex0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarizedField(k, x, y, 0.0)
Ey0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarizedField(k, x, y, 0.0)
Ez0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarizedField(k, x, y, 0.0)

Nx = 100
Ny = 100
xdom = collect(LinRange(-5,5,Nx))
ydom = collect(LinRange(-5,5,Ny))

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coords.(vertices(GRID))

XDOM, YDOM = first.(COORDS), last.(COORDS)

GC.gc()
Ex, Ey, Ez = BoundaryWall.boundaryWallVec(waveVector, xm, ym, x, y, XDOM, YDOM, σ, rij, ds, length(ds), N, SVector(Ex0, Ey0, Ez0))

Etot = mapreduce(x->x.^2, +, [Ex, Ey, Ez])

let
f,ax=heatmap(xdom, ydom, abs2.(reshape(Etot, Nx, Ny)), interpolate=false,colormap=:turbo)
lines!(ax, x, y, color=:black)
ax.aspect=DataAspect()
f
end