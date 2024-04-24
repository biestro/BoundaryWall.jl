using GLMakie
using StaticArrays
using Meshes
using QuadGK
using LinearAlgebra

include("../src/BoundaryWall.jl")
include("../src/GeometryUtils.jl")

N = 250
# circ = [GeometryUtils.createEllipse(-2.5, 7.0, LinRange(pi/2, 3pi/2, N), 0.0, SVector(x0,0.0)) for x0 in [0.0, -0.5]]
circ = [GeometryUtils.createEllipse(-1.0, 1.5, LinRange(-pi, pi, N), 0.0, SVector(0.0,0.0))]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = GeometryUtils.calcDistances(xm,ym)

lines(x,y)

Nx, Ny = 100,100
xdom = LinRange(-1.5,1.5, Nx)
ydom = LinRange(-2.0,2.0, Ny)
MESH = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(MESH), MESH.topology)
# COORDS = coordinates.(centroid.(MESH))
COORDS = SVector.(coordinates.(MESH.vertices))
XDOM, YDOM = first.(COORDS), last.(COORDS)

planeIncidence(k::SVector, x::Float64, y::Float64) = exp(- 1im  *  (k[1] * x + k[2] * y ));
planeIncidence(k::SVector, x::T, y::T) where T <:SArray{Tuple{}, Float64, 0, 1} = exp.(- 1im  *  (k[1] * x + k[2] * y ));


shapedKernel(k::Float64, r::SVector{2, Float64}, θ::Float64,β::Float64, ω::Float64) = exp( -((θ-β)/ω) ^2) * planeIncidence(k*SVector(cos(θ), sin(θ)), r[1], r[2])
gaussianKernel(k::Float64, r::SVector{2, Float64}, θ::Float64, β::Float64, ω::Float64) = pi * ω^2 * exp(-(ω*(θ-β)/2)^2) * planeIncidence(k*SVector(cos(θ), sin(θ)), r[1], r[2])

using HCubature
"""
from the right
"""
# shapedIncidence(k::Float64, r::SVector{2,Float64}, r0::SVector{2,Float64}, β::Float64, ω::Float64) = 1/2π * first(quadgk(t -> shapedKernel(k, r - r0, t, β , ω), -π/2+β, π/2+β))
shapedIncidence(k::Float64, r::SVector{2,Float64}, r0::SVector{2,Float64}, β::Float64, ω::Float64) = 1/2π * first(hquadrature(t -> shapedKernel(k, r - r0, t, β , ω), -π/2+β, π/2+β))
gaussianIncidence(k::Float64, r::SVector{2,Float64}, r0::SVector{2,Float64}, β::Float64, ω::Float64) = 1/2π * first(hquadrature(t -> gaussianKernel(k, r - r0, t, β , ω), -π/2+β, π/2+β))

# shapedIncidence(20.0, COORDS, SVector(-1.0, 0.0), pi/3, 1.0)
using BenchmarkTools

waveVec = 5.1735*SVector(cosd(45), sind(45))
incidentWave(k::SVector{2, Float64}, x::Float64, y::Float64) = shapedIncidence(norm(k), SVector(x,y), SVector(0.0, 0.0), pi/2 - atan(k...), 0.5)
gaussianWave(k::SVector{2, Float64}, x::Float64, y::Float64) = gaussianIncidence(norm(k), SVector(x,y), SVector(0.0, 0.0), pi/2-atan(k...), 10.14)


incident = [gaussianWave(waveVec, r[1], r[2]) for r in COORDS]
let
GC.gc();
f,ax=viz(MESH, color=real(incident))
lines!(ax,x,y)
xlims!(ax, minimum(XDOM), maximum(XDOM))
ylims!(ax, minimum(YDOM), maximum(YDOM))
ax.aspect=DataAspect()
f
end


banded = 10
wave = BoundaryWall.boundaryWallWave(
  waveVec,
  incidentWave,
  x,
  y,
  xm,
  ym,
  XDOM, 
  YDOM,
  -0.25im,
  ds,
  rij,
  length(ds),
  N,
  banded,
  6.0
)

function plotWave(_m::SimpleMesh,_wave::Union{Vector{ComplexF64}, Vector{Float64}}, _c::Symbol)
    _fig = Figure()    
    _ax  = Axis(_fig[1,1], xticksmirrored=true, yticksmirrored=true, 
    xminorticksvisible=true, yminorticksvisible=true,
    # xtickalign=1, ytickalign=1,
    # xminortickalign=1, yminortickalign=1,
    xgridvisible=false,ygridvisible=false)
    viz!(_ax, _m; color=_wave,colormap=_c)
    Colorbar(_fig[1,2], colorrange=(minimum(_wave), maximum(_wave)), colormap=_c)
    _ax.aspect=DataAspect()
    return _fig

end


let
  
f=plotWave(MESH, abs2.(wave), :turbo)

# scatter!(x,y, color=(:white,1.0), markersize=5)
# [lines!(x[n:n+N], y[n:n+N], color=:black) for n in [1,N+2]]
lines!(x,y,color=:white)
xlims!(minimum(XDOM), maximum(XDOM))
ylims!(minimum(YDOM), maximum(YDOM))
f
end

lines(x[1:N+2], y[1:N+2])