# TODO: broadcast the calculation of greens tensor (i.e. remove the if else)

using LinearAlgebra
using StaticArrays: SVector
using SpecialFunctions
using GLMakie
using CircularArrays

import BoundaryWall as BWM

using Meshes

function polarizedField(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64, ϕ::Float64)
  # φ: arbitrary phase
  return exp(im*ϕ)*exp(-im * (k[1] * x + k[2] * y + k[3] * z))
end

Nx = 150
Ny = 150
Nz = 3
xdom = collect(LinRange(-25,25,Nx))
ydom = collect(LinRange(-20,20,Ny))

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coords.(vertices(GRID))

XDOM, YDOM = map(i -> getindex.(COORDS,i), 1:2)

begin
  N = 151
  R = 1.0
  σ = -0.25im
  T = LinRange(0, 2pi, N)
  centers = BWM.buildGrid(BWM.SquareGrid((0.0, 0.0),5.0), 3,3)
  Δr = maximum(centers) ./ 2
  map!(c->c .- Δr, centers, centers)
  indices_to_delete = sort(sample(eachindex(centers), replace=false, length(centers)÷4))
  deleteat!(centers, (5))
  
  
  
  circ = [BWM.createEllipse(R,R, T, pi/2, c) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  x  = vcat(getindex.(circ, 1)...)
  y  = vcat(getindex.(circ, 2)...)
  xm = vcat(getindex.(circ, 3)...)
  ym = vcat(getindex.(circ, 4)...)
  ds = vcat(getindex.(circ, 5)...)
  rij = BWM.calcDistances(xm, ym)
  
  f, ax = scatter(xm, ym)
  ax.aspect=DataAspect()
  f
end 

fieldEx(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64) = polarizedField(k, x, y, z, 0.0)#gaussianWave(k, SVector(x,y,z), 5.0; abstol=1e-8)
fieldEy(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64) = polarizedField(k, x, y, z, 0.0)
fieldEz(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64) = 0.0polarizedField(k, x, y, z, 0.0)

# waveNumber = 10.1735# 3.8317
waveNumber = 1.5
kRho       = SVector(0.5, 0.0)
kz         = 1.0

waveVector = waveNumber * SVector(kRho..., kz) / norm(SVector(kRho..., kz))



permittivity_strength = -2.0


Ex, Ey, Ez = BWM.boundaryWallVec(waveVector, 
                        xm,ym,
                        x,y,
                        XDOM, 
                        YDOM,
                        0.0, 
                        -0.25im, 
                        rij,
                        ds, 
                        length(ds), 
                        N,
                        SVector(fieldEx,fieldEy,fieldEz),
                        permittivity_strength)

# plotting
# fig,ax,hm=contour(xdom, ydom, reshape(abs2.(Ez), Nx, Ny), levels=LinRange(0,3,10), colormap=:turbo);lines!(ax, x,y, color=:white); ax.aspect=DataAspect();display(fig)
# Etot = mapreduce(x -> abs2.(x), +, [Ey])

Etot = mapreduce(x -> abs2.(x), +, [Ex, Ey,Ez])
fig,ax,hm=heatmap(xdom, ydom, reshape(abs2.(Ez), Nx, Ny), colormap=:turbo);scatter!(ax, x,y, color=:white); ax.aspect=DataAspect();display(fig)

fig, ax = contour(xdom, ydom, reshape(Etot, Nx, Ny), levels=LinRange(0,20,10), colormap=:turbo); scatter!(ax, x,y,color=:black);fig