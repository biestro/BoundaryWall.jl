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

Nx = 50
Ny = 50
Nz = 3
xdom = collect(LinRange(-2,2,Nx))
ydom = collect(LinRange(-4,4,Ny))

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coordinates.(vertices(GRID))

XDOM, YDOM = map(i -> getindex.(COORDS,i), 1:2)


N = 150
R  = 1/sqrt(2)
# x,y,xm,ym,rij,ds = BWM.createConfocalBilliard(3.0, 2.0, N)
TH = LinRange(-pi, pi, N)
# circ = [BWM.createEllipse(R, R, TH, 1.0pi, SVector(0.0, y0)) for y0 in [-2.0, 2.0]]
circ = [BWM.createEllipse(0.0, 3R, TH .* 0.5, pi/4, SVector(0.0, 0.0))]
# circ = GeometryUtils.buildConfocalBilliard(3.0, 2.0, N)
x  = vcat(getindex.(circ,1)...)
y  = vcat(getindex.(circ,2)...)
xm = vcat(getindex.(circ,3)...)
ym = vcat(getindex.(circ,4)...)
ds = vcat(getindex.(circ,5)...)
rij = BWM.calcDistances(xm, ym)

scatter(x,y)

fieldEx(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64) = gaussianWave(k, SVector(x,y,z), 5.0; abstol=1e-8)
fieldEy(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64) = polarizedField(k, x, y, z, 0.0)
fieldEz(k::SVector{3, Float64}, x::Float64, y::Float64, z::Float64) = polarizedField(k, x, y, z, 0.0)

# waveNumber = 10.1735# 3.8317
waveNumber = 4sqrt(2.637553)
kRho       = SVector(1.0, 0.0)
kz         = 0.0

waveVector = waveNumber * SVector(kRho..., kz) / norm(SVector(kRho..., kz))



permittivity_strength = -100.0


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
Etot = mapreduce(x -> abs2.(x), +, [Ex, Ey, Ez])
fig,ax,hm=heatmap(xdom, ydom, reshape(Etot, Nx, Ny), colormap=:turbo);scatter!(ax, x,y, color=:white); ax.aspect=DataAspect();display(fig)