using StatsBase: sample
using GLMakie  
using Meshes: RectilinearGrid, coordinates, SimpleMesh, vertices, viz!
using StaticArrays
using LinearAlgebra: diagind

using BoundaryWall

begin
N = 20
R = 2.0
σ = -0.25im
T = LinRange(0, 2pi, N)
centers = buildGrid(HexagonalGrid((0.0, 0.0),5.0), 8,8)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)
indices_to_delete = sort(sample(eachindex(centers), replace=false, length(centers)÷4))
# deleteat!(centers, tuple(indices_to_delete...))



circ = [createEllipse(R,R*1/3, T, pi/2, c) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = calcDistances(xm, ym)

f, ax = scatter(xm, ym)
ax.aspect=DataAspect()
f
end


Nx, Ny = 250,150
xdom = LinRange(floor(minimum(x))-2,ceil(maximum(x))+2, Nx)
ydom = LinRange(floor(minimum(y))-2,ceil(maximum(y))+2, Ny)
MESH = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(MESH), MESH.topology)
# COORDS = coordinates.(centroid.(MESH))
COORDS = coordinates.(MESH.vertices)
XDOM, YDOM = first.(COORDS), last.(COORDS)


waveNumber = 1.8
th = 0
waveVector = waveNumber*SVector(cosd(th), sind(th))


GC.gc()

banded = 2

@time wave = boundaryWallWave(waveVector, planeWave, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, banded, Inf);

let
  fig = Figure()
  ax  = Axis(fig[1,1], title="banded (using |i-j|<$banded)")
  viz!(ax, MESH, color=abs2.(wave), colormap=:linear_kbgyw_5_98_c62_n256)
  ax.aspect=DataAspect() 
  # [lines!(ax, GeometryUtils.createCircle(R, LinRange(-pi, pi, 50), c)[1:2]..., color=:black) for ax in ax, c in centers]
  [lines!(ax, createEllipse(R,R*2/3, T, -i*pi/2, c)[1:2]..., color=:white) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  
  fig
end

# 3d image
surface(xdom, ydom, abs2.(reshape(wave, Nx, Ny)), colormap=:amp)

let 
  fig = Figure()
  ax = Axis3(fig[1,1])
  # heatmap!(ax, xdom, ydom, abs2.(reshape(wave, Nx, Ny)), colormap=:amp, interpolate=true)

  band!(ax,Makie.Point3f.(xdom[1], ydom, 0.0), Makie.Point3f.(xdom[1], ydom, abs2.(reshape(wave, Nx, Ny)[1,:])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom[end], ydom, 0.0), Makie.Point3f.(xdom[end], ydom, abs2.(reshape(wave, Nx, Ny)[end,:])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom, ydom[end], 0.0), Makie.Point3f.(xdom, ydom[end], abs2.(reshape(wave, Nx, Ny)[:,end])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom, ydom[1], 0.0), Makie.Point3f.(xdom, ydom[1], abs2.(reshape(wave, Nx, Ny)[:,1])), fxaa=true)
  [lines!(ax, Makie.Point3f.(createEllipse(R,R*2/3, T, -i*pi/2, c)[1:2]...,0.0), color=:black, fxaa=true) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]

  fig
end
