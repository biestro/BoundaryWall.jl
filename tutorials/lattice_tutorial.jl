using StatsBase: sample
using GLMakie  
using Meshes: RectilinearGrid, coordinates, SimpleMesh, vertices, viz!
using StaticArrays
using LinearAlgebra: diagind

using BoundaryWall

rand_angle = rand(6*4)*2pi

begin
N = 40
R = 1.5
σ = -0.25im
T = LinRange(0, 2pi, N)
centers = buildGrid(SquareGrid((0.0, 0.0),5.0), 6,4)
# centers = buildGrid(HexagonalGrid((0.0, 0.0),5.0), 8,8)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)
indices_to_delete = sort(sample(eachindex(centers), replace=false, length(centers)÷5))
deleteat!(centers, tuple(indices_to_delete...))
deleteat!(rand_angle, tuple(indices_to_delete...))

# indices_to_delete = [6,7,10,11]
# deleteat!(centers, tuple(indices_to_delete...))



circ = [createEllipse(R,2R/3, T, rand_angle[i], c) for (i,c) in enumerate(centers)]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = calcDistances(xm, ym)

f, ax = scatter(xm, ym)
ax.aspect=DataAspect()
f


Nx, Ny = 150,150
xdom = LinRange(floor(minimum(x))-2,ceil(maximum(x))+2, Nx)
ydom = LinRange(floor(minimum(y))-2,ceil(maximum(y))+2, Ny)
MESH = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(MESH), MESH.topology)
# COORDS = coordinates.(centroid.(MESH))
COORDS = coordinates.(MESH.vertices)
XDOM, YDOM = first.(COORDS), last.(COORDS)


waveNumber = 2.8
th = 0
waveVector = waveNumber*SVector(cosd(th), sind(th))


GC.gc()

banded = 2

@time wave = boundaryWallWave(waveVector, planeWave, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, banded, -10.0);
end
let
  fig = Figure()
  ax  = Axis(fig[1,1], title="Lattice without defects")
  # viz!(ax, MESH, color=abs2.(wave), colormap=:linear_kbgyw_5_98_c62_n256)
  Z = abs2.(wave)
  heatmap!(ax, xdom, ydom, reshape(Z, Nx, Ny), colormap=:linear_kbgyw_5_98_c62_n256)
  ax.aspect=DataAspect() 
  # [lines!(ax, GeometryUtils.createCircle(R, LinRange(-pi, pi, 50), c)[1:2]..., color=:black) for ax in ax, c in centers]
  # [lines!(ax, createEllipse(R,R, T, pi/2, c)[1:2]..., color=:white) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  # [lines!(ax, createEllipse(R,R, T, -i*pi/2, c)[1:2]..., color=:white) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  scatter!(ax, x,y,color=:white, markersize=5)
  save("lattice.png", fig,px_per_unit=4)
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
