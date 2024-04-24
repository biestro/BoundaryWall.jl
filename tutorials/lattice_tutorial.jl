using StatsBase: sample
using GLMakie  
using Meshes: RectilinearGrid, coordinates, SimpleMesh, vertices, viz!
using StaticArrays
using LinearAlgebra: diagind

include("../src/BWM_geometry.jl")
include("../src/BoundaryWall.jl")

begin
N = 20
R = 2.0
σ = -0.25im
T = LinRange(0, 2pi, N)
centers = GeometryUtils.buildGrid(GeometryUtils.RectangularGrid((0.0, 0.0),5.0, 5.0), 8,8)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)
indices_to_delete = sort(sample(eachindex(centers), replace=false, length(centers)÷4))
# deleteat!(centers, tuple(indices_to_delete...))



circ = [GeometryUtils.createEllipse(R,R*1/3, T, -i*pi/2, c) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)

f, ax = scatter(xm, ym)
ax.aspect=DataAspect()
f
end




Nx, Ny = 250,150
xdom = LinRange(floor(minimum(x))-1,ceil(maximum(x))+1, Nx)
ydom = LinRange(floor(minimum(y))-1,ceil(maximum(y))+1, Ny)
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

@time wave = BoundaryWall.boundaryWallWave(waveVector, BoundaryWall.planeWave, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, banded, Inf);

let
  fig = Figure()
  ax  = Axis(fig[1,1], title="banded (using |i-j|<$banded)")
  viz!(ax, MESH, color=abs2.(wave), colormap=:amp)
  ax.aspect=DataAspect() 
  # [lines!(ax, GeometryUtils.createCircle(R, LinRange(-pi, pi, 50), c)[1:2]..., color=:black) for ax in ax, c in centers]
  [lines!(ax, GeometryUtils.createEllipse(R,R*2/3, T, -i*pi/2, c)[1:2]..., color=:black) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  
  fig
end

# 3d image
surface(xdom, ydom, abs2.(reshape(wave, Nx, Ny)), colormap=:amp)

let 
  fig = Figure()
  ax = Axis3(fig[1,1])
  heatmap!(ax, xdom, ydom, abs2.(reshape(wave, Nx, Ny)), colormap=:amp, interpolate=true)

  band!(ax,Makie.Point3f.(xdom[1], ydom, 0.0), Makie.Point3f.(xdom[1], ydom, abs2.(reshape(wave, Nx, Ny)[1,:])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom[end], ydom, 0.0), Makie.Point3f.(xdom[end], ydom, abs2.(reshape(wave, Nx, Ny)[end,:])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom, ydom[end], 0.0), Makie.Point3f.(xdom, ydom[end], abs2.(reshape(wave, Nx, Ny)[:,end])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom, ydom[1], 0.0), Makie.Point3f.(xdom, ydom[1], abs2.(reshape(wave, Nx, Ny)[:,1])), fxaa=true)
  [lines!(ax, Makie.Point3f.(GeometryUtils.createEllipse(R,R*2/3, T, -i*pi/2, c)[1:2]...,0.0), color=:black, fxaa=true) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]

  fig
end

let 
  GC.gc()
  myred = Makie.to_colormap(:amp)[end÷2]
  fig = Figure(size=(600,500))
  ga = fig[1,1] = GridLayout()
  # gc = fig[2,1] = GridLayout()
  ax,_=heatmap(ga[1,1], xdom, ydom, abs2.(reshape(wave_full, Nx, Ny)), colormap=:amp)

  # surface!(ax, xdom, ydom, abs2.(reshape(wave_full, Nx, Ny)), colormap=:amp)

  abs_wave_in = abs2.(reshape(wave_full,Nx,Ny)[end,:])
  abs_wave_out = abs2.(reshape(wave_full,Nx,Ny)[1,:])
  abs_wave_top = abs2.(reshape(wave_full,Nx,Ny)[:,end])
  abs_wave_bot = abs2.(reshape(wave_full,Nx,Ny)[:,1])
  max_wave = 1.0 # maximum(abs_wave_in)

  hidedecorations!(ax)
  ax_top,_=  lines(ga[0,1], xdom, abs_wave_top,color=:black)
  ax_bot,_=  lines(ga[2,1], xdom, abs_wave_bot  ,color=:black)
  ax_bot.yreversed=true
  ax_left,_= lines(ga[1,0], abs_wave_out ./ max_wave, ydom  ,color=:black)
  ax_left.xreversed=true
  ax_right,_=lines(ga[1,2], abs_wave_in ./ max_wave, ydom,color=:black)

  [ax.xlabel=L"x" for ax in [ax_top, ax_bot]]
  [ax.ylabel=L"y" for ax in [ax_left, ax_right]]
  [ax.ylabelrotation=0 for ax in [ax_left, ax_right, ax_top, ax_bot]]
  [ax.xlabelrotation=0 for ax in [ax_left, ax_right]]
  [ax.xgridvisible=false for ax in [ax_left, ax_right, ax_top, ax_bot]]
  [ax.ygridvisible=false for ax in [ax_left, ax_right, ax_top, ax_bot]]

  ax_right.yaxisposition=:right
  ax_top.xaxisposition=:top
  ax_left.xaxisposition=:top
  ax_top.yaxisposition=:right
  # hidedecorations!(ax_left)
  # hidedecorations!(ax_right)
  # hidedecorations!(ax_top)
  # hidedecorations!(ax_bot)
  ax.aspect=DataAspect()
  
  # hidexdecorations!(ax_left)
  # colgap!(ga, 1, Relative(0.15))
  # rowgap!(ga, 1, Relative(0.05))
  # hidexdecorations!(ax_right)

  [lines!(ax, GeometryUtils.createEllipse(R,R*2/3, T, -i*pi/2, c)[1:2]..., color=:black) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  colgap!(ga, 1, -10.5)
  colgap!(ga, 2, -10.5)
  rowgap!(ga, 1, -10.5)
  rowgap!(ga, 2, -10.5)
  # ax = [ax0, ax1, ax2]
  # [ax.aspect=DataAspect() for ax in ax]

  rowsize!(ga, 1, Aspect(1, 1))
  # colsize!(ga, 1, Aspect(1, 1))
  # colsize!(ga, 1, Aspect(1, 1))
  # rowsize!(ga, 1, Aspect(1, 0.5))
  # band!(ax,Makie.Point3f.(xdom[1], ydom, 0.0), Makie.Point3f.(xdom[1], ydom, abs2.(reshape(wave_full, Nx, Ny)[1,:])), color=(myred, 0.8), fxaa=true)
  # band!(ax,Makie.Point3f.(xdom[end], ydom, 0.0), Makie.Point3f.(xdom[end], ydom, abs2.(reshape(wave_full, Nx, Ny)[end,:])), color=(myred, 0.8), fxaa=true)
  # zlims!(ax, 0, )
  # ax.xgridvisible=false
  # ax.ygridvisible=false
  fig
end
