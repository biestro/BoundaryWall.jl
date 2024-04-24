using StatsBase: sample
using GLMakie  
using Meshes: RectilinearGrid, coordinates, SimpleMesh, vertices, viz!, viz
using StaticArrays
using LinearAlgebra: diagind

include("../src/GeometryUtils.jl")
include("../src/BoundaryWall.jl")





begin
N = 20
n_grid = 10
R = 1.0/n_grid/2
σ = -0.25im
T = LinRange(0, 2pi, N)
# centers = GeometryUtils.buildGrid(GeometryUtils.HoneyLattice((0.0, 0.0), 5.0), 8)
centers = GeometryUtils.buildGrid(GeometryUtils.TriangularGrid((0.0, 0.0), 4R, 4R*sqrt(3)/2, 2R), n_grid)

Δx = maximum(first.(centers)) ./ 2
Δy = 0.0
map!(c->c .- (Δx,Δy), centers, centers)
# centers = SVector.(first.(centers), last.(centers))


circ = [GeometryUtils.createEllipse(R,R, T,0.0, c) for (c,i) in zip(centers, LinRange(0,1, length(centers)))]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)

f, ax = scatter(xm, ym)
source_point = coordinates(centroid(Triangle(Tuple.([centers[end-n_grid-1], centers[end], centers[end-1]])...)))
# f, ax = scatter(centers)
scatter!(ax,source_point)
# scatter!(ax, )
ax.aspect=DataAspect()
f
end

Nx, Ny = 100,100

using Meshes


begin
  GC.gc()

refine_iterative(x::SimpleMesh,n::Int, method::RefinementMethod) = foldl(|>, fill(i->refine(i, method),n), init=x)

# foldl

tri = discretize(PolyArea((minimum(x)-3R, minimum(y)-R), 
                          (0.0, maximum(y)+3R), 
                          (maximum(x)+3R,minimum(y)-R)), 
                          FIST())
mesh = refine_iterative(tri, 7, TriSubdivision())
mesh = refine(mesh, TriRefinement())
dom = coordinates.(mesh.vertices)
XDOM = first.(dom)
YDOM = last.(dom)
# mesh = refine(mesh, TriRefinement())

# mesh= refine(mesh,QuadRefinement())
f=Figure()
ax = Axis(f[1,1])

viz!(ax,mesh, showsegments=true)
scatter
scatter!(ax, xm, ym)
# scatter!(ax, XDOM, YDOM)

f
end

waveNumber = 45.0
th = 135
waveVector = waveNumber*SVector(cosd(th), sind(th))

GC.gc()

banded = 2
@time wave_aprx = BoundaryWall.boundaryWallWave(waveVector, xm, ym, XDOM, YDOM, σ, rij, ds, length(ds), diagind(size(rij)...));
@time wave_band = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, banded, Inf);
@time wave_full = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, length(ds), N, Inf);

using CairoMakie


t, m = GeometryUtils.getPlottableMesh(mesh)
let
  GC.gc()
  fig = Figure(figure_padding=0,fonts = (; regular = "Tex Gyre TermesX"), fontsize=24)
  ax = Axis(fig[1,1])
  # viz!(ax, mesh, color=abs2.(wave_band), colormap=:amp)
  # viz!(ax, mesh, color=real(BoundaryWall.incidentWave.(Ref(waveVector),XDOM, YDOM)), colormap=:roma, shading=NoShading)
  Makie.mesh!(ax, t,m, color=abs2.(wave_band), colormap=:batlow)
  # scatter!(ax, xm, ym, color=:white)
  [lines!(ax,GeometryUtils.createEllipse(R,R, T,0.0, c)[1:2]..., color=(:white,0.5)) for c in centers]
  # Colorbar(fig[1,2], labelrotation=0,colorrange=(0,1),colormap=:batlow, ticks = (0:1, ["Min", "Max"]))
  text!(ax, -0.85, 1.25, text=L"|\psi|^2", fontsize=40)
  ax.aspect=DataAspect()
  hidedecorations!(ax)
  hidespines!(ax)
  # save("wave.pdf", fig)
  fig
end

let
  fig = Figure()
  ax = [Axis(fig[1,1], title="Riemann Sum"), 
        Axis(fig[1,2], title="banded (using |i-j|<$banded)"),
        Axis(fig[1,3], title="Fully Integrated")]
  viz!(ax[1], MESH, color=abs2.(wave_aprx), colormap=:amp)
  viz!(ax[2], MESH, color=abs2.(wave_band), colormap=:amp)
  viz!(ax[3], MESH, color=abs2.(wave_full), colormap=:amp)
  [ax.aspect=DataAspect() for ax in ax]
  # [lines!(ax, GeometryUtils.createCircle(R, LinRange(-pi, pi, 50), c)[1:2]..., color=:black) for ax in ax, c in centers]
  # [lines!(ax, GeometryUtils.createEllipse(R,R, T, -i*pi/2, c)[1:2]..., color=:black) for ax in ax, (c,i) in zip(centers, LinRange(0,1, length(centers)))]
  [scatter!(ax, x, y) for ax in ax]
  fig
end

# 3d ima(ge
# surface(xdom, ydom, abs2.(reshape(wave_full, Nx, Ny)), colormap=:amp)

let 
  fig = Figure()
  ax = Axis3(fig[1,1])
  heatmap!(ax, xdom, ydom, abs2.(reshape(wave_full, Nx, Ny)), colormap=:amp, interpolate=true)

  band!(ax,Makie.Point3f.(xdom[1], ydom, 0.0), Makie.Point3f.(xdom[1], ydom, abs2.(reshape(wave_full, Nx, Ny)[1,:])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom[end], ydom, 0.0), Makie.Point3f.(xdom[end], ydom, abs2.(reshape(wave_full, Nx, Ny)[end,:])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom, ydom[end], 0.0), Makie.Point3f.(xdom, ydom[end], abs2.(reshape(wave_full, Nx, Ny)[:,end])), fxaa=true)
  band!(ax,Makie.Point3f.(xdom, ydom[1], 0.0), Makie.Point3f.(xdom, ydom[1], abs2.(reshape(wave_full, Nx, Ny)[:,1])), fxaa=true)
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