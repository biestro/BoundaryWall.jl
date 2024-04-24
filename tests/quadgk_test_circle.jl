
# testing quadrature for TMATRIX

using GLMakie, StaticArrays, LinearAlgebra
using Meshes
using StatsBase: sample
using Distributions: Uniform
using CircularArrays
using QuadGK
using ThreadsX
using Meshes
using CircularArrays



include("../src/GeometryUtils.jl")
include("../src/BoundaryWall.jl")

# """
# Modify this for banded Mij
# """
# LinearAlgebra.diagind(a::Int64, b::Int64, c::UnitRange{Int64}) = mapreduce(i->diagind(a,b,i), vcat, c)


# for n in 10:100
begin
N = 150
# parab = GeometryUtils.confocalBilliard(1.0, 2.0/3.0,N)
parab = GeometryUtils.createEllipse(1.0, 0.5, LinRange(-pi,pi, N), 0.0,(0.0, 0.0))

x       = getindex(parab, 1)
y       = getindex(parab, 2)
x,y     = GeometryUtils.divideCurve(x,y, N)
xm, ym  = GeometryUtils.calcMidpoints(x,y)
ds      = GeometryUtils.calcArcLength(x,y)
rij     = GeometryUtils.calcDistances(xm, ym)

rijMid     = rij
arcLengths = ds
dindex  = diagind(size(rijMid)...)

σ          = -0.25im

# positions as SVectors
rMid = SVector.(xm, ym)
rPos = CircularArray(SVector.(x,y)[1:end-1])
# arr  = [vectorPos(rPos[j], rPos[j+1] - rPos[j], t) for j in 1:N, t in 0:0.2:1]
# scatter(arr[:], color=(:blue, 0.2))
lines(x,y)
end

T = -inv([BoundaryWall.calcMmatrix(waveNumber, σ, rMid[i], rPos[j], rPos[j+1]) for i in 1:N, j in 1:N])

heatmap( circshift(abs2.(T),( N÷2, N÷2)), colorscale=log10, colormap=:grayC)

Nx, Ny = 150,150

MESH = GeometryUtils.confocalBilliardMesh(1.0, 2.0/3.0, 20, Nx, Ny)
Z    = coordinates.(MESH.vertices)
XDOM, YDOM = first.(Z), last.(Z)

xdom = LinRange(floor(minimum(x))-1.2,ceil(maximum(x))+1.2, Nx)
ydom = LinRange(floor(minimum(y))-1.2,ceil(maximum(y))+1.2, Ny)
MESH = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(MESH), MESH.topology)
# COORDS = coordinates.(centroid.(MESH))
COORDS = coordinates.(MESH.vertices)

XDOM, YDOM = first.(COORDS), last.(COORDS)

# viz(MESH, showsegments=true)

GC.gc()
waveNumber = 7.0156	
th = 45
waveVector = waveNumber*SVector(cosd(th), sind(th))

wave = BoundaryWall.boundaryWallWave(waveVector, xm, ym, XDOM, YDOM, σ, rij, ds, length(ds), diagind(size(rij)...))

wave2 = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, 10, Inf)
# using BenchmarkTools

# @benchmark BoundaryWall.calcMmatrix(waveNumber, σ, rij, ds, 6, rMid, rPos)
# @benchmark BoundaryWall.calcMmatrix(waveNumber, σ, rij, ds, dindex)
# @benchmark Mij  = [BoundaryWall.calcMmatrix(waveNumber, σ, rMid[i], rPos[j], rPos[j+1]) for i in 1:N, j in 1:N]

let
  T = 1
  fig = Figure(size=(700,600))
  ax = Axis3(fig[1,1], elevation=pi/2 * 0.5, azimuth = -3pi/4 + 0.2, title="ψ₁")
#   hidedecorations!(ax)
#   hidespines!(ax)
#   ax.aspect=DataAspect()
  zlims!(ax, -2,2)
  ax.zlabel="Re ψ(t)"
  ax.zlabelrotation=0
#   for t in 0:0.01:T
  ax.aspect = (1.0, 1.0, 0.5)
  t = 0.0; dt = 0.01
#   Colorbar(fig[1,2], colormap=:balance, colorrange=(-1,1))
  for (i,t) in enumerate(0.0:dt:T-dt)
  empty!(ax)
  
  println("$t")
#   viz!(ax, MESH, color=real(wave2 .* exp(2pi/T*im*t)), colormap=:berlin)
  surface!(ax, xdom, ydom, reshape(real(wave2 .* exp(-2pi*im*t)), Nx, Ny), colormap=:balance, overdraw=false, fxaa=true)
#   text!(ax, -1.0, -1.0, -1.0, text="t=$t")
  lines!(ax, Makie.Point3f.(x, y, 0.01), color=(:black,1.0), linewidth=2, fxaa=true, overdraw=false)
  display(fig)
  sleep(0.001)
#   save("./images/circle/fig"*lpad(i,4,'0')*".png", fig, px_per_unit=2)
#   t += dt
  end
#   fig
end


viz(MESH, color=real(BoundaryWall.incidentWave.(Ref(waveVector), XDOM, YDOM) .* exp(im*1)), colormap=:amp)


let
PSI = abs2.(wave)
PSI2 = abs2.(wave2)
f = Figure()
ax = [Axis(f[1,1],title="Fully integrated"),
      Axis(f[1,2],title="Band integrated, |i-j|<$band")]
# hm=heatmap!(ax,xdom, ydom, reshape(PSI, Nx, Ny),colormap=:dense)
viz!(ax[1], MESH, color=(PSI), colormap=:turbo)
viz!(ax[2], MESH, color=(PSI2), colormap=:turbo)
lines!(ax[1],x,y,color=:black)
lines!(ax[2],x,y,color=:black)
[ax.aspect=DataAspect() for ax in ax]
# cb=Colorbar(f[1,2], hm)
f
end


