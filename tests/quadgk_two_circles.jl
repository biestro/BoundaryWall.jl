

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

N = 10
R = 1.0
σ = -0.25im
T = LinRange(-pi, pi, N+1)
centers = GeometryUtils.buildGrid(GeometryUtils.HexagonalGrid((0.0, 0.0),  3.0), 2,2)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)
INDICES = sort(sample(eachindex(centers), replace=false, length(centers)÷3))
deleteat!(centers, tuple(INDICES...))



circ = [GeometryUtils.createCircle(R, T, c) for c in centers]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)




lines(x,y)

Nx, Ny = 250,250
xdom = LinRange(floor(minimum(x))-1,ceil(maximum(x))+1, Nx)
ydom = LinRange(floor(minimum(y))-1,ceil(maximum(y))+1, Ny)
MESH = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(MESH), MESH.topology)
# COORDS = coords.(centroid.(MESH))
COORDS = coords.(MESH.vertices)
XDOM, YDOM = first.(COORDS), last.(COORDS)


waveNumber = 5.0
th = 0
waveVector = waveNumber*SVector(cosd(th), sind(th))

rMid = SVector.(xm, ym)
rPos = CircularArray(SVector.(x,y)[1:end])

Mij_aprx = BoundaryWall.calcMmatrix(waveNumber, σ, rij, ds, diagind(size(rij)...))
Mij_band = BoundaryWall.calcMmatrix(waveNumber, σ, rij, ds, 1, rMid, rPos, N)
Mij_inte = [BoundaryWall.calcMmatrix(waveNumber, σ, rMid[i], rPos[j], rPos[j+1]) for i in axes(rMid,1), j in BoundaryWall.segmentIterator(rPos, N)]

Mij_inte

T_aprx = -inv(Mij_aprx)
T_band = -inv(Mij_band)
T_inte = -inv(Mij_inte)



let

  fig = Figure()
  ax = Axis(fig[1,1], xgridvisible=false, ygridvisible=false)
  
  sl = SliderGrid(fig[2, 1],
    (label="i", range = eachindex(rMid), startvalue = 1,),
    (label="j", range = eachindex(rMid), startvalue = 1 ),
    (label="i2", range = BoundaryWall.segmentIterator(rPos, N), startvalue = 1,))

  pointj = lift(sl.sliders[1].value) do _j
    Makie.Point2f(rMid[_j])
  end

  pointi = lift(sl.sliders[2].value) do _i
    Makie.Point2f(rMid[_i])
  end

  linei = lift(sl.sliders[3].value) do _i
    [BoundaryWall.r(rPos[_i], rPos[_i+1], t) for t in 0:0.1:1]
  end

  r_i = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    "|r-r'|="*string(round(norm(rMid[_i] - rMid[_j]), digits=4))
  end

  line_ij = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    [rMid[_i], rMid[_j]]
  end
  
  # M_ij = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
  #   # string(round(M(rMid[_i], rPos[_i], rPos[_i+1]), digits=4))
  #   "1st Ord.: "*string(round(G0(waveNumber, rMid[_i], rMid[_j])*norm(rPos[_i] - rPos[_i+1]), digits=4))
  # end

  # M_ij_quad = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
  #   "QuadGK: "*string(round(M(waveNumber, rMid[_i], rPos[_j], rPos[_j+1]), digits=4))
  #   # string(round(G0(waveNumber, rMid[_i], rMid[_j])*ds[_j], digits=4))
  # end

  # lines!(ax, x2, y2, color=(:green, 0.5), linewidth=2)
  lines!(ax, x, y, color=(:black, 0.5), linewidth=2)
  lines!(ax, line_ij, color=:black, linewidth=2)
  lines!(ax, linei, color=:red, linewidth=4)
  scatter!(ax, pointj, color=:blue, markersize=10)
  scatter!(ax, pointi, color=:red, markersize=10)
  # text!(ax, -2.8,-1., text=r_i)
  # text!(ax, -2.8,-1.5, text=M_ij)
  # text!(ax, -2.8,-2.0, text=M_ij_quad)
  ax.aspect=DataAspect()
  fig
end


let
  fig = Figure()
  ax = [Axis(fig[1,1], title="fully-approximated"),
  Axis(fig[1,2], title="band-integrated"),
  Axis(fig[1,3], title="integrated")]
  heatmap!(ax[1], abs2.(T_aprx), colorrange=(0,1), colormap=:amp)
  heatmap!(ax[2], abs2.(T_band), colorrange=(0,1), colormap=:amp)
  heatmap!(ax[3], abs2.(T_inte), colorrange=(0,1), colormap=:amp)
  [ax.aspect=DataAspect() for ax in ax]
  fig
end

GC.gc()

banded = round(Int, size(rij,1) * 5/100)
@time wave_aprx = BoundaryWall.boundaryWallWave(waveVector, xm, ym, XDOM, YDOM, σ, rij, ds, length(ds), diagind(size(rij)...));
@time wave_band = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, banded, Inf);
@time wave_full = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, length(ds), N, Inf);

let
  fig = Figure()
  ax = [Axis(fig[1,1], title="discretized"), 
        Axis(fig[1,2], title="banded (using |i-j|<$banded)"),
        Axis(fig[1,3], title="integrated")]
  viz!(ax[1], MESH, color=abs2.(wave_aprx), colormap=:balance)
  viz!(ax[2], MESH, color=abs2.(wave_band), colormap=:balance)
  viz!(ax[3], MESH, color=abs2.(wave_full), colormap=:balance)
  [ax.aspect=DataAspect() for ax in ax]
  # [lines!(ax, GeometryUtils.createCircle(0.5, LinRange(-pi/2, pi/2, 50), c)[1,2]) for (ax, c) in zip(ax, centers)]
  [
    [
      lines!(ax, getindex.(Ref(GeometryUtils.createCircle(R, LinRange(-pi, pi, 150), c)), (1, 2))..., color=(:black,0.0), linewidth=1.5)
      for c in centers
    ] 
    for ax in ax
  ]
  fig
end