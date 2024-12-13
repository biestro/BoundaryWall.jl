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

begin
N = 10
R = 1.0
σ = -0.25im
T = LinRange(-pi, pi, N+1)
centers = GeometryUtils.buildGrid(GeometryUtils.RectangularGrid((0.0, 0.0),5.0, 5.0), 5,5)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)

# rotate centers
R_mat(θ::Float64) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
ϕ = pi/5
RMAT = SMatrix{2,2}(R_mat(ϕ))
centers = [RMAT * SVector(_x, _y) for (_x, _y) in zip(first.(centers), last.(centers))]



INDICES = sort(sample(eachindex(centers), replace=false, length(centers)÷5))
deleteat!(centers, tuple(INDICES...))
centers_3d = map(x -> [x[1], x[2], 0.0], centers)



circ = [GeometryUtils.createCircle(R, T, c) for c in centers]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)


# rotated = [R_mat(ϕ) * [_x; _y] for (_x, _y) in zip(x, y)]
# rotated_mid = [R_mat(ϕ) * [_x; _y] for (_x, _y) in zip(xm, ym)]

# x, y = first.(rotated), last.(rotated)
# xm, ym = first.(rotated_mid), last.(rotated_mid)

f, ax = scatter(xm, ym)
ax.aspect=DataAspect()
f
end




Nx, Ny = 250,250
xdom = LinRange(floor(minimum(x))-1,ceil(maximum(x))+1, Nx)
ydom = LinRange(floor(minimum(y))-1,ceil(maximum(y))+1, Ny)
MESH = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(MESH), MESH.topology)
# COORDS = coords.(centroid.(MESH))
COORDS = coords.(MESH.vertices)
XDOM, YDOM = first.(COORDS), last.(COORDS)


waveNumber = 2.0
th =0
waveVector = waveNumber*SVector(cosd(th), sind(th))


GC.gc()

banded = round(Int, size(rij,1) * 5/100)
@time wave_aprx = BoundaryWall.boundaryWallWave(waveVector, xm, ym, XDOM, YDOM, σ, rij, ds, length(ds), diagind(size(rij)...), 5.0);
@time wave_band = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, rij, length(ds), N, banded, 5.0);
@time wave_full = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, XDOM, YDOM, σ, ds, length(ds), N, 5.0);

let
  fig = Figure(fontsize=24)
  ax = [Axis(fig[1,1], title="discretized"), 
        Axis(fig[1,2], title="banded (using |i-j|<$banded)"),
        Axis(fig[1,3], title="integrated")]
  viz!(ax[1], MESH, color=abs2.(wave_aprx), colormap=:turbo)
  viz!(ax[2], MESH, color=abs2.(wave_band), colormap=:turbo)
  viz!(ax[3], MESH, color=abs2.(wave_full), colormap=:turbo)
  [ax.aspect=DataAspect() for ax in ax]
  # [lines!(ax, GeometryUtils.createCircle(0.5, LinRange(-pi/2, pi/2, 50), c)[1,2]) for (ax, c) in zip(ax, centers)]
  [
    [
      lines!(ax, getindex.(Ref(GeometryUtils.createCircle(R, LinRange(-pi, pi, 150), c)), (1, 2))..., color=(:white,1.0), linewidth=1.5)
      for c in centers
    ] 
    for ax in ax
  ]
  fig
end

_abaqus_mesh = GeometryUtils.circularLatticeMesh(centers_3d, R, x, y, (100,100))

_coords = coords.(_abaqus_mesh.vertices)
_x, _y = first.(_coords), last.(_coords)

@time wave_band_mesh = BoundaryWall.boundaryWallWave(waveVector, x, y, xm, ym, _x, _y, σ, ds, rij, length(ds), N, banded, 5.0);


Ex0(k::SVector{2, Float64}, x::Float64, y::Float64) = BoundaryWall.polarizedField(k, x, y, 0.0)
Ey0(k::SVector{2, Float64}, x::Float64, y::Float64) = BoundaryWall.polarizedField(k, x, y, 0.0)
Ez0(k::SVector{2, Float64}, x::Float64, y::Float64) = BoundaryWall.polarizedField(k, x, y, 0.0)

Ex, Ey, Ez  = BoundaryWall.boundaryWallVec(waveVector, xm, ym, x, y, _x, _y, -0.25im, rij, ds, length(ds), N, SVector(Ex0, Ey0, Ez0) )

let 
  fig = Figure(fontsize=24)
  ax = [Axis(fig[1,1], title="Meshed domain (band integrated)"), 
        Axis(fig[1,2], title="Band integrated"), 
        Axis(fig[1,3], title="Riemann sum")]
  viz!(ax[1], _abaqus_mesh, color=abs2.(wave_band_mesh))
  viz!(ax[2], MESH, color=abs2.(wave_band))
  [lines!(ax[2], getindex.(Ref(GeometryUtils.createCircle(R, LinRange(-pi, pi, 150), c)), (1, 2))..., color=(:white,1.0), linewidth=1.5)
      for c in centers] 
  [ax.aspect=DataAspect() for ax in ax]
  heatmap!(ax[3], xdom, ydom, abs2.(reshape(wave_aprx, Nx, Ny)))
  [lines!(ax[3], getindex.(Ref(GeometryUtils.createCircle(R, LinRange(-pi, pi, 150), c)), (1, 2))..., color=(:white,1.0), linewidth=1.5)
      for c in centers] 
  xlims!(ax[1],minimum(_x), maximum(_x))
  ylims!(ax[1],minimum(_y), maximum(_y))
  xlims!(ax[2],minimum(XDOM), maximum(XDOM))
  ylims!(ax[2],minimum(YDOM), maximum(YDOM))
  fig
end

stokes = BoundaryWall.calcStokes.(Ex, Ey)

GLMakie.activate!()
let 
  fig = Figure(fontsize=24)
  ax = [Axis(fig[1,1], title="Ex"), 
        Axis(fig[1,2], title="Ey"), 
        Axis(fig[1,3], title="Ez")]
  viz!(ax[1], _abaqus_mesh, color=abs2.(Ex), colormap=:turbo, colorrange=(minimum(abs2.(Ex)), maximum(abs2.(Ex))))
  viz!(ax[2], _abaqus_mesh, color=abs2.(Ey), colormap=:turbo)
  viz!(ax[3], _abaqus_mesh, color=abs2.(Ez), colormap=:turbo)
  [ax.aspect=DataAspect() for ax in ax]
  fig
end


_c, _t = GeometryUtils.getPlottableMesh(_abaqus_mesh)

using GLMakie
GLMakie.activate!()
let fig=Figure(fontsize=14,size=(700,600), theme=theme_latexfonts())
  ga = fig[1,1] = GridLayout()
  gb = fig[1,2] = GridLayout()
  gc = fig[2,1] = GridLayout()
  gd = fig[2,2] = GridLayout()
  ax,hm=Makie.mesh(ga[1,1],_c, _t, color=getindex.(stokes,1),shading=NoShading, colormap=:amp,                            ); Colorbar(ga[1,2], hm)
  ax.title=L"S_0"; ax.aspect=DataAspect(); rowsize!(ga, 1, Aspect(1, 1))
  ax,hm=Makie.mesh(gb[1,1],_c, _t, color=getindex.(stokes,2),shading=NoShading, colormap=:balance, colorrange=(-0.04,0.04)); Colorbar(gb[1,2], hm, ticks=[-0.04, 0.0, 0.04])
  ax.title=L"S_1/S_0"; ax.aspect=DataAspect(); rowsize!(gb, 1, Aspect(1, 1))
  ax,hm=Makie.mesh(gc[1,1],_c, _t, color=getindex.(stokes,3),shading=NoShading, colormap=:amp,                            ); Colorbar(gc[1,2], hm)
  ax.title=L"S_2/S_0"; ax.aspect=DataAspect(); rowsize!(gc, 1, Aspect(1, 1))
  ax,hm=Makie.mesh(gd[1,1],_c, _t, color=getindex.(stokes,4),shading=NoShading, colormap=:balance, colorrange=(-0.02,0.02)); Colorbar(gd[1,2], hm, ticks=[-0.02, 0.0, 0.02])
  ax.title=L"S_3/S_0"; ax.aspect=DataAspect(); rowsize!(gd, 1, Aspect(1, 1))
  
  
  save("./images2/stokes_linear.png", fig, px_per_unit = 3)
  fig
end