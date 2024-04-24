# TODO: update boundaryWallWave to take r::Vector{SVector{2,Float64}} instead of two vectors x,y

using GLMakie
using StaticArrays: SVector
using LinearAlgebra
using Meshes

include("../src/BWM_geometry.jl")
include("../src/BWM_main.jl")

begin # definitions
HBAR        = 1.0
MASS        = HBAR/2
SIGMA       = (2*MASS/HBAR^2)*(1/4*im)
N           = 100
NX          = 150
NY          = 150
ϕ           = 155
waveVector  = 5*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
end
################################################################################
# # parametric curve 
# g(t::Float64, a::Float64) = (a+sin(4*t)) .* (cos(t), sin(t))
# g(t::Float64, a::Float64) = (a * cos(t), a * sin(t))

# θ = LinRange(-pi, pi, N)

# z = g.(θ, 4.0)
# x, y            = first.(z), last.(z)
# x, y            = GeometryUtils.divideCurve(x, y, N)
# xm, ym          = GeometryUtils.calcMidpoints(x,y)
# arc_lengths     = GeometryUtils.calcArcLength(x,y)
# distance_matrix = GeometryUtils.calcDistances(xm,ym)
################################################################################

# builtin geometry (parabolic billiard)
centers = ((LinRange(0.0, 1.0, 4)).^0.5) * 3.0
barrier = [GeometryUtils.createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), 0.0, SVector(x0, 0.0)) for x0 in centers]
x  = vcat(getindex.(barrier, 1)...)
y  = vcat(getindex.(barrier, 2)...)
xm = vcat(getindex.(barrier, 3)...)
ym = vcat(getindex.(barrier, 4)...)
ds = vcat(getindex.(barrier, 5)...)
rij = GeometryUtils.calcDistances(xm, ym)

lines(x,y)

# domain
x0, xf = (-5, 5)
y0, yf = (-5, 5)
xdom = LinRange(x0, xf, NX)
ydom = LinRange(y0, yf, NY)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coordinates.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)

banded = 2


@time wave = BoundaryWall.boundaryWallWave(waveVector, (k,r)->BoundaryWall.gaussianWave(k,r - SVector(0.0, 0.0), 10.0; abstol=1e-7), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, 10.0);

let 
  fig = Figure(size=(600,600))
  ax = Axis(fig[1,1],
        xtickalign=1,ytickalign=1,
        xticksmirrored=1,yticksmirrored=1,
        xgridvisible=false, ygridvisible=false,
        # xticklabelsvisible=false,yticklabelsvisible=false,
        xminorticksvisible=true,yminorticksvisible=true,
        xminortickalign=1,
        yminortickalign=1,
        xminorticks=IntervalsBetween(2),
        yminorticks=IntervalsBetween(2),
        )
  # heatmap!(ax, xdom[:], ydom[:], real(wave), colormap=:dense)
  # scatter!(ax, XDOM, YDOM, color=real)
  ax.title = "Band integrated |i - j|<$banded"
  viz!(ax, MESH, color=abs2.(wave), shading=NoShading, colormap=:linear_kbgyw_5_98_c62_n256)
  # viz!(ax, MESH, color=abs2.(BoundaryWall.gaussianWave.(Ref(waveVector), COORDS, 10.0; abstol=1e-4)), shading=NoShading, colormap=:linear_kbgyw_5_98_c62_n256)
  # viz!(ax, MESH, color=real(wave), shading=NoShading, colormap=:Spectral)
  
  [lines!(ax, GeometryUtils.createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), 0.0, SVector(x0, 0.0))[1:2]...,color=(:white,0.5)) for x0 in centers]
  ax.aspect=DataAspect() 
  # save("")
  fig
end

# quiver plot
begin

xp, yp, u, v = BoundaryWall.gradient(xdom, ydom, wave)

strength = sqrt.(u.^2 + v.^2)

  k = 5

  XP = xp[1:k:end, 1:k:end]
  YP = yp[1:k:end, 1:k:end]
  U  = u[1:k:end, 1:k:end]
  V  = v[1:k:end, 1:k:end]
  S  = strength[1:k:end, 1:k:end]

  # f,ax=contour(xs,y)
  f = Figure(theme=theme_black())
  ax = Axis(f[1,1])
  
  arrows!(ax, Makie.Point2f.(XP, YP)[:], 
              Makie.Point2f.(U,V)[:], 
              align=:center,
              normalize=false,
              lengthscale=0.01,
              color=S[:],
              colormap=:cividis)
              
  # scatter!(ax, XP[:], YP[:])
  ax.aspect=DataAspect()
  f
end

# animation
let 
  fig = Figure(size=(600,600))
  ax = Axis(fig[1,1],
        xtickalign=1,ytickalign=1,
        xticksmirrored=1,yticksmirrored=1,
        xgridvisible=false, ygridvisible=false,
        # xticklabelsvisible=false,yticklabelsvisible=false,
        xminorticksvisible=true,yminorticksvisible=true,
        xminortickalign=1,
        yminortickalign=1,
        xminorticks=IntervalsBetween(2),
        yminorticks=IntervalsBetween(2),
        )
  # heatmap!(ax, xdom[:], ydom[:], real(wave), colormap=:dense)
  # scatter!(ax, XDOM, YDOM, color=real)
  ax.title = "Band integrated |i - j|<$banded"
  # viz!(ax, MESH, color=abs2.(wave), shading=NoShading, colormap=:linear_kbgyw_5_98_c62_n256)
  viz!(ax, MESH, color=real(wave .* exp(-im)), shading=NoShading, colormap=:Spectral)
  
  [lines!(ax, GeometryUtils.createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), 0.0, SVector(x0, 0.0))[1:2]...,color=(:black,0.5)) for x0 in centers]
  ax.aspect=DataAspect() 
  # save("")
  fig
end



