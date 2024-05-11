# TODO: update boundaryWallWave to take r::Vector{SVector{2,Float64}} instead of two vectors x,y

using GLMakie
using StaticArrays: SVector
using LinearAlgebra
using Meshes

using BoundaryWall

begin # definitions
HBAR        = 1.0
MASS        = HBAR/2
SIGMA       = (2*MASS/HBAR^2)*(1/4*im)
N           = 100
NX          = 235
NY          = 235
ϕ           = 180
R           = 4.0
waveVector  = 30*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
end

centers = [SVector(0.0, 0.0), SVector(4R, 0.0), SVector(0.0,-4R)]
# barrier = [createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), th, SVector(x0, 0.0)) for (x0, th) in zip([0.0, 2.0], [pi/6, pi/4])]
barrier = [createEllipse(0.0, R, LinRange(-pi/2,pi/2,N), th, x0) for (x0,th) in zip(centers, [pi/4, 1pi/8, 3pi/8])]
x  = vcat(getindex.(barrier, 1)...)
y  = vcat(getindex.(barrier, 2)...)
xm = vcat(getindex.(barrier, 3)...)
ym = vcat(getindex.(barrier, 4)...)
ds = vcat(getindex.(barrier, 5)...)
rij = calcDistances(xm, ym)

lines(x,y)

# domain
x0, xf = (-5, 20)
y0, yf = (-20, 5)
xdom = LinRange(x0, xf, NX)
ydom = LinRange(y0, yf, NY)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coordinates.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)

banded = 5


@time wave = boundaryWallWave(waveVector, (k,r)->gaussianWave(k,r - SVector(0.0, 0.0), 45.0; abstol=1e-8), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, 
             Inf);

# wave path
xp, yp, u, v = BoundaryWall.gradient(xdom, ydom, wave)


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
  ax.title = "Band integrated |i - j|<$banded"
  # viz!(ax, MESH, color=abs2.(wave), shading=NoShading, colormap=:turbo)
  heatmap!(ax,xdom, ydom, abs2.(reshape(wave, NX, NY)), interpolate=false, colormap=:balance)

  k = 8
  

  XP = xp[1:k:end, 1:k:end]
  YP = yp[1:k:end, 1:k:end]
  U  = u[1:k:end, 1:k:end]
  V  = v[1:k:end, 1:k:end]
  # arrows!(ax, Makie.Point2f.(XP, YP)[:],Makie.Point2f.(U,V)[:], normalize=true,lengthscale=1,color=(:white, 0.5))
  # [lines!(ax, createEllipse(0.0, 7.0, LinRange(-pi/2,pi/2,N), pi/4, SVector(x0, 0.0))[1:2]...,color=(:white,1.0)) for x0 in centers[1]]
  # [lines!(ax, createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), th, SVector(x0, 0.0))[1:2]..., color=:white) for (x0, th) in zip([0.0, 2.0], [pi/6, pi/4])]
  # scatter!(ax, x, y, color=:white)
  # xlims!(ax, xdom[1], xdom[end])
  # ylims!(ax, ydom[1], ydom[end])
  ax.aspect=DataAspect() 
  # save("")
  fig
end

# quiver plot
θ = LinRange(-pi, pi/2, 300)
x,y = cos.(θ),sin.(θ)

x,y = divideResonatorCurve(x,y, 300)
xm,ym = calcMidpoints(x, y)
ds = calcArcLength(x,y)
rij = calcDistances(xm, ym)
lines(ds)
