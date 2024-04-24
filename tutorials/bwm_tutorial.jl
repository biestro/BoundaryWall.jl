using CairoMakie
using StaticArrays: SVector
using LinearAlgebra
using Meshes

include("../src/BWM_geometry.jl")
include("../src/BoundaryWall.jl")

begin # definitions
HBAR        = 1.0
MASS        = HBAR/2
SIGMA       = (2*MASS/HBAR^2)*(1/4*im)
N           = 100
NDOM        = 150
ϕ           = 135
waveVector  = sqrt(1.365)*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
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
y, x,ym, xm, distance_matrix, arc_lengths = BWM_geometry.createConfocalBilliard(2.0, 3.0, N)

# domain
x0, xf = (-8.5, 8.5)
y0, yf = (-7.5, 5)
xdom = LinRange(x0, xf, NDOM)
ydom = LinRange(y0, yf, NDOM)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coordinates.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)

banded = 2


@time wave = BoundaryWall.boundaryWallWave(waveVector, (k,r)->BoundaryWall.planeWave(k,r), x, y, xm, ym, XDOM, YDOM, SIGMA, arc_lengths, distance_matrix, length(arc_lengths), N, banded, Inf);

# xp, yp, u, v = BoundaryWall.gradient(xdom, ydom, wave)


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
  contourf!(ax, xdom, ydom, (abs2.(reshape(wave, NDOM, NDOM))).^0.4, levels=4, colormap=:julia_colorscheme)
  # viz!(ax[3], MESH, color=abs2.(wave), shading=false, colorscheme=:dense)
  
  
  lines!(ax, x, y,color=:white)
  ax.aspect=DataAspect() 
  save("./docs/src/assets/logo.pdf", fig)
  fig
end
