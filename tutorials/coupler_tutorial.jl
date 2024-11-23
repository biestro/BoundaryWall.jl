
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
NX          = 150
NY          = 150
ϕ           = 0.0
waveVector  = -8*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
TH          = LinRange(0,pi/2,N)
end

circle = [createEllipse(4.5, 5.0,TH,0.0, SVector(0.0, 0.0)); createEllipse(4.0, 4.0,reverse(TH),0.0, SVector(0.0, 0.0))]


x  = vcat( getindex.(circle, 1)...)
y  = vcat( getindex.(circle, 2)...)
xm = vcat( getindex.(circle, 3)...)
ym = vcat( getindex.(circle, 4)...)
ds = vcat( getindex.(circle, 5)...)
rij = calcDistances(xm, ym)

f,ax= scatter(x,y, color=1:N,colormap=:thermal); ax.aspect=DataAspect(); f

# domain
x0, xf = (-8, 8)
y0, yf = (-8, 8)
xdom = LinRange(x0, xf, NX)
ydom = LinRange(y0, yf, NY)
COORDS = [(x,y) for x  in xdom, y in ydom]
XDOM, YDOM = first.(COORDS)[:], last.(COORDS)[:]

banded = 5


@time wave = boundaryWallWave(SVector(-10.0, 0.0), (k,r)->gaussianWave(k, r - SVector(-0.0, 4.5), 5.0; abstol=1e-8), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, -15.0);

plotFieldDensity(abs2.(wave), xdom, ydom, :turbo)


# wave path
xp, yp, u, v = BoundaryWall.gradient(xdom, ydom, wave)


function plotFieldDensity(f::Vector{Float64}, xdom, ydom, cmap)
  fig = Figure(size=(600,600))
  ax = Axis(fig[1,1],
        xtickalign=1,ytickalign=1,
        xticksmirrored=1,yticksmirrored=1,
        xgridvisible=false, ygridvisible=false,
        # xticklabelsvisible=false,yticklabelsvisible=false,
        xminorticksvisible=true,yminorticksvisible=true,
        xminortickalign=0,
        yminortickalign=0,
        xminorticks=IntervalsBetween(2),
        yminorticks=IntervalsBetween(2),
        )
  ax.title = "Band integrated |i - j|<$banded"
  # viz!(ax, MESH, color=abs2.(wave), shading=NoShading, colormap=:turbo)
  # heatmap!(ax,xdom, ydom, real.(reshape(wave, NX, NY)), interpolate=true, colormap=:balance)
  heatmap!(ax,xdom, ydom, reshape(f, length(xdom), length(ydom)), interpolate=false, colormap=cmap)
  # heatmap!(ax,xdom, ydom, reshape(f, length(xdom), length(ydom)), interpolate=true, colormap=:linear_kbgyw_5_98_c62_n256)

  scatter!(ax, xm, ym, color=:white)

  xlims!(ax, xdom[1], xdom[end])
  ylims!(ax, ydom[1], ydom[end])
  ax.aspect=DataAspect() 
  # save("")
  fig
end


