# TODO: update BWMWave to take r::Vector{SVector{2,Float64}} instead of two vectors x,y

using CairoMakie
CairoMakie.activate!()
using StaticArrays: SVector
using LinearAlgebra


import BoundaryWall as BWM

begin # definitions
HBAR        = 1.0
MASS        = HBAR/2
SIGMA       = (2*MASS/HBAR^2)*(1/4*im)
N           = 100
NX          = 150
NY          = 150
ϕ           = 180.0
waveVector  = 5*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
end

centers = ((LinRange(0.0, 1.0, 4)).^0.5) * 3.0 # barrier = [createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), th, SVector(x0, 0.0)) for (x0, th) in zip([0.0, 2.0], [pi/6, pi/4])]
barrier = [createEllipse(0.0, 8.0, LinRange(-pi/2,pi/2,N), pi/4, SVector(x0, 0.0)) for x0 in centers[1]]
x  = vcat(getindex.(barrier, 1)...)
y  = vcat(getindex.(barrier, 2)...)
xm = vcat(getindex.(barrier, 3)...)
ym = vcat(getindex.(barrier, 4)...)
ds = vcat(getindex.(barrier, 5)...)
rij = calcDistances(xm, ym)

lines(x,y)

# domain
x0, xf = (-8, 8)
y0, yf = (-8, 8)
xdom = LinRange(x0, xf, NX)
ydom = LinRange(y0, yf, NY)
COORDS = [(x,y) for x  in xdom, y in ydom]
XDOM, YDOM = first.(COORDS)[:], last.(COORDS)[:]

banded = 5


@time wave = BWM.boundaryWallWave(waveVector, (k,r)->BWM.gaussianWave(k,r - SVector(0.0, 0.0), 10.0; abstol=1e-8), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, -7.0);

# wave path
xp, yp, u, v = BWM.gradient(xdom, ydom, wave)


begin
  fig = Figure(theme=BWM.theme, size=(600,300))
  gl  = fig[1,1] = GridLayout()
  # ax = [Axis(fig[1,1]), Axis(fig[1,2])]
  wave_density = abs2.(reshape(wave, NX, NY))
  # ax.title = "Beam splitter simulation using negative potential"
  # viz!(ax, MESH, color=abs2.(wave), shading=NoShading, colormap=:turbo)
  # heatmap!(ax,xdom, ydom, real.(reshape(wave, NX, NY)), interpolate=true, colormap=:balance)
  ax1,_=heatmap(gl[1,1],xdom, ydom, wave_density, interpolate=true, colormap=BWM.cmap)
  # heatmap!(ax[2],xdom, ydom, real(reshape(wave, NX, NY)), interpolate=true, colormap=BWM.cmap)
  # arrows!(ax[2], xp, yp, u, v)
  ax1.xautolimitmargin = (0.05f0, 0.05f0)
  ax1.yautolimitmargin = (0.05f0, 0.05f0)
  
  
  k = 6

  XP = xp[1:k:end, 1:k:end]
  YP = yp[1:k:end, 1:k:end]
  U  = u[1:k:end, 1:k:end]
  V  = v[1:k:end, 1:k:end]
  C  = wave_density[1:k:end, 1:k:end]

  ax2,_=arrows(gl[1,2],Makie.Point2f.(XP, YP)[:],Makie.Point2f.(U,V)[:], normalize=true, align=:center,lengthscale=1,color=C[:], colormap=BWM.cmap)
  ax2.xlabel="x"
  ax1.ylabelrotation=0
  ax1.xlabel="x"
  ax1.ylabel="y"
  [lines!(ax1, createEllipse(0.0, 7.0, LinRange(-pi/2,pi/2,N), pi/4, SVector(x0, 0.0))[1:2]...,color=(:black,1.0)) for x0 in centers[1]]
  [lines!(ax2, createEllipse(0.0, 7.0, LinRange(-pi/2,pi/2,N), pi/4, SVector(x0, 0.0))[1:2]...,color=(:black,1.0)) for x0 in centers[1]]
  # [lines!(ax, createEllipse(0.0, 5.0, LinRange(-pi/2,pi/2,N), th, SVector(x0, 0.0))[1:2]..., color=:white) for (x0, th) in zip([0.0, 2.0], [pi/6, pi/4])]
  # xlims!(ax[1], xdom[1], xdom[end])
  # ylims!(ax[1], ydom[1], ydom[end])
  # [ax.aspect=DataAspect()  for ax in [ax]
  # colsize!(gl,1, Aspect(1, 5))
  resize_to_layout!(fig)
  [ax.aspect=DataAspect() for ax in [ax1, ax2]]
  hidedecorations!(ax1, ticks=false,minorticks=false, label=false)
  hidedecorations!(ax2, ticks=false,minorticks=false, label=false)
  save("beam_splitter.svg",fig)
  fig
end