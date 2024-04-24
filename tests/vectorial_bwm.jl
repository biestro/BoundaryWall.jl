# TODO: broadcast the calculation of greens tensor (i.e. remove the if else)

using LinearAlgebra
using StaticArrays: SVector
using SpecialFunctions
using GLMakie

include("./src/GeometryUtils.jl")
include("./src/BoundaryWall.jl")




using Meshes

Nx = 100
Ny = 50
xdom = collect(LinRange(-2,2,Nx))
ydom = collect(LinRange(-2,2,Ny))

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coordinates.(vertices(GRID))

XDOM, YDOM = first.(COORDS), last.(COORDS)

xDomain = XDOM;
yDomain = YDOM;


N = 500
R  = 1.0
TH = LinRange(-pi, pi, N+1)
circ = GeometryUtils.createEllipse(R, R, TH, 1.0pi, (0.0, 0.0))
# circ = GeometryUtils.buildConfocalBilliard(3.0, 2.0, N)
x  = getindex(circ,1)
y  = getindex(circ,2)
xm = getindex(circ,3)
ym = getindex(circ,4)
ds = getindex(circ,5)
# ds = getindex(circ,6)
rij = GeometryUtils.calcDistances(xm, ym)

# waveNumber = 10.1735# 3.8317
waveNumber = 3.8317
waveAngle  = pi/2
waveVector = waveNumber * SVector(cos(waveAngle), sin(waveAngle))

include("./src/BoundaryWall.jl");

function polarizedField(k::SVector{2, Float64}, x::Float64, y::Float64, ϕ::Float64)
  return exp(im*ϕ)*exp(-im * (k[1] * x + k[2] * y))
end


# circular polarization
Ex02(k::SVector{2, Float64}, x::Float64, y::Float64) = polarizedField(k, x, y, 0.0)
Ey02(k::SVector{2, Float64}, x::Float64, y::Float64) = polarizedField(k, x, y, 0.0)
Ez02(k::SVector{2, Float64}, x::Float64, y::Float64) = polarizedField(k, x, y, 0.0)


Ex, Ey, Ez = BoundaryWall.boundaryWallVec(waveVector, 
                                          xm,
                                          ym,
                                          XDOM, 
                                          YDOM, 
                                          -0.25im, 
                                          rij,ds, 
                                          length(ds), 
                                          diagind(size(rij)...),
                                          SVector(Ex02, Ey02, Ez02))


# heatmap(xdom, ydom, abs2.(reshape(Ey,Nx,Ny)), colormap=:turbo, colorrange=(0,2))

stokes = BoundaryWall.calcStokes.(Ex, Ey)

s0 = getindex.(stokes,1)
s1 = getindex.(stokes,2)
s2 = getindex.(stokes,3)
s3 = getindex.(stokes,4)
heatmap(xdom, ydom, reshape(abs2.(Ey), Nx, Ny), colorrange=(0,2))
let 
  # splot[abs.(splot) .> 2.0] .= NaN
  fig = Figure()
  ga = fig[1,1] = GridLayout()
  gb = fig[1,2] = GridLayout()
  gc = fig[2,1] = GridLayout()
  gd = fig[2,2] = GridLayout()
  aa, ha=heatmap(ga[1,1], xdom, ydom, reshape(s0, Nx, Ny), interpolate=true, colormap=:amp, colorrange=(0,3))
  ab, hb=heatmap(gb[1,1], xdom, ydom, reshape(s1, Nx, Ny), interpolate=true, colormap=:balance, colorrange=(-1,1) .* (maximum(s1) - minimum(s1)))
  ac, hc=heatmap(gc[1,1], xdom, ydom, reshape(s2, Nx, Ny), interpolate=true, colormap=:balance, colorrange=(-1,1) .* (maximum(s2) - minimum(s2)))
  ad, hd=heatmap(gd[1,1], xdom, ydom, reshape(s3, Nx, Ny), interpolate=true, colormap=:balance, colorrange=(-1,1) .* (maximum(s3) - minimum(s3)))
  Colorbar(ga[1,0], ha, flipaxis=false)
  Colorbar(gb[1,2], hb, flipaxis=true)
  Colorbar(gc[1,0], hc, flipaxis=false)
  Colorbar(gd[1,2], hd, flipaxis=true)
  
  
  [lines!(ax, xm, ym, color=:black) for ax in [aa, ab, ac, ad]]
  [ax.aspect=DataAspect() for ax in [aa, ab, ac, ad]]
  ab.xtrimspine=true
  # Colorbar(fig[1,2], hm)
  fig
end

using CSV, DataFrames

df_stokes = DataFrame([s0 s1 s2 s3], [:s0, :s1, :s2,:s3])

CSV.write("stokes.csv", df_stokes)