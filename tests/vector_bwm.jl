using GLMakie, StaticArrays, LinearAlgebra
using CircularArrays
import BoundaryWall as BWM

N = 100
R = 1.0
σ = -0.25im
T = LinRange(pi/2, 3pi/2, N)
# circ = [BWM.createCircle(R, T, SVector(0.0, i)) for i in [0.0, 3.0, 6.0]]
circ = [BWM.createCircle(R, T, SVector(0.0, 0.0))]
circ = [BWM.createEllipse(R,R*0.3, T,0.0, SVector(0.0, 0.0))]
x  = vcat(getindex.(circ, 1)...)
y  = vcat(getindex.(circ, 2)...)
xm = vcat(getindex.(circ, 3)...)
ym = vcat(getindex.(circ, 4)...)
ds = vcat(getindex.(circ, 5)...)
rij = BWM.calcDistances(xm, ym)

y,x,ym,xm,rij,ds=createConfocalBilliard(3.0, 2.0, N)


rMid = SVector.(xm, ym)
rPos = CircularArray(SVector.(x,y)[1:end-1])

lines(x,y)

waveNumber = sqrt(1.802926)
th = 0.0
waveVector = SVector(cosd(th), sind(th)) * waveNumber

_g = BoundaryWall.calcGreenTensor(waveVector,-0.25im,rij,ds,4,rMid,rPos,length(ds)) # greens tensor

γ = 10000.0
greenMatrix = [I(N) - γ*_g[1] -γ*_g[2] -γ*_g[3]; 
               -γ*_g[2] I(N) - γ*_g[4] -γ*_g[5]; 
               -γ*_g[3] -γ*_g[5] I(N) - γ*_g[6]]
heatmap(abs2.(-inv(greenMatrix)), colorscale=log10)

Ex0(k::SVector{2, Float64}, x::Float64, y::Float64) = BWM.polarizedField(k, x, y, 0.0)
Ey0(k::SVector{2, Float64}, x::Float64, y::Float64) = BWM.polarizedField(k, x, y, 0.0)
Ez0(k::SVector{2, Float64}, x::Float64, y::Float64) = BWM.polarizedField(k, x, y, 0.0)

Nx = Ny = 100

using Meshes
xdom = collect(LinRange(-6,6,Nx))
ydom = collect(LinRange(-3,5.5,Ny))

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coordinates.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)

banded = 2

@time Ex, Ey, Ez = BWM.boundaryWallVec(waveVector, xm, ym, x, y, XDOM, YDOM, σ, rij, ds, length(ds), N, SVector(Ex0, Ey0, Ez0));

Etot = mapreduce(x->x.*conj(x), +, [Ex, Ey, Ez])
TE = mapreduce(x->abs2.(x), +, [Ex, Ey])
# TM = mapreduce(x->x.*conj(x), +, [Ex, Ey])
TM = abs2.(Ez)

let 
  fig = Figure()
  ax = [Axis(fig[1,1]), Axis(fig[1,2])]
  heatmap!(ax[1],xdom, ydom, reshape(TE, Nx, Ny), colorrange=(1.991, 2.02),colormap=Reverse(:grayC), highclip = :red, lowclip = :blue)
  heatmap!(ax[2],xdom, ydom, reshape(TM, Nx, Ny), colormap=Reverse(:grayC))
  [ax.aspect=DataAspect() for ax in ax]
  fig
end
S = BWM.calcStokes.(Ex, Ey)

SIGMA = -0.25im

S1 = getindex.(S,2)

heatmap(xdom, ydom, reshape(S1, Nx, Ny))