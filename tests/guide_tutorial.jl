"""
Tutorial for a port-like system
"""

using GLMakie, StaticArrays, LinearAlgebra
using Meshes
using StatsBase: sample
using Distributions: Uniform

include("./src/GeometryUtils.jl")
begin
N = 201
# verts = [(4.0, -4.0), (-4.0, -4.0), (-8.0, 1.0)]
verts = [-3.0, -5.0] # bottom corner of square box

opening = [-1.0, 1.0]
c_length = 4.0

theta = LinRange(-pi,pi , N)
x,y   = GeometryUtils.createCircle(0.5, theta, (-3.0, -3.0))
xm, ym = GeometryUtils.calcMidpoints(x,y)
ds = GeometryUtils.calcArcLength(x,y)
# x1,y1   = GeometryUtils.createCircle(0.5, theta, (1.0, 2.0))
# xm1, ym1 = GeometryUtils.calcMidpoints(x1,y1)
# ds1 = GeometryUtils.calcArcLength(x1,y1)
# x = [x; x1]; y = [y; y1]
# xm= [xm; xm1]; ym = [ym; ym1]
# ds = [ds; ds1]

y_barrier = collect(opening[1]:-ds[1]:verts[2])
# x_barrier = ones(length(y_barrier)) * verts[1]
x_barrier = collect(LinRange(verts[1]-c_length,verts[1], length(y_barrier)))

# lower wall
x_barrier_low = collect(x_barrier[end]:ds[1]:-verts[1])
y_barrier_low = ones(length(x_barrier_low)) * y_barrier[end]

# right wall
# y_barrier_right = collect(y_barrier_low[1])
# x_barrier_right = 

x_barrier = [x_barrier;x_barrier_low]
y_barrier = [y_barrier;y_barrier_low]

x_barrier, y_barrier = GeometryUtils.divideResonatorCurve(x_barrier,y_barrier, ds[1])

ds_barrier = GeometryUtils.calcArcLength(x_barrier, y_barrier)

xm_barrier, ym_barrier = GeometryUtils.calcMidpoints(x_barrier, y_barrier)

DS = [ds; ds_barrier; ds_barrier]
XM, YM = [xm; reverse(xm_barrier); xm_barrier], [ym; reverse(ym_barrier); -ym_barrier]

# find distances
# scatter(XM, YM, color=eachindex(ds))
RIJ = GeometryUtils.calcDistances(XM, YM)


fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
scatter!(ax1,XM,YM, color=eachindex(DS))
lines!(ax2, diff(DS))
ax1.aspect=DataAspect()
# ax2.aspect=DataAspect()
fig
end

# lines(diff(DS))

include("./src/BoundaryWall.jl")

begin
  x = LinRange(-2,2,100)
  y = LinRange(-2,2,100)
  z = BoundaryWall.incidentWave.(Ref(SVector(10.0, 0.0)), x, y')
  heatmap(x,y,real(z), colormap=:balance)
end
using Meshes
begin
Nx, Ny = 70,70
# xdom = LinRange(floor(minimum(XM))-1,ceil(maximum(XM))+1, Nx)
# ydom = LinRange(floor(minimum(YM))-1,ceil(maximum(YM))+1, Ny)
xdom = LinRange(-2,2,Nx)
ydom = LinRange(-2,2,Ny)
# xdom = [centers[end][2]+2R, centers[end÷2][2], centers[1][2]-2R]
# ydom = [-1.0, 21.0]

X,Y = GeometryUtils.createCircle(1.0, LinRange(-pi,pi, 500), (0.0,0.0))
XM,YM = GeometryUtils.calcMidpoints(X,Y)
DS = GeometryUtils.calcArcLength(X,Y)
RIJ = GeometryUtils.calcDistances(XM,YM)

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coords.(vertices(GRID))

XDOM, YDOM = first.(COORDS), last.(COORDS)
end
psi  = BoundaryWall.boundaryWallWave(-3.8317*SVector(cosd(0), sind(0)), XM, YM, XDOM, YDOM, -0.25im, RIJ, DS, length(DS), diagind(size(RIJ)...))
psi = reshape(psi, Nx, Ny)

heatmap(abs2.(psi))

using Interpolations

xs = xdom
ys = ydom
zs  = psi
ms = [Makie.Point2f(x,y) for x in xs, y in ys]
# heatmap(v)
itp = linear_interpolation((xs,ys),zs)
grad_psi = [Interpolations.gradient(itp, idx...) for idx in knots(itp)]
dv       = imag( conj(psi) .* grad_psi )
strength = vec(@. sqrt(first(dv) ^ 2 + last(dv)^ 2))



begin
  k = 1
xp = first.(ms[:])[1:k:end]
yp =  last.(ms[:])[1:k:end]
zp =  abs2.(zs[:])[1:k:end]
us = first.(dv[:])[1:k:end]
vs =  last.(dv[:])[1:k:end]
cs =      strength[1:k:end]

# f,ax=contour(xs,y)
f = Figure()
ax = Axis(f[1,1])
arrows!(ax, Makie.Point2f.(xp, yp), Makie.Point2f.(us,vs), normalize=true,lengthscale=0.1,color=zp)
f
end