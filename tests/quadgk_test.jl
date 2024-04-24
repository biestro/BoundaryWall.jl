
# testing quadrature for TMATRIX

using GLMakie, StaticArrays, LinearAlgebra
using Meshes
using StatsBase: sample
using Distributions: Uniform
using CircularArrays
using QuadGK


include("../src/GeometryUtils.jl")
include("../src/BoundaryWall.jl")

function r(rMinus::SVector{2, Float64}, rMidpoint::SVector{2, Float64}, rPlus::SVector{2, Float64}, t::Float64)
  if t < 0.5
  return rMinus + t * (rMidpoint - rMinus) 
  end
  return rMidpoint + t * (rPlus - rMidpoint) 
end

r(r0::SVector{2, Float64},r1::SVector{2, Float64}, t::Float64) = r0 + t * (r1 - r0) 
G0(k::Float64, ri::SVector{2, Float64}, rj::SVector{2,Float64}) = BoundaryWall.calcGreenFun(k, norm(ri-rj))
M(k::Float64, ri::SVector{2,Float64}, r0::SVector{2,Float64}, r1::SVector{2, Float64}) = first(quadgk(t -> G0(k, ri, r(r0, r1, t)), 0.0, 0.5, 1.0)) * norm(r1 - r0)  # scaling
# approx_M(ri::Int, j::Int) = G0(waveNumber, ri, rMid[j]) * norm(rPos[j+1] - rPos[j])


# for n in 10:100
begin
N = 15
N2 = 200
ellipse = GeometryUtils.createEllipse(1.0, 1.0, LinRange(-pi, pi, N), pi/4, (0.0, 0.0))

x2, y2,_,_,_ = GeometryUtils.createEllipse(1.0, 2.0, LinRange(-pi, pi, N2), pi/4, (0.0, 0.0))
x       = getindex(ellipse, 1)
y       = getindex(ellipse, 2)
xm      = getindex(ellipse, 3)
ym      = getindex(ellipse, 4)
ds      = getindex(ellipse, 5)
rij     = GeometryUtils.calcDistances(xm, ym)

rijMid     = rij
arcLengths = ds
dindex  = diagind(size(rijMid)...)
waveVector = SVector(2.0, 0.0)
waveNumber = norm(waveVector)
σ          = -0.25im

# positions as SVectors
rMid = SVector.(xm, ym)
rPos = CircularArray(SVector.(x,y)[1:end-1])
# arr  = [vectorPos(rPos[j], rPos[j+1] - rPos[j], t) for j in 1:N, t in 0:0.2:1]
# scatter(arr[:], color=(:blue, 0.2))

end

let


  fig = Figure()
  ax = Axis(fig[1,1], xgridvisible=false, ygridvisible=false)
  
  sl = SliderGrid(fig[2, 1],
    (label="i",  range = 1:N, startvalue = 1,),
    (label="j", range = 1:N, startvalue = 1))

  pointj = lift(sl.sliders[1].value) do _j
    Makie.Point2f(rMid[_j])
  end

  pointi = lift(sl.sliders[2].value) do _i
    Makie.Point2f(rMid[_i])
  end

  linei = lift(sl.sliders[2].value) do _i
    [r(rPos[_i], rPos[_i+1], t) for t in 0:0.1:1]
  end

  r_i = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    "|r-r'|="*string(round(norm(rMid[_i] - rMid[_j]), digits=4))
  end

  line_ij = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    [rMid[_i], rMid[_j]]
  end
  
  M_ij = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    # string(round(M(rMid[_i], rPos[_i], rPos[_i+1]), digits=4))
    "1st Ord.: "*string(round(G0(waveNumber, rMid[_i], rMid[_j])*norm(rPos[_i] - rPos[_i+1]), digits=4))
  end

  M_ij_quad = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    "QuadGK: "*string(round(M(waveNumber, rMid[_i], rPos[_j], rPos[_j+1]), digits=4))
    # string(round(G0(waveNumber, rMid[_i], rMid[_j])*ds[_j], digits=4))
  end

  # lines!(ax, x2, y2, color=(:green, 0.5), linewidth=2)
  lines!(ax, x, y, color=(:black, 0.5), linewidth=2)
  lines!(ax, line_ij, color=:black, linewidth=2)
  lines!(ax, linei, color=:red, linewidth=4)
  scatter!(ax, pointj, color=:blue, markersize=10)
  scatter!(ax, pointi, color=:red, markersize=10)
  text!(ax, -2.8,-1., text=r_i)
  text!(ax, -2.8,-1.5, text=M_ij)
  text!(ax, -2.8,-2.0, text=M_ij_quad)
  ax.aspect=DataAspect()
  fig
end

# for n in 10:100
begin
  N = 205
  parabol = GeometryUtils.confocalBilliard(3.0, 2.0, N)
  
  x       = getindex(parabol, 1)
  y       = getindex(parabol, 2)
  xm, ym  = GeometryUtils.calcMidpoints(x, y)
  ds      = GeometryUtils.calcArcLength(x, y)
  rij     = GeometryUtils.calcDistances(xm, ym)
  
  dindex  = diagind(size(rij)...)
  waveVector = SVector(5.0, 0.0)
  waveNumber = norm(waveVector)
  σ          = -0.25im
  
  # positions as SVectors
  rMid = SVector.(xm, ym)
  rPos = CircularArray(SVector.(x,y)[1:end-1])
  # arr  = [vectorPos(rPos[j], rPos[j+1] - rPos[j], t) for j in 1:N, t in 0:0.2:1]
  # scatter(arr[:], color=(:blue, 0.2))
  
end


# calculate M matrix
# Mij          = @. σ * BoundaryWall.calcGreenFun(waveNumber, rijMid) * arcLengths # multiplies cols by ds[j]
# Mij[dindex] .= @. σ * BoundaryWall.calcGreenFun(waveNumber, arcLengths/2) * arcLengths

# T = -inv(Mij) # calc t matrix
let


  fig = Figure()
  ax = Axis(fig[1,1], xgridvisible=false, ygridvisible=false)
  
  sl = SliderGrid(fig[2, 1],
    (label="i",  range = 1:N, startvalue = 1,),
    (label="j", range = 1:N, startvalue = 1))

  pointj = lift(sl.sliders[1].value) do _j
    Makie.Point2f(rMid[_j])
  end

  pointi = lift(sl.sliders[2].value) do _i
    Makie.Point2f(rMid[_i])
  end

  linei = lift(sl.sliders[2].value) do _i
    [r(rPos[_i], rPos[_i+1], t) for t in 0:0.1:1]
  end

  r_i = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    "|r-r'|="*string(round(norm(rMid[_i] - rMid[_j]), digits=4))
  end

  line_ij = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    [rMid[_i], rMid[_j]]
  end
  
  M_ij = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    # string(round(M(rMid[_i], rPos[_i], rPos[_i+1]), digits=4))
    "1st Ord.: "*string(round(G0(waveNumber, rMid[_i], rMid[_j])*norm(rPos[_i] - rPos[_i+1]), digits=4))
  end

  M_ij_quad = lift(sl.sliders[2].value, sl.sliders[1].value) do _i, _j
    "QuadGK: "*string(round(M(waveNumber,rMid[_i], rPos[_j], rPos[_j+1]), digits=4))
    # string(round(G0(waveNumber, rMid[_i], rMid[_j])*ds[_j], digits=4))
  end

  # lines!(ax, x2, y2, color=(:green, 0.5), linewidth=2)
  lines!(ax, x, y, color=(:black, 0.5), linewidth=2)
  lines!(ax, line_ij, color=:black, linewidth=2)
  lines!(ax, linei, color=:red, linewidth=4)
  scatter!(ax, pointj, color=:blue, markersize=10)
  scatter!(ax, pointi, color=:red, markersize=10)
  text!(ax, -5.8,-4., text=r_i)
  text!(ax, -5.8,-4.5, text=M_ij)
  text!(ax, -5.8,-5.0, text=M_ij_quad)
  ax.aspect=DataAspect()
  fig
end


begin 
fig = Figure()
waveVector = SVector(0.0, sqrt(2.070366))
waveNumber = norm(waveVector)
ax = [Axis(fig[1,1], yreversed=true, title="Fully integrated"),
      Axis(fig[1,2], yreversed=true, title="Approximated")]
Mmat = [M(waveNumber, rMid[_i], rPos[_j], rPos[_j+1]) for _i in 1:N-1, _j in 1:N-1]

T = -inv(Mmat)

T_orig = -inv(BoundaryWall.calcTMatrix(waveVector, -0.25im, rij, ds, diagind(size(rij)...)))

heatmap!(ax[1], abs2.(T),colorscale=log10, colormap=:dense)
heatmap!(ax[2], abs2.(T_orig), colorscale=log10, colormap=:dense)
[ax.aspect=DataAspect() for ax in ax]
fig
end


function boundaryWallWave(
  waveVector::SVector{2,Float64},
  xMid::Vector,
  yMid::Vector,
  xDomain::Vector, 
  yDomain::Vector,
  σ::ComplexF64,
  Mij::Matrix,
  arcLengths::Vector,
  numSegments::Int,
  )

# calculates boundary wall method
waveNumber = norm(waveVector)
waveAtBoundary = BoundaryWall.incidentWave.(Ref(waveVector), xMid, yMid);
waveAtDomain   = BoundaryWall.incidentWave.(Ref(waveVector), xDomain, yDomain)

# calculate M matrix
# Mij          = @. σ * calcGreenFun(waveNumber, rijMid) * arcLengths[1]#  arcLengths # multiplies cols by ds[j]
# Mij[dindex] .= @. σ * calcGreenFun(waveNumber, arcLengths/2) * arcLengths[1] # arcLengths;

# T          = potentialFunction .* inv(I(numSegments) .- potentialFunction .* Mij);
# TPHI = T * waveAtBoundary;
# TPHI         = -(Mij \ waveAtBoundary)
# TPHI = -IterativeSolvers.gmres(Mij, waveAtBoundary)

# TPHI = -solve(LinearProblem(Mij, waveAtBoundary),IterativeSolversJL_GMRES)
TPHI = -Mij\ waveAtBoundary

# potentialContribution = mapreduce( (j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1] , +, 1:numSegments)
# wave =  waveAtDomain + ThreadsX.mapreduce((j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1],+, 1:numSegments)
wave = waveAtDomain
@inbounds for j in 1:numSegments
  RR = @.  hypot(xDomain  - (@view xMid[j]), yDomain - (@view yMid[j]))
  wave += @. σ * BoundaryWall.calcGreenFun(waveNumber,RR) * (@view arcLengths[j]) * (@view TPHI[j])
end

# wave =  + potentialContribution; # Huygens-like principle
# return -inv(Mij)
# T = -inv(Mij)
# return T[id]
return wave
end

xdom = LinRange(-2,2,100)
ydom = LinRange(-2,2,100)

X = repeat(xdom, inner = length(ydom))
Y = repeat(ydom, inner = length(xdom))

wave = boundaryWallWave(waveVector, xm, ym, X,Y,-0.25im, Mmat, ds, length(ds))

heatmap(abs2.(reshape(wave, 100, 100)))