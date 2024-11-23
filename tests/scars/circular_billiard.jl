using GLMakie
using BoundaryWall
using StaticArrays
using Meshes
using FunctionZeros
using ThreadsX

using CircularArrays
import BoundaryWall as BWM

begin # defining boundary
SIGMA        = -0.25im
R            = 1.0
N            = 5000
TH           = LinRange(-pi,pi, N)
x,y,xm,ym,ds = BWM.createCircle(R, TH, SVector(0.0, 0.0))
rij          = BWM.calcDistances(xm, ym)

incidence    = 0.0
end

rMid = SVector.(xm,ym)
rPos = CircularArray(SVector.(x,y))


begin # domain
nx,ny = 200,200
x0, xf = (-2R, 2R)
y0, yf = (-2R, 2R)
xdom = LinRange(x0, xf, nx)
ydom = LinRange(y0, yf, ny)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coords.(vertices(MESH)))

XDOM, YDOM = first.(COORDS), last.(COORDS)
end


rDom = SVector.(XDOM, YDOM)

function buildDegenerate(
  nearlyDegenerateM::Int64,
  quantumNumberM0::Int64,
  quantumNumberN0::Int64,
  turningPointsQ::Int64,
  windingsP::Int64,
  phi_0::Float64,
  waveFun::Function,
  xBoundary::Vector,
  yBoundary::Vector,
  xMid::Vector,
  yMid::Vector,
  xDomain::Vector, 
  yDomain::Vector,
  σ::ComplexF64,
  arcLengths::Vector,
  rijMid::Matrix,
  numSegments::Int,
  numSegmentsMod::Int,
  band::Int,
  potentialStrength::Union{Float64, ComplexF64, Vector})

# nearlyDegenerateM = 2
# quantumNumberM0   = 10
# quantumNumberN0   = 2
GC.gc()

wave = Vector{ComplexF64}[]#  zeros(ComplexF64, length(xDomain))

ThreadsX.foreach(-nearlyDegenerateM:nearlyDegenerateM) do _K
  _m_prime      = quantumNumberM0 + turningPointsQ*_K
  _n_prime      = quantumNumberN0 - windingsP*_K

  _bessel_roots = FunctionZeros.besselj_zero(_m_prime, _n_prime) # find roots
  
  _wave_number  = _bessel_roots/R
  println("k_{$_m_prime, $_n_prime}/R: $_wave_number")
  push!(wave, exp(im*_K*turningPointsQ*phi_0)*BWM.boundaryWallWave(_wave_number*SVector(1.0, 0.0), waveFun,
           xBoundary, yBoundary, xMid, yMid, xDomain, yDomain, σ, arcLengths, rijMid, numSegments, numSegmentsMod, band, potentialStrength))

end
return wave
end





# scarring
