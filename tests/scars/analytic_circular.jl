using GLMakie
using SpecialFunctions
using FunctionZeros
using StaticArrays
using FunctionZeros
using Meshes

import BoundaryWall as BWM


function boundaryWallDegenerate(
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

function omegaN(n::Int64, γ::Float64, σ::ComplexF64,k::Float64, R::Float64)
    _top = 2π*γ*σ*R*besselj(n,k*R)*hankelh1(n,k*R)
    return _top/(1-_top)
end

function uN(n::Int64, γ::Float64, σ::ComplexF64,k::Float64, R::Float64)
    return omegaN(n, γ,σ,k,R)*besselj(n,k*R)/hankelh1(n,k*R)
end

# valFun(_k::Float64, _r::Float64,_ω::ComplexF64) = besselj(0,_k*_r)*(1+_ω)

function circleWave(k::Float64, γ::Float64, σ::ComplexF64, R::Float64, r::SVector{2,Float64},α::Float64, _N_sum::Int64)

    if r[1]<R
        _val= besselj(0,k*r[1])*(1+omegaN(0,γ,σ,k,R))
        _sum = sum([im^n*besselj(n,k*r[1])*(1+omegaN(n,γ,σ,k,R))*cos(n*(r[2]+(-1)^n*α)) for n in 1:_N_sum])
        return _val + 2*_sum
    end
    
    # else
    _val = besselj0(k*r[1]) + uN(0,γ,σ,k,R) * hankelh1.(0,k*r[1])
    _sum = sum([im^n*(besselj.(n,k*r[1]) + uN(n,γ,σ,k,R)*hankelh1(n,k*r[1]))*cos(n*(r[2]+(-1)^n*α)) for n in 1:_N_sum])
    return _val + 2*_sum
end


import BoundaryWall as BWM

begin # defining boundary
SIGMA        = -0.25im
R            = 1.0
N            = 200
TH           = LinRange(-pi,pi, N)
x,y,xm,ym,ds = BWM.createCircle(R, TH, SVector(0.0, 0.0))
rij          = BWM.calcDistances(xm, ym)

incidence    = 0.0
end

begin # domain
nx,ny = 100,100
x0, xf = (-2R, 2R)
y0, yf = (-2R, 2R)
xdom = LinRange(x0, xf, nx)
ydom = LinRange(y0, yf, ny)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coordinates.(vertices(MESH)))

car2pol(x::Float64, y::Float64) = SVector(hypot(x,y), atan(y,x))
polar_coords = [car2pol(_r...) for _r in COORDS]

end





begin
nearlyDegenerateM = 2
quantumNumberM0   = 25
quantumNumberN0   = 1
turningPointsQ    = 2
windingsP         = 1
phi_0 = 0.0
wave = zeros(length(polar_coords))
nSum = 50

# using ThreadsX
ThreadsX.foreach(-nearlyDegenerateM:nearlyDegenerateM) do _K
_m_prime      = quantumNumberM0 + turningPointsQ*_K
_n_prime      = quantumNumberN0 - windingsP*_K

_bessel_roots = FunctionZeros.besselj_zero(_m_prime, _n_prime) #

println("Wave number: $_K")
_wave = @inbounds [circleWave(_bessel_roots/R, 10000.0,SIGMA,R,_r, 1.0pi, nSum) for _r in polar_coords]
global wave += @inbounds exp(im*_K*turningPointsQ*phi_0) * _wave * 1/sqrt(2*nearlyDegenerateM+1)
println("Done!")
end

println("calculating BWM")

wave_numeric = boundaryWallDegenerate(nearlyDegenerateM, quantumNumberM0, quantumNumberN0, turningPointsQ, windingsP, 0.0, BWM.planeWave, x,y,xm, ym,first.(COORDS), last.(COORDS), SIGMA, ds, rij, length(ds), N, 2, Inf)

end




let fig = Figure()
    ax = [Axis(fig[1,1]), Axis(fig[1,2])]
    # viz!(ax, MESH, color=angle.(wave),colormap=:twilight)
    viz!(ax[1], MESH, color=abs.(wave),colormap=:turbo);viz!(ax[2], MESH, color=abs.(sum(wave_numeric)), colormap=:turbo)
    # viz!(ax[1], MESH, color=angle.(wave),colormap=:twilight);viz!(ax[2], MESH, color=angle.(sum(wave_numeric)), colormap=:twilight)
    [lines!(ax, x, y) for ax in ax]
    [ax.aspect=DataAspect() for ax in ax]
    fig
end