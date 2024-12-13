using GLMakie
using SpecialFunctions
using FunctionZeros
using StaticArrays
using FunctionZeros
using Meshes
using GSL



sf_hankelh1_n(n::Any, x::Any) = GSL.sf_bessel_Jn(n, x) + im * GSL.sf_bessel_Yn(n,x)


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
  
#   ThreadsX.foreach(-nearlyDegenerateM:nearlyDegenerateM) do _K
    for _K in -nearlyDegenerateM:nearlyDegenerateM
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
    _top = 2π*γ*σ*R*GSL.sf_bessel_Jn(n,k*R)*sf_hankelh1_n(n,k*R)
    return _top/(1-_top)
end

# z_p,m, ω_±p = 0
"""
Omega for infinite potential
# Arguments
`_m`: order of Bessel root
`_p`: eigenstate index
"""
function omegaN(_m::Int64,_p::Int64)
    if _m == _p
        return 0.0
    else
        return -1.0
    end
end

# function uN(n::Int64, p::Int64, γ::Float64, σ::ComplexF64,k::Float64, R::Float64)
#     # return omegaN(n, γ,σ,k,R)*besselj(n,k*R)/hankelh1(n,k*R)
#     if isinf(γ)
#         return omegaN(n, p)*GSL.sf_bessel_Jn(n,k*R)/sf_hankelh1_n(n,k*R)
#     else
#         return omegaN(n, γ,σ,k,R)*GSL.sf_bessel_Jn(n,k*R)/sf_hankelh1_n(n,k*R)
#     end
# end

function uN(n::Int64, p::Int64, k::Float64, R::Float64)
    return omegaN(n, p)*GSL.sf_bessel_Jn(n,k*R)/sf_hankelh1_n(n,k*R)
end

# valFun(_k::Float64, _r::Float64,_ω::ComplexF64) = besselj(0,_k*_r)*(1+_ω)

"""
    Circular billiard (gamma = Inf)
"""
function circleWave(k::Union{Float64, BigFloat}, 
                    γ::Union{Float64, BigFloat}, 
                    σ::Union{ComplexF64,Complex{BigFloat}}, 
                    R::Union{Float64, BigFloat}, 
                    r::Union{SVector{2,Float64}, SVector{2,BigFloat}},
                    p::Int64,
                    # r::Vector{SVector{2,Float64}},
                    α::Union{Float64,BigFloat}, 
                    _N_sum::Int64)

    # _params = [SVector(im^n,omegaN(n,p,γ,σ,k,R),cos(n*(r[2]+(-1)^n*α)),uN(n,γ,σ,k,R)) for n in 1:_N_sum]
    _params = [SVector(im^n,omegaN(n,p),cos(n*(r[2]+(-1)^n*α)),uN(n,p,k,R)) for n in 1:_N_sum]
    _Jn = GSL.sf_bessel_Jn_array(1,_N_sum,k*r[1])
    
    if r[1]<R
        # _val= GSL.sf_bessel_J0(k*r[1])*(1+omegaN(0,γ,σ,k,R))
        _val= GSL.sf_bessel_J0(k*r[1])*(1+omegaN(0,p))

        _sum = sum(@. getindex(_params,1)*_Jn*(1+getindex(_params,2))*getindex(_params,3))

        return _val + 2*_sum
    end
    
    # # else
    _val = GSL.sf_bessel_J0(k*r[1]) + uN(0,p,k,R) * (GSL.sf_bessel_J0(k*r[1])+im*GSL.sf_bessel_Y0(k*r[1]))
    _Yn  = GSL.sf_bessel_Yn_array(1,_N_sum,k*r[1])
    _sum = sum(
        @. getindex(_params,1) * (_Jn + getindex(_params,4) * (_Jn + im*_Yn)) * getindex(_params,3)
        )
    return _val + 2*_sum   
    
end


function circleDaLuz(
    k::Union{Float64,BigFloat},
    R::Union{Float64, BigFloat},
    r::Union{SVector{2,Float64}, SVector{2,BigFloat}},
    n::Int64
    )
    if r[1]<R
        # not assuming k_m,n other than eigenvalues
        # there is always tunneling
        return GSL.sf_bessel_Jn(n,k*r[1])*exp(im*n*r[2])
    else    

    _Yn  = GSL.sf_bessel_Yn(n,k*R)
    _Jn  = GSL.sf_bessel_Jn(n,k*R)
    return im * exp(im*n*r[2]) * (GSL.sf_bessel_Jn(n,k*r[1]) - _Jn * (GSL.sf_bessel_Jn(n,k*r[1])+GSL.sf_bessel_Yn(n,k*r[1]))/(_Jn + im*_Yn))
    end

end




import BoundaryWall as BWM

begin # defining boundary
SIGMA        = -0.25im
R            = 1.0
N            = 500
TH           = LinRange(-pi,pi, N)
x,y,xm,ym,ds = BWM.createCircle(R, TH, SVector(0.0, 0.0))
rij          = BWM.calcDistances(xm, ym)

incidence    = 0.0
end

begin # domain
setrounding(BigFloat, RoundDown)
setprecision(BigFloat, 256)  # Precision is set in bits
nx,ny = 450,450
x0, xf = (-1.0R, 1.0R)
y0, yf = (-1.0R, 1.0R)
xdom = LinRange(x0, xf, nx)
ydom = LinRange(y0, yf, ny)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coords.(vertices(MESH)))
XDOM = BigFloat.(first.(COORDS))
YDOM = BigFloat.(last.(COORDS))

car2pol(x::Float64, y::Float64) = SVector(hypot(x,y), atan(y,x))
car2pol(x::BigFloat, y::BigFloat) = SVector(hypot(x,y), atan(y,x))
polar_coords = [car2pol(_r...) for _r in COORDS]
# polar_coords = [car2pol(_r...) for _r in SVector.(XDOM, YDOM)]

end


using .Threads
begin
GC.gc()
nearlyDegenerateM = 2
quantumNumberM0   = 100
quantumNumberN0   = 50
turningPointsQ    = 5
windingsP         = 2
phi_0 = 0.0
nSum = 50
wave = zeros(length(polar_coords))
wave_maioli = zeros(length(polar_coords))

# using ThreadsX
# ThreadsX.foreach(-nearlyDegenerateM:nearlyDegenerateM) do _K
@inbounds for _K in -nearlyDegenerateM:nearlyDegenerateM
_m_prime      = quantumNumberM0 + turningPointsQ*_K
_m_prime
_n_prime      = quantumNumberN0 - windingsP*_K

_bessel_roots = _k = FunctionZeros.besselj_zero(_m_prime, _n_prime)/R #

println("Wave number: $_K, \t $_bessel_roots  ")

_wave = @inbounds [(besselj(_m_prime,_k*_r[1])*exp(im*_m_prime*_r[2])/(sqrt(π)*R*besselj(_m_prime-1,_k*R)))*(_r[1]<R) for _r in polar_coords]
_wave = _wave + _wave .* exp.(-im * _m_prime * last.(polar_coords))

global wave += exp(im*_K*turningPointsQ*phi_0)*_wave*1/sqrt(2*nearlyDegenerateM+1)

# _wave = @inbounds [(circleWave(_k, Inf,SIGMA,R,_r,_m_prime, 1.0pi, nSum)) *(_r[1]<R) for _r in polar_coords]

# global wave_maioli += exp(im*_K*turningPointsQ*phi_0)*(_wave)*1/sqrt(2*nearlyDegenerateM+1)


# global wave += GSL.sf_bessel_Jn.(_m_prime,_k*_r).*exp.(-im*_m_prime*_θ)/(sqrt(π)*R*GSL.sf_bessel_Jn(_m_prime-1,_k*R)).*(_r.<R)

println("Done!")
end
end
wave_numeric = boundaryWallDegenerate(nearlyDegenerateM, quantumNumberM0, quantumNumberN0, turningPointsQ, windingsP, 0.0, BWM.planeWave, x,y,xm, ym,first.(COORDS), last.(COORDS), SIGMA, ds, rij, length(ds), N, 2, Inf)

_m, _n = (1,1)
_k = FunctionZeros.besselj_zero(_m,_n)
wave_daluz = [circleDaLuz(_k, R, _r, _m) for _r in polar_coords]
wave_maioli = [(circleWave(_k/R, Inf,SIGMA,R,_r,_m, 1.0pi, nSum)abs(circleWave(_k/R, Inf,SIGMA,R,_r,_m, 1.0pi, nSum))*exp(im*_m*_r[2]))*(_r[1]<R) for _r in polar_coords]

# wave_scars  = [GSL.sf_bessel_Jn(_m,_k*_r[1])*exp(im*_m*_r[2])/(sqrt(π)*R*besselj(_m-1,_k*R))*(_r[1]<R) for _r in polar_coords]
wave_scars  = [exp(-im*pi/2)*(GSL.sf_bessel_Jn(_m,_k*_r[1])*exp(im*_m*_r[2]))*(_r[1]<R) for _r in polar_coords]

let fig = Figure()
    ax = [Axis(fig[1,1])]
    # heatmap!(ax[1], xdom, ydom, reshape(real.(wave_scars), nx, ny), colormap=:linear_kbgyw_5_98_c62_n256)
    # heatmap!(ax[2], xdom, ydom, reshape(real.(wave_maioli), nx, ny), colormap=:linear_kbgyw_5_98_c62_n256)

    # heatmap!(ax[1], xdom, ydom, reshape(abs2.(wave_maioli), nx, ny))
    heatmap!(ax[1], xdom, ydom, reshape(abs2.(wave), nx, ny), colormap=:linear_kbgyw_5_98_c62_n256, interpolate=true)
    [lines!(ax, x, y) for ax in ax]
    [ax.aspect=DataAspect() for ax in ax]
    fig
end
