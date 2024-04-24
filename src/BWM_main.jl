"""
Developed by Alberto RB for use in Billiard problems.
If you use this, thank me at least!!

Main functions used for the Boundary Wall method.
"""

module BWM_main
using StaticArrays, LinearAlgebra
using SpecialFunctions
using ThreadsX
using IterativeSolvers
using LinearSolve
# using HYPRE
using QuadGK
using FLoops
using LoopVectorization
using CircularArrays
using HCubature

using Interpolations

export boundaryWallWave, boundaryWallVec, planeWave, gaussianWave, shapedWave, gradient

global eGamma = 0.57721566490153286060

"""
  `segmentIterator(v)`

Returns a modulo type filter for iterating each segment. Needed when dealing with
line integrated (quadgk routine) methods for calculating M matrix for discontinuous 
curves. This way, the method does not take into account the vector position between
disjointed curves as a line integral.

# Arguments:
- `v::CircularVector`: Array with vector positions of each vertex
- `numSegments::Int64`: Number of segments per element (not all disjoint elements)
"""
segmentIterator(v::CircularVector,numSegments::Int64) = filter( x -> mod(x,numSegments+1)!=0, eachindex(v))


 # for gaussian wave
Base.atan(v::SVector{2,Float64}) = atan(v[2],v[1])


"""simple plane wave"""
planeWave(k::SVector{2, Float64}, r::SVector{2, Float64}) = exp(-1im * dot(k,r))

"""Willis et al (2008). Amplified total internal reflection: theory, analysis, and demonstration of existence via FDTD. Opt Lett.

Be mindful of the convergence, as some expansions might not converge fast enough (e.g. very long gaussian pulses)

"""
gaussianWave(k::SVector{2, Float64}, r::SVector{2, Float64}, ω::Float64; abstol=0) = 1/2pi * first(hquadrature(t -> pi * ω^2 * exp(-(ω*(t-atan(k))/2)^2) * planeWave(norm(k)*SVector(cos(t),sin(t)),r), -pi/2+atan(k), pi/2+atan(k); atol=abstol))

"""Vaishnav et al. (2007). Matter-wave scattering and guiding by atomic arrays. PRA, doi:10.1103/physreva.76.013620"""
shapedWave(k::SVector{2, Float64}, r::SVector{2, Float64}, ω::Float64; abstol=0) = 1/2pi * first(hquadrature(t -> exp(-((t-atan(k))/ω)^2) * planeWave(norm(k)*SVector(cos(t),sin(t)),r), -pi/2+atan(k), pi/2+atan(k); atol=abstol))


function heavyside(x::Float64, xp::Float64)
  return sign(x-xp)
end

function window(x::Float64,a::Float64,b::Float64)
  return heavyside(x, a) - heavyside(x,b)
end

"""
  `calcMmatrix(waveNumber,σ,rijMid,arcLengths,dindex)`

Computes the first order discretization of an M-matrix used 
for the Boundary Wall Method. (see MGE da Luz, 1997).

Useful for testing linear solve schemes and comparing band 
integration methods.

# Arguments:
- `waveNumber::Float64`: self explanatory
- `σ::ComplexF64`: Free particle Green's function parameter 
- `rijMid::Matrix`: distance matrix between segments
- `arcLengths::Vector`: self explanatory
- `dindex::StepRange`: indexes for diagonal, supplied beforehand
"""
function calcMmatrix(waveNumber::Float64,
                     σ::ComplexF64,
                     rijMid::Matrix,
                     arcLengths::Vector,
                     dindex::StepRange)

  Mij          = @. σ * BoundaryWall.calcGreenFun(waveNumber, rijMid) * arcLengths[1] #  arcLengths # multiplies cols by ds[j]
  Mij[dindex] .= @. σ * BoundaryWall.calcGreenFun(waveNumber, arcLengths/2) * arcLengths[1] # arcLengths;
  return Mij
end



"""
  `calcTPHI(Mij, waveAtBoundary, gamma::Union{Float64, ComplexF64})`

Calculates TΦ depending for an infinite or finite potential (constant though).
# Arguments: 
"""
function calcTPHI(Mij::Matrix{ComplexF64}, waveAtBoundary::Vector{ComplexF64}, gamma::Union{Float64, ComplexF64})
  # -solve(LinearProblem(Mij, waveAtBoundary),KrylovJL_GMRES())
  if (Base.isinf(gamma))
    # return -solve(Mij, waveAtBoundary)
    println("Using infinite approximation...")
    return -Mij\waveAtBoundary
  end
  return gamma *  solve(LinearProblem(LinearAlgebra.I(size(Mij, 1)) .- gamma * Mij, waveAtBoundary))
end

# could, in principle, replace it with the method above
function calcTPHI(Mij::Matrix{ComplexF64}, waveAtBoundary::Vector{ComplexF64}, gamma::Vector)
  # -solve(LinearProblem(Mij, waveAtBoundary),KrylovJL_GMRES())
 
  return gamma .*  solve(LinearProblem(LinearAlgebra.I(size(Mij, 1)) .- gamma .* Mij, waveAtBoundary))
end


"""
  `calcMmatrix(waveNumber,σ,rijMid,arcLengths,dindex)`

Computes the first order discretization of an M-matrix used 
for the Boundary Wall Method. (see MGE da Luz, 1997).

Useful for testing linear solve schemes and comparing band 
integration methods.

# Arguments:
- `waveVector::Float64`: self explanatory
- `σ::ComplexF64`: Free particle Green's function parameter 
- `ri::SVector{2,Float64}`: observation point (midpoint)
- `r0::SVector{2,Float64}`: lower endpoint of curve segment
- `r1::SVector{2,Float64}`: upper endpoint of curve segment
"""
function calcMmatrix(waveNumber::Float64, σ::ComplexF64, ri::SVector{2,Float64}, r0::SVector{2,Float64}, r1::SVector{2, Float64})
  # G0(k::Float64, ri::SVector{2, Float64}, rj::SVector{2,Float64}) = BoundaryWall.calcGreenFun(k, norm(ri-rj))
  # M(k::Float64, ri::SVector{2,Float64}, r0::SVector{2,Float64}, r1::SVector{2, Float64}) = first(quadgk(t -> G0(k, ri, r(r0, r1, t)), 0.0, 0.5, 1.0)) * norm(r1 - r0)  # scaling
  
  # Integrate( G₀(k, r(t)) * dt, t∈(0, 1))
  return σ * first(quadgk(t -> calcGreenFun(waveNumber, ri, r(r0, r1, t)), 0.0, 0.5, 1.0)) * norm(r1 - r0)  # scaling
end

function calcMmatrix(type::Symbol, θ::Float64, waveNumber::Float64, σ::ComplexF64, ri::SVector{2,Float64}, r0::SVector{2,Float64}, r1::SVector{2, Float64})
  # G0(k::Float64, ri::SVector{2, Float64}, rj::SVector{2,Float64}) = BoundaryWall.calcGreenFun(k, norm(ri-rj))
  # M(k::Float64, ri::SVector{2,Float64}, r0::SVector{2,Float64}, r1::SVector{2, Float64}) = first(quadgk(t -> G0(k, ri, r(r0, r1, t)), 0.0, 0.5, 1.0)) * norm(r1 - r0)  # scaling
  
  # Integrate( G₀(k, r(t)) * dt, t∈(0, 1))
  return σ * first(quadgk(t -> calcGreenFun(type, θ, waveNumber, norm(ri - r(r0, r1, t))), 0.0, 0.5, 1.0)) * norm(r1 - r0)  # scaling
end

"""
  `calcMmatrix(waveNumber,σ,rijMid,arcLengths,dindex)`

Band integrated M-matrix (see MGE da Luz, 1997). Banded M matrix.

# Arguments:
- `waveNumber::Float64`: self explanatory
- `σ::ComplexF64`: Free particle Green's function parameter 
- `rijMid::Matrix`: distance matrix between segments
- `arcLengths::Vector`: self explanatory
- `dindex::StepRange`: indexes for diagonal, supplied beforehand
- `band::Int64`,
"""
function calcMmatrix(waveNumber::Float64,
                     σ::ComplexF64,
                     rijMid::Matrix,
                     arcLengths::Vector,
                     band::Int64,
                     _rm::Vector{SVector{2, Float64}}, # rMid in tests
                     _r::CircularVector{SVector{2, Float64}, Vector{SVector{2, Float64}}},
                     numSegments) # rPos in tests
  
  @assert ((band > 0) && band < size(rijMid,1))
  Mij          = @. σ * BoundaryWall.calcGreenFun(waveNumber, rijMid) * arcLengths[1] #  arcLengths # multiplies cols by ds[j]
  [Mij[i,j] = calcMmatrix(waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Mij,1), (j,_j) in zip(axes(Mij, 2), segmentIterator(_r, numSegments)) if abs(i - j) < band]
  
  return Mij
end


# abstract type IncidentWave end;
# struct


"""
  `incidentWave(k, x, y)`

Returns a plane wave (or any other wave) travelling in direction `k/norm(k)`, 
calculated at a point (x,y).

# Arguments:
- `k::SVector`: wave Vector
- `x::Float64`: x coord
- `y::Float64`: y coord
"""
function incidentWave(k::SVector, x::Float64, y::Float64)
  @warn "Deprecated in favour of planeWave"
  return exp(- 1im  *  (k[1] * x + k[2] * y ));
        #  2*exp(- 1im  *  (k[1] * x + 0.25*k[1] * y)) +
        #  2*exp(- 1im  *  (k[1] * x - 0.25*k[1] * y));  
  # return exp(-im * norm(k)*hypot(x, y - 1.6165807537309522))#/hypot(x, y - 25)
  # return gaussian(norm(k), x, y)
end



"""
  `boundaryWallWave(waveVector,xBoundary,yBoundary,xMid,yMid,xDomain,yDomain,σ,arcLengths,numSegments,numSegmentsMod,band,potentialStrength)`

Computes the scattered field by a wall using a band integrated (quadgk) T-matrix.

# Arguments:
- `waveVector::SVector{2,Float64}`: incident wave vector
- `xBoundary::Vector`: boundary points
- `yBoundary::Vector`: boundary points
- `xMid::Vector`: boundary midpoints
- `yMid::Vector`: boundary midpoints
- `xDomain::Vector,`: sample points
- `yDomain::Vector`: sample points
- `σ::ComplexF64`: green's function constant
- `arcLengths::Vector`: self explanatory
- `rijMid::Matrix`: distance matrix
- `numSegments::Int`: number of **all** segments
- `numSegmentsMod::Int`: number of segments per element (assumes homogeneity)
- `band::Int`: diagonal band length 
- `potentialStrength::Union{Float64, ComplexF64, Vector}`: γ or k²Δϵ, depending on the kernel
"""
function boundaryWallWave(
  waveVector::SVector{2,Float64},
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
  potentialStrength::Union{Float64, ComplexF64, Vector}
  )

# calculates boundary wall method
waveNumber = norm(waveVector)
# waveAtBoundary = incidentWave.(Ref(waveVector), xMid, yMid);
# waveAtDomain   = incidentWave.(Ref(waveVector), xDomain, yDomain)

rMid = SVector.(xMid, yMid)
rDom = SVector.(xDomain, yDomain)
rPos = CircularArray(SVector.(xBoundary, yBoundary)[1:end])

waveAtBoundary = waveFun.(Ref(waveVector), rMid)
waveAtDomain   = waveFun.(Ref(waveVector), rDom)

# midpoints and circular position vector


Mij = calcMmatrix(waveNumber, σ, rijMid, arcLengths, band, rMid, rPos, numSegmentsMod)

TPHI = calcTPHI(Mij, waveAtBoundary, potentialStrength)

# potentialContribution = mapreduce( (j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1] , +, 1:numSegments)
# wave =  waveAtDomain + ThreadsX.mapreduce((j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1],+, 1:numSegments)
wave = waveAtDomain
@inbounds for j in 1:numSegments
  RR = @.  hypot(xDomain  - (@view xMid[j]), yDomain - (@view yMid[j]))
  wave += @. σ * calcGreenFun(waveNumber,RR) * (@view arcLengths[j]) * (@view TPHI[j])
end

return wave
end

"""
  `boundaryWallWave(waveVector,xBoundary,yBoundary,xMid,yMid,xDomain,yDomain,σ,arcLengths,numSegments,numSegmentsMod,potentialStrength)`

Computes the scattered field by a wall using a fully integrated (quadgk) T-matrix.

# Arguments:
- `waveVector::SVector{2,Float64}`: incident wave vector
- `xBoundary::Vector`: boundary points
- `yBoundary::Vector`: boundary points
- `xMid::Vector`: boundary midpoints
- `yMid::Vector`: boundary midpoints
- `xDomain::Vector,`: sample points
- `yDomain::Vector`: sample points
- `σ::ComplexF64`: green's function constant
- `arcLengths::Vector`: self explanatory
- `numSegments::Int`: number of **all** segments
- `numSegmentsMod::Int`: number of segments per element (assumes homogeneity)
- `potentialStrength::Union{Float64, ComplexF64, Vector}`: γ or k²Δϵ, depending on the kernel
"""
function boundaryWallWave(
  waveVector::SVector{2,Float64},
  waveFun::Function,
  xBoundary::Vector,
  yBoundary::Vector,
  xMid::Vector,
  yMid::Vector,
  xDomain::Vector, 
  yDomain::Vector,
  σ::ComplexF64,
  arcLengths::Vector,
  numSegments::Int,
  numSegmentsMod::Int,
  potentialStrength::Union{Float64, ComplexF64, Vector}
  )

# calculates boundary wall method
waveNumber = norm(waveVector)


# calculate M matrix
# Mij          = @. σ * calcGreenFun(waveNumber, rijMid) * arcLengths[1]#  arcLengths # multiplies cols by ds[j]
# Mij[dindex] .= @. σ * calcGreenFun(waveNumber, arcLengths/2) * arcLengths[1] # arcLengths;

# midpoints and circular position vector
rMid = SVector.(xMid, yMid)
rDom = SVector.(xDomain, yDomain)
rPos = CircularArray(SVector.(xBoundary, yBoundary)[1:end])

waveAtBoundary = waveFun.(Ref(waveVector), rMi1d);
waveAtDomain   = waveFun.(Ref(waveVector), rDom)


# Mij = calcMmatrix(waveNumber, σ, rijMid, arcLengths, 3, rMid, rPos)
# Mij = calcMmatrix(waveNumber, σ, rijMid, arcLengths, band, rMid, rPos, numSegmentsMod)
Mij = [calcMmatrix(waveNumber, σ, rMid[i], rPos[j], rPos[j+1]) for i in axes(rMid,1), j in segmentIterator(rPos, numSegmentsMod)]

# T          = potentialFunction .* inv(I(numSegments) .- potentialFunction .* Mij);
# TPHI = T * waveAtBoundary;
# TPHI         = -(Mij \ waveAtBoundary)
# TPHI = -IterativeSolvers.gmres(Mij, waveAtBoundary)

# TPHI = -solve(LinearProblem(Mij, waveAtBoundary),IterativeSolversJL_GMRES)
# TPHI = -solve(LinearProblem(Mij, waveAtBoundary),KrylovJL_GMRES())
TPHI = calcTPHI(Mij, waveAtBoundary, potentialStrength)

# potentialContribution = mapreduce( (j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1] , +, 1:numSegments)
# wave =  waveAtDomain + ThreadsX.mapreduce((j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1],+, 1:numSegments)
wave = waveAtDomain
@inbounds for j in 1:numSegments
  RR = @.  hypot(xDomain  - (@view xMid[j]), yDomain - (@view yMid[j]))
  wave += @. σ * calcGreenFun(waveNumber,RR) * (@view arcLengths[j]) * (@view TPHI[j])
end

# wave =  + potentialContribution; # Huygens-like principle
# return -inv(Mij)
# T = -inv(Mij)
# return T[id]
return wave
end


"""
  `boundaryWallWave(waveVector,xMid, yMid, xDomain, yDomain, σ, rijMid, arcLengths, numSegments, dindex)`

Computes the scattered field by a wall (aka the Boundary Wall Method), assuming
an infinite potential strength.

# Arguments:
- `waveVector::SVector{2,Float64}`: self explanatory
- `xMid::Vector`: x coordinate of segment's midpoints
- `yMid::Vector`: y coordinate of segment's midpoints
- `xDomain::Vector,`: x coords for domain of evaluation
- `yDomain::Vector`: y coords for domain of evaluation
- `σ::ComplexF64`: free particle Green's function parameter
- `rijMid::Matrix`: distance matrix between segments
- `arcLengths::Vector`: length of segments
- `numSegments::Int`: number of segments
- `dindex::StepRange`: diagonal index supplied beforehand
"""
function boundaryWallWave(
  waveVector::SVector{2,Float64},
  xMid::Vector,
  yMid::Vector,
  xDomain::Vector, 
  yDomain::Vector,
  σ::ComplexF64,
  rijMid::Matrix,
  arcLengths::Vector,
  numSegments::Int,
  dindex::StepRange,
  )

@warn "This method is deprecated in favour of newer ones that allow band integrated matrices."

# calculates boundary wall method
waveNumber = norm(waveVector)
waveAtBoundary = incidentWave.(Ref(waveVector), xMid, yMid);
waveAtDomain   = incidentWave.(Ref(waveVector), xDomain, yDomain)

# calculate M matrix
# Mij          = @. σ * calcGreenFun(waveNumber, rijMid) * arcLengths[1]#  arcLengths # multiplies cols by ds[j]
# Mij[dindex] .= @. σ * calcGreenFun(waveNumber, arcLengths/2) * arcLengths[1] # arcLengths;

Mij = calcMmatrix(waveNumber, σ, rijMid, arcLengths, dindex)

# T          = potentialFunction .* inv(I(numSegments) .- potentialFunction .* Mij);
# TPHI = T * waveAtBoundary;
# TPHI         = -(Mij \ waveAtBoundary)
# TPHI = -IterativeSolvers.gmres(Mij, waveAtBoundary)

# TPHI = -solve(LinearProblem(Mij, waveAtBoundary),IterativeSolversJL_GMRES)
TPHI = -solve(LinearProblem(Mij, waveAtBoundary),KrylovJL_GMRES())

# potentialContribution = mapreduce( (j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1] , +, 1:numSegments)
# wave =  waveAtDomain + ThreadsX.mapreduce((j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1],+, 1:numSegments)
wave = waveAtDomain
@inbounds for j in 1:numSegments
  RR = @.  hypot(xDomain  - (@view xMid[j]), yDomain - (@view yMid[j]))
  wave += @. σ * calcGreenFun(waveNumber,RR) * (@view arcLengths[j]) * (@view TPHI[j])
end

# wave =  + potentialContribution; # Huygens-like principle
# return -inv(Mij)
# T = -inv(Mij)
# return T[id]
return wave
end

"""
  `boundaryWallWave(..., potentialFunction::Float64)`

Computes scattered field by a wall at domain points with a functional potential strength.
"""
function boundaryWallWave(
  waveVector::SVector{2,Float64},
  xMid::Vector,
  yMid::Vector,
  xDomain::Vector, 
  yDomain::Vector,
  σ::ComplexF64,
  rijMid::Matrix,
  arcLengths::Vector,
  numSegments::Int,
  dindex::StepRange,
  potentialFunction::Vector,
  )

@warn "This method is deprecated in favour of newer ones that allow band integrated matrices."
# calculates boundary wall method
waveNumber = norm(waveVector)
waveAtBoundary = incidentWave.(Ref(waveVector), xMid, yMid);
waveAtDomain   = incidentWave.(Ref(waveVector), xDomain, yDomain)

# calculate M matrix
Mij          = @. σ * calcGreenFun(waveNumber, rijMid) * arcLengths[1] # arcLengths # multiplies cols by ds[j]
Mij[dindex] .= @. σ * calcGreenFun(waveNumber, arcLengths/2) * arcLengths[1]# arcLengths;

TPHI         = inv(LinearAlgebra.I(numSegments) .- potentialFunction .* Mij)
TPHI         = potentialFunction .* TPHI * waveAtBoundary
# TPHI         = T * waveAtBoundary;
# TPHI         = -(Mij \ waveAtBoundary)

# potentialContribution = mapreduce( (j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1] , +, 1:numSegments)
# potentialContribution = ThreadsX.mapreduce((j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1],+, 1:numSegments)
potentialContribution = zeros(ComplexF64, length(xDomain));
@inbounds for j in 1:numSegments
  potentialContribution += @. σ * calcGreenFun(waveNumber, hypot(xDomain  - (@view xMid[j]), yDomain - (@view yMid[j]))) * (@view arcLengths[j]) * (@view TPHI[j])
end

wave = waveAtDomain + potentialContribution; # Huygens-like principle
# return -inv(Mij)
# T = -inv(Mij)
# return T[id]
return wave
end

"""
  `boundaryWallWave(..., potentialFunction::Float64)`

Computes scattered field by a wall at domain points with a scalar potential strength.
"""
function boundaryWallWave(
  waveVector::SVector{2,Float64},
  xMid::Vector,
  yMid::Vector,
  xDomain::Vector, 
  yDomain::Vector,
  σ::ComplexF64,
  rijMid::Matrix,
  arcLengths::Vector,
  numSegments::Int,
  dindex::StepRange,
  potentialStrength::Union{ComplexF64,Float64},
  )

@warn "This method is deprecated in favour of newer ones that allow band integrated matrices."

# calculates boundary wall method
waveNumber     = norm(waveVector)
waveAtBoundary = incidentWave.(Ref(waveVector), xMid, yMid);
waveAtDomain   = incidentWave.(Ref(waveVector), xDomain, yDomain)

# calculate M matrix
Mij          = @. σ * calcGreenFun(waveNumber, rijMid) * arcLengths[1] # arcLengths # multiplies cols by ds[j]
Mij[dindex] .= @. σ * calcGreenFun(waveNumber, arcLengths/2) * arcLengths[1] # arcLengths;

# TPHI         = IterativeSolvers.gmres(A, waveAtBoundary)
TPHI         = potentialStrength *  solve(LinearProblem(LinearAlgebra.I(numSegments) .- potentialStrength * Mij, waveAtBoundary))

# TPHI          = potentialStrength * inv(LinearAlgebra.I(numSegments) .- potentialStrength * Mij)  * waveAtBoundary ;

potentialContribution = zeros(ComplexF64, length(xDomain));
@inbounds for j in 1:numSegments
  # @inbounds potentialContribution += @. σ * calcGreenFun(waveNumber, hypot(xDomain  - (@view xMid[j]), yDomain - (@view yMid[j]))) * (@view arcLengths[j]) * (@view TPHI[j])
  Gzz = [calcGreenFun(waveNumber, hypot(x  .- (@view xMid[j]), y .- (@view yMid[j]))) for (x,y) in zip(xDomain, yDomain)]
  potentialContribution += σ * (Gzz .* arcLengths[j] .* TPHI[j])
end

wave = waveAtDomain + potentialContribution; # Huygens-like principle
return wave
end



"""
  `calcGreenFun(k, r)`

Returns the 2D free-particle Green function (hankel function order 0)

# Arguments:
- `k::Float64`: waveNumber
- `r::Float64`: distance between observation and source, |r - r'|
"""
function calcGreenFun(k::Float64, r::Float64)
  hankelh1(0,k*r)
end 

"""
  `calcGreenFun(k, ri, rj, t)`

Function shortcut for SVectors, might not be needed but still have it
for clarity

"""
calcGreenFun(waveNumber::Float64, ri::SVector{2, Float64}, rj::SVector{2,Float64}) = calcGreenFun(waveNumber, norm(ri-rj))


"""
  `r(r0, r1, t)`

Equation of a line parametrized in t ∈ (0,1)
"""
function r(r0::SVector{2, Float64},r1::SVector{2, Float64}, t::Float64)
  @assert ((t <= 1.0) &&  (t >= 0.0)) "`t` must be within (0,1)"

  return r0 + t * (r1 - r0) 
end


"""
  divideMeshBoundary(x, y, numSegments)

Divides a boundary for a mesh (not really sure what I used it for)
# Arguments:
- `x::Vector{Float64} `: x coord
- `y::Vector{Float64}`: y coord
- `numSegments::Int`: number of segments
"""
function divideMeshBoundary(x::Vector{Float64}, y::Vector{Float64}, numSegments::Int)
  @assert isodd(numSegments) "n must be odd"
  # calculate arc length T (dt=each segment)
  dx = diag(x .- x',1) 
  dy = diag(y .- y',1)
  dt = @. sqrt(dx^2+dy^2)
  dt = [0.0; dt]

  cumLengths = cumsum(dt)
  totalLength = cumLengths[end]
  segmentLength = totalLength / numSegments
  
  #xDivided = zeros(numSegments+1)
  #yDivided = zeros(numSegments+1)
  xDivided = zeros(numSegments)
  yDivided = zeros(numSegments)
  xDivided[1] = x[1]
  yDivided[1] = y[1]
  
  #for segment in 1:numSegments-1
  for segment in 1:numSegments-1
      targetLength = segmentLength * segment
      
      # Find the index of the point corresponding to the target length
      index = searchsortedfirst(cumLengths, targetLength)  # Returns the index of the first value in a greater than or equal to x
      
      # Interpolate to find the (x, y) coordinates at the target length
      fraction = (targetLength - cumLengths[index - 1]) / (cumLengths[index] - cumLengths[index - 1])
      xSegment = x[index - 1] + fraction * (x[index] - x[index - 1])
      ySegment = y[index - 1] + fraction * (y[index] - y[index - 1])
      
      xDivided[segment+1] = xSegment
      yDivided[segment+1] = ySegment
  end
  # return 
  return [xDivided; x[end]], [yDivided; y[end]]
end

"""
  `calcGreenFun(type::Symbol, θ::Float64, k::Float64, r::Float64)`

Calculates Green's dyads for parallel incidence.

# Arguments:
- `type`: A symbol, either :xx, :xy, :yy, or :zz. It denotes the dyadic of Green's 
        tensor to calculate.
- `θ`:    Incident angle.
- `k`:    Wave number.
- `r`:    Distance |r - r'|
"""
function calcGreenFun(type::Symbol, θ::Float64, k::Float64, r::Float64)
  # make this broadcasted 
  
  if type == :xx
    # if k*r < 1.0
    #   println("approximation used")
    #   return im/pi * (
    #          sin(θ)^2  * (2*eGamma - im*pi + 2log(k*r/2) - 0.25*(-2+2eGamma-im*pi+2log(0.5*k*r))*(k*r)^2) +
    #          cos(2θ)/k * (-2(k*r)^(-2) + 0.5*(-1+2eGamma-im*pi+2log(0.5*k*r)) + 0.03125*(5 - 4eGamma+2im*pi-4log(0.5*k*r))*(k*r)^2)
    #          )
    # end
    return sin(θ)^2 * hankelh1(0, k*r) + im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r)
  elseif type==:xy
    # if k*r < 1.0
    #   println("approximation used")
    #   return im/pi * sin(2θ) * 0.5 * (-4 * (k*r)^(-2) - 1  + 0.0625*(-3+4eGamma-2im*pi+4log(0.5*k*r))*(k*r)^2)
    # end
    return sin(2θ)/2 * hankelh1(2, k*r)
  elseif type==:yy
    # if k*r < 1.0
    #   println("approximation used")

    #   return im/pi * (
    #          cos(θ)^2  * (2*eGamma - im*pi + 2log(k*r/2)- (-2+2eGamma-im*pi+2log(k*r/2))/4*(k*r)^2) - 
    #          cos(2θ)/k * (-2(k*r)^(-2) + 0.5*(-1+2eGamma-im*pi+2log(0.5*k*r)) + 0.03125*(5 - 4eGamma+2im*pi-4log(0.5*k*r))*(k*r)^2)
    #          )
    # end
    return (cos(θ)^2 * hankelh1(0, k*r) - im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r))
  elseif type==:zz
    return hankelh1(0, k*r)
  end
end

# function calcGreenFun(type::Symbol, θ::Float64, k::Float64, r::Float64)
#   # make this broadcasted 
#   # this is still faster than the multiple dispatch (using a struct)
#   if type == :xx
#     return sin(θ)^2 * hankelh1(0, k*r) + im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r)
#   elseif type==:xy
#     return sin(2θ)/2 * hankelh1(2, k*r)
#   elseif type==:yy
#     return (cos(θ)^2 * hankelh1(0, k*r) - im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r))
#   elseif type==:zz
#     return hankelh1(0, k*r)
#   end
# end

"""
  `calcGreenTensor(waveNumber::Float64,θ::Float64,rijMid::Matrix,arcLengths::Vector,dindex::StepRange)`

Constructs a block matrix representation for Green's tensor, considering parallel
incidence to the system.

Returns the dyads we need for its construction.

# Arguments: 
- `waveNumber`: wave number
- `θ`: angle of wave incidence (radians)
- `rijMid`: distance matrix
- `arcLengths`: segment length vector (should be homogeneous)
"""
function calcGreenTensor(waveNumber::Float64,
                         θ::Float64,
                         rijMid::Matrix,
                         arcLengths::Vector,
                         dindex::StepRange)
  Gxx          = calcGreenFun.(:xx, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  Gxy          = calcGreenFun.(:xy, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  Gyy          = calcGreenFun.(:yy, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  Gzz          = calcGreenFun.(:zz, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  # Gzz          = @. im / 4 * BoundaryWall.calcGreenFun(waveNumber, rijMid) * arcLengths
  # Gxx          = @. im / 4 * calcGreenFun.(Ref(Xx()), θ, waveNumber, rijMid) 
  
  Gxx[dindex] .= calcGreenFun.(:xx, θ, waveNumber, arcLengths/2) * arcLengths[1]#  arcLengths;
  Gxy[dindex] .= calcGreenFun.(:xy, θ, waveNumber, arcLengths/2) * arcLengths[1]#  arcLengths;
  Gyy[dindex] .= calcGreenFun.(:yy, θ, waveNumber, arcLengths/2) * arcLengths[1]#  arcLengths;
  Gzz[dindex] .= calcGreenFun.(:zz, θ, waveNumber, arcLengths/2) * arcLengths[1]#  arcLengths;

  # return [ones(size(Gxx))-view(Gxx,:,:) view(Gxy,:,:);view(Gxy,:,:) ones(size(Gxx))-view(Gyy,:,:)]
  return Gxx, Gxy, Gyy, Gzz
  # return Gxx, Gxy, Gyy, Gzz
end

function cauchy_quadgk(g, a, b; kws...)
  a < 0.5 < b || throw(ArgumentError("domain must include 0.5"))
  g₀ = g(0.0)
  g₀int = b == -a ? zero(g₀) : g₀ * log(abs(b/a)) / (b - a)
  return quadgk_count(x -> (g(x)-g₀)/x + g₀int, a, 0.5, b; kws...)
end


"""
  `calcGreenTensor(waveNumber::Float64,θ::Float64,rijMid::Matrix,arcLengths::Vector,dindex::StepRange)`

Constructs a block matrix representation for Green's tensor, considering parallel
incidence to the system.

Returns the dyads we need for its construction.

# Arguments: 
- `waveNumber`: wave number
- `θ`: angle of wave incidence (radians)
- `rijMid`: distance matrix
- `arcLengths`: segment length vector (should be homogeneous)
- `band`: band to use for banded matrix
"""
function calcGreenTensor(waveNumber::Float64,
                         θ::Float64,
                         σ::ComplexF64,
                         rijMid::Matrix,
                         arcLengths::Vector,
                         band::Int64,
                         _rm::Vector{SVector{2, Float64}},
                         _r::CircularVector{SVector{2, Float64},Vector{SVector{2, Float64}}},
                         numSegments::Int64,)
  Gxx          = σ * calcGreenFun.(:xx, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  Gxy          = σ * calcGreenFun.(:xy, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  Gyy          = σ * calcGreenFun.(:yy, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  Gzz          = σ * calcGreenFun.(:zz, θ, waveNumber, rijMid) * arcLengths[1] #arcLengths
  # Gzz          = @. im / 4 * BoundaryWall.calcGreenFun(waveNumber, rijMid) * arcLengths
  
  # band elements
  # use _j to index positions
  [Gxx[i,j] = calcMmatrix(:xx, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gxx,1), (j,_j) in zip(axes(Gxx, 2), segmentIterator(_r, numSegments)) if (abs(i - j) < band  &&  i != j)]
  [Gxy[i,j] = calcMmatrix(:xy, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gxy,1), (j,_j) in zip(axes(Gxy, 2), segmentIterator(_r, numSegments)) if (abs(i - j) < band  &&  i != j)]
  [Gyy[i,j] = calcMmatrix(:yy, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gyy,1), (j,_j) in zip(axes(Gyy, 2), segmentIterator(_r, numSegments)) if (abs(i - j) < band  &&  i != j)]

  # [Gxx[i,j] = calcMmatrix(:xx, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gxx,1), (j,_j) in zip(axes(Gxx, 2), segmentIterator(_r, numSegments)) if (abs(i - j) < band)]
  # [Gxy[i,j] = calcMmatrix(:xy, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gxy,1), (j,_j) in zip(axes(Gxy, 2), segmentIterator(_r, numSegments)) if (abs(i - j) < band)]
  # [Gyy[i,j] = calcMmatrix(:yy, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gyy,1), (j,_j) in zip(axes(Gyy, 2), segmentIterator(_r, numSegments)) if (abs(i - j) < band)]
  
  # diagonal elements
  # [Gxx[i,i] = σ * first(quadgk(t -> calcGreenFun(:xx, θ, waveNumber, norm(_rm[i] - r(_r[i], _r[i+1], t))), 0.0,  0.5, 1.0, rtol=1e-3)) * norm(_r[i+1] - _r[i]) for i in axes(Gxx,2)] # scaling
  [Gxx[j,j] = σ * (
                   first(quadgk(t -> calcGreenFun(:xx, θ, waveNumber, norm(_rm[j] - r(_r[_j], _r[_j+1], t))), 0.0,  0.499)) + 
                   first(quadgk(t -> calcGreenFun(:xx, θ, waveNumber,norm(_rm[j] - r(_r[_j], _r[_j+1], t))), 0.501,  1.0))
                  ) * norm(_r[_j+1] - _r[_j]) for (j,_j) in zip(axes(Gxx, 2), segmentIterator(_r, numSegments))]
  
  [Gxy[j,j] = σ * (
                    first(quadgk(t -> calcGreenFun(:xy, θ, waveNumber, norm(_rm[j] - r(_r[_j], _r[_j+1], t))), 0.0,  0.499)) + 
                    first(quadgk(t -> calcGreenFun(:xy, θ, waveNumber, norm(_rm[j] - r(_r[_j], _r[_j+1], t))), 0.501,  1.0))
                   ) * norm(_r[_j+1] - _r[_j]) for (j,_j) in zip(axes(Gxx, 2), segmentIterator(_r, numSegments))]
  
  [Gyy[j,j] = σ * (
                    first(quadgk(t -> calcGreenFun(:yy, θ, waveNumber, norm(_rm[j] - r(_r[_j], _r[_j+1], t))), 0.0,  0.499)) + 
                    first(quadgk(t -> calcGreenFun(:yy, θ, waveNumber, norm(_rm[j] - r(_r[_j], _r[_j+1], t))), 0.501,  1.0))
                   ) * norm(_r[_j+1] - _r[_j]) for (j,_j) in zip(axes(Gxx, 2), segmentIterator(_r, numSegments))]

  # dindex = diagind(size(rijMid)...)
  # [Gxx[i,i] = calcGreenFun(:xx, θ, waveNumber, arcLengths[i]/2) * arcLengths[i] for i in eachindex(arcLengths)]
  # [Gxy[i,i] = calcGreenFun(:xy, θ, waveNumber, arcLengths[i]/2) * arcLengths[i] for i in eachindex(arcLengths)]
  # [Gyy[i,i] = calcGreenFun(:yy, θ, waveNumber, arcLengths[i]/2) * arcLengths[i] for i in eachindex(arcLengths)]
  # [Gyy[i,j] = calcMmatrix(:yy, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gyy,1), (j,_j) in zip(axes(Gyy, 2), segmentIterator(_r, numSegments)) if abs(i - j) < band]
  
  [Gzz[i,j] = calcMmatrix(:zz, θ, waveNumber, σ, _rm[i], _r[_j], _r[_j+1]) for i in axes(Gzz,1), (j,_j) in zip(axes(Gzz, 2), segmentIterator(_r, numSegments)) if abs(i - j) < band]
  
  # return [ones(size(Gxx))-view(Gxx,:,:) view(Gxy,:,:);view(Gxy,:,:) ones(size(Gxx))-view(Gyy,:,:)]
  return Gxx, Gxy, Gyy, Gzz
  # return Gxx, Gxy, Gyy, Gzz
end


"""
  `calcGreenFun(type::Symbol, θ::Float64, k::Float64, r::Float64)`

Calculates Green's dyads for any incidence angle.

# Arguments:
- `type`: A symbol, either :xx, :xy, :yy, or :zz. It denotes the dyadic of Green's 
        tensor to calculate.
- `θ`:    Incident angle.
- `kB`:   3d wave number.
- `kRho`: 
- `r`:    Distance |r - r'|
"""
function calcGreenFun(type::Symbol, θ::Float64, kB::Float64, kRho::Float64, r::Float64)
  # make this broadcasted 
  # this is still faster than the multiple dispatch (using a struct)
  if type == :xx
    return sin(θ)^2 * hankelh1(0, k*r) + im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r)
  elseif type==:xy
    return sin(2θ)/2 * hankelh1(2, k*r)
  elseif type==:yy
    return (cos(θ)^2 * hankelh1(0, k*r) - im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r))
  elseif type==:xy
    return ()
  elseif type==:zz
    return hankelh1(0, k*r)
  end
end


"""
  `boundaryWallWave(waveVector,xMid, yMid, xDomain, yDomain, σ, rijMid, arcLengths, numSegments, dindex)`

Computes the scattered, vectorial field by a wall (aka the Boundary Wall Method),
assuming an infinite potential strength. Assumes incident field travels parallel
to the x-y plane (SVector{2, Float64})

# Arguments:
- `waveVector::SVector{2,Float64}`: self explanatory
- `xMid::Vector`: x coordinate of segment's midpoints
- `yMid::Vector`: y coordinate of segment's midpoints
- `xDomain::Vector,`: x coords for domain of evaluation
- `yDomain::Vector`: y coords for domain of evaluation
- `σ::ComplexF64`: free particle Green's function parameter
- `rijMid::Matrix`: distance matrix between segments
- `arcLengths::Vector`: length of segments
- `numSegments::Int`: number of segments
- `dindex::StepRange`: diagonal index supplied beforehand
- `incidentVector`: vector with functions for x,y, and z fields, 
                    parameters: k::SVector, x::Float64, y::Float64
"""
function boundaryWallVec(
  waveVector::SVector{2,Float64},
  xMid::Vector{Float64},
  yMid::Vector{Float64},
  xBoundary::Vector{Float64},
  yBoundary::Vector{Float64},
  xDomain::Vector{Float64}, 
  yDomain::Vector{Float64},
  σ::ComplexF64,
  rijMid::Matrix{Float64},
  arcLengths::Vector{Float64},
  numSegments::Int64,
  numSegmentsMod::Int64,
  incidentWave::SVector{3, Function}
  )

# calculates boundary wall method
waveNumber      = norm(waveVector)
waveAngle       = atan(waveVector...)

waveAtBoundaryX = [incidentWave[1](waveVector,x,y) for (x,y) in zip(xMid,yMid)];
waveAtBoundaryY = [incidentWave[2](waveVector,x,y) for (x,y) in zip(xMid,yMid)];
waveAtBoundaryZ = [incidentWave[3](waveVector,x,y) for (x,y) in zip(xMid,yMid)];

waveAtDomainX   = [incidentWave[1](waveVector,x,y) for (x,y) in zip(xDomain, yDomain)];
waveAtDomainY   = [incidentWave[2](waveVector,x,y) for (x,y) in zip(xDomain, yDomain)];
waveAtDomainZ   = [incidentWave[3](waveVector,x,y) for (x,y) in zip(xDomain, yDomain)];

rMid = SVector.(xMid, yMid)
rPos = CircularArray(SVector.(xBoundary, yBoundary)[1:end])
# calculate M matrix

Mxx,Mxy,Myy,Mzz = calcGreenTensor(waveNumber, waveAngle, σ, rijMid, arcLengths, 10, rMid, rPos, numSegmentsMod)
println("new method")
dindex = diagind(size(rijMid)...)
# Mxx,Mxy,Myy,Mzz = σ .* calcGreenTensor(waveNumber, waveAngle, rijMid, arcLengths, dindex)

greensTensor = [Mxx Mxy; Mxy Myy] # ./ arcLengths[1] # normalization stuff

TPHI_TE = - greensTensor \ [waveAtBoundaryX; waveAtBoundaryY]

TPHI_X = TPHI_TE[1:numSegments]
TPHI_Y = TPHI_TE[numSegments+1:end]

TPHI_TM = -Mzz \ waveAtBoundaryZ

# potentialContribution = mapreduce( (j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1] , +, 1:numSegments)
# wave =  waveAtDomain + ThreadsX.mapreduce((j) -> σ * calcGreenFun.(waveNumber, hypot.(xDomain.-view(xMid,j)[1], yDomain.-view(yMid,j)[1])) * view(arcLengths,j)[1] * view(TPHI,j)[1],+, 1:numSegments)
waveX = waveAtDomainX
waveY = waveAtDomainY
waveZ = waveAtDomainZ

# RR = zeros()
# Gxy = zeros(length(xDomain)) or add local to Gxy
@inbounds for j in 1:numSegments
  RR = [hypot(x  .- (@view xMid[j]), y .- (@view yMid[j])) for (x,y) in zip(xDomain, yDomain)]

  Gxy = calcGreenFun.(:xy, waveAngle, waveNumber, RR)

  
  waveX += @inbounds @. σ * (calcGreenFun(:xx, waveAngle, waveNumber, RR) * (@view TPHI_X[j]) + Gxy  * (@view TPHI_Y[j])) * (@view arcLengths[j]) 
  waveY += @inbounds @. σ * (Gxy * (@view TPHI_X[j]) + calcGreenFun(:yy, waveAngle, waveNumber, RR)  * (@view TPHI_Y[j])) * (@view arcLengths[j]) 
  waveZ += @inbounds @. σ *  calcGreenFun(:zz, waveAngle, waveNumber, RR) * (@view arcLengths[j]) * (@view TPHI_TM[j])
end

# wave =  + potentialContribution; # Huygens-like principle
# return -inv(Mij)
# T = -inv(Mij)
# return T[id]
return waveX, waveY, waveZ
end


"""
  `calcStokes(Ex::ComplexF64, Ey::ComplexF64)`

Calculate stokes parameters for a TE field
"""
function calcStokes(Ex::ComplexF64, Ey::ComplexF64)
  s0 = abs2(Ex) + abs2(Ey)
  return s0, (abs2(Ex) - abs2(Ey)) / s0, 2*real(conj(Ex)*Ey)/s0, 2*imag(conj(Ex)*Ey)/s0
end

"""
  `calcGreenFun(type::Symbol, θ::Float64, k::Float64, r::Float64)`

Calculates Green's dyads for parallel incidence.

# Arguments:
- `type`: A symbol, either :xx, :xy, :yy, or :zz. It denotes the dyadic of Green's 
        tensor to calculate.
- `θ`:    Incident angle.
- `k`:    Wave number.
- `r`:    Distance |r - r'|
"""
function calcGreenFunSmall(type::Symbol, θ::Float64, k::Float64, r::Float64)
  # make this broadcasted 
  
  if type == :xx
    if r < 0.1
      println("approximation used")
      return im/pi * (
             sin(θ)^2  * (2*eGamma - im*pi + 2log(k*r/2) - 0.25*(-2+2eGamma-im*pi+2log(0.5*k*r))*(k*r)^2) +
             cos(2θ)/k * (-2(k*r)^(-2) + 0.5*(-1+2eGamma-im*pi+2log(0.5*k*r)) + 0.03125*(5 - 4eGamma+2im*pi-4log(0.5*k*r))*(k*r)^2)
             )
    end
    return sin(θ)^2 * hankelh1(0, k*r) + im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r)
  elseif type==:xy
    if r < 0.1
      println("approximation used")
      return im/pi * sin(2θ) * 0.5 * (-4 * (k*r)^(-2) - 1  + 0.0625*(-3+4eGamma-2im*pi+4log(0.5*k*r))*(k*r)^2)
    end
    return sin(2θ)/2 * hankelh1(2, k*r)
  elseif type==:yy
    if r < 0.1
      println("approximation used")

      return im/pi * (
             cos(θ)^2  * (2*eGamma - im*pi + 2log(k*r/2)- (-2+2eGamma-im*pi+2log(k*r/2))/4*(k*r)^2) - 
             cos(2θ)/k * (-2(k*r)^(-2) + 0.5*(-1+2eGamma-im*pi+2log(0.5*k*r)) + 0.03125*(5 - 4eGamma+2im*pi-4log(0.5*k*r))*(k*r)^2)
             )
    end
    return (cos(θ)^2 * hankelh1(0, k*r) - im/4 * cos(2θ)/(k*r)*hankelh1(1, k*r))
  elseif type==:zz
    return hankelh1(0, k*r)
  end
end

function polarizedField(k::SVector{2, Float64}, x::Float64, y::Float64, ϕ::Float64)
  return exp(im*ϕ)*exp(-im * (k[1] * x + k[2] * y))
end

"""
My own gradient function for 2d rectangular meshes.

In principle, the gradient can be calculated using the BWM
"""
function gradient(xdom::LinRange, ydom::LinRange, z::Vector{ComplexF64})
  @assert length(z)==length(xdom)*length(ydom) "Array's length don't match"
  
  ms = [SVector(x,y) for x in xdom, y in ydom]
  itp = linear_interpolation((xdom,ydom),reshape(z,length(xdom), length(ydom)))
  grad = [Interpolations.gradient(itp, idx...) for idx in knots(itp)]
  dv       = imag( conj(reshape(z, length(xdom), length(ydom))) .* grad )
  strength = vec(@. sqrt(first(dv) ^ 2 + last(dv)^ 2))
  
  xs = first.(ms)
  ys = last.(ms)
  
  us = first.(dv)
  vs = last.(dv)
  
  return xs, ys, us, vs
end

end
