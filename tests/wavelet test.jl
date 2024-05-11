# wavelet test
using Wavelets
using GLMakie
using StaticArrays
using LinearAlgebra
using BoundaryWall
using CircularArrays
using LinearSolve
using FFTW

import BoundaryWall as BWM

J = 10
N = 2^J
T = LinRange(-pi, pi, N)
waveVector = SVector(5.0, 1.0)

x,y,xm,ym,ds = createEllipse(1.0,2.0, T,0.0, SVector(0.0, 0.0))
rij = calcDistances(xm,ym)
# waveAtBoundary = [gaussianWave(waveVector, r, 10.0, abstol=1e-5) for r in SVector.(xm, ym)]
waveAtBoundary = [planeWave(waveVector, r,) for r in SVector.(xm, ym)]


Mij = BWM.calcMmatrix(norm(waveVector),
                     -0.25im,
                     rij,
                     ds,
                     6,
                     SVector.(xm, ym),
                     CircularArray(SVector.(x,y)[1:end]),
                     N)



function calcFFTSystem(M::Matrix{ComplexF64}, v::Vector{ComplexF64})
  M_fft = fft(M)
  M_fft[abs2.(M_fft).<0.01] .= 0
  k = solve(LinearProblem(M_fft,fft(v)), KrylovJL_GMRES()).u
  # y_inv = idwt(M_fft \ dwt(v, wt2), wt2)
  ifft!(k)
  return k  
end


@allocated y=Mij\waveAtBoundary
@allocated yf=calcFFTSystem(Mij, waveAtBoundary)

lines(abs2.(y))


L = LinearSolve.LinearProblem(Mij, waveAtBoundary)
S = solve(L, LinearSolve.QRFactorization)
S.u

begin
# ax = Axis(fig[1,1])n

n = size(Mij,1)

M = Mij
b = waveAtBoundary
wt =  wavelet(WT.cdf97, WT.Lifting)
M_wavelet = fft(M)
M_wavelet[abs2.(M_wavelet) .< 0.05] .= 0
b_wavelet = fft(b)

x = M\b
x_wavelet = inv(M_wavelet)*b_wavelet
x_new = ifft(x_wavelet)

let
fig = Figure()

ga = fig[1,1]=GridLayout()
heatmap(ga[1,1], abs2.(M))
heatmap(ga[1,2], abs2.(M_wavelet))  
lines(ga[2,1],abs2.(x))
lines(ga[2,2],abs2.(reverse(x_new)))
# ax.
# lines(ga[2,2], x)
fig
end
end

using BenchmarkTools
