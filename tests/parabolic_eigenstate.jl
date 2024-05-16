using Meshes
using GLMakie
using StaticArrays

using SpecialFunctions


using HypergeometricFunctions

parabolicEven(a::Union{ComplexF64, Float64}, z::Union{ComplexF64, Float64}) = exp(-z^2/4) * pFq((0.5 * a + 0.25,), (0.5,), 0.5*z^2)
parabolicOdd(a::Union{ComplexF64, Float64}, z::Union{ComplexF64, Float64}) = exp(-z^2/4) * pFq((0.5 * a + 0.75,), (1.5,), 0.5*z^2)

function cart2par(x::Float64,y::Float64)
  # convert cartesian to parabolic coordinates
  # convention to have the positive x axis as the 
  # branch cut
  z = -x + 1im * y; # convert to imaginary unit
  z = sqrt(2*z);
  _xi = imag(z);
  _eta = real(z);
  return SVector(_xi, _eta)
end 

x = LinRange(-2.0, 2.0, 100)
y = LinRange(-2.0, 2.0, 100)

waveNumber = k = 1.75
a = 5.5
psiEven(_a::Float64, _k::Float64, r::SVector{2,Float64}) = parabolicEven(_a, sqrt(2_k)*r[1]) * parabolicEven(-_a, sqrt(2_k)*r[2])
wave = [psiEven(a,k,cart2par(_x,_y)) for _x in x, _y in y]

let 
  y1 = [parabolicOdd(1.0, _x) for _x in x]
  fig = Figure()
  ax = Axis(fig[1,1])
  # lines!(ax, x, real(y1))
  # contour!(ax, abs2.(wave))
  heatmap!(ax, x, y, abs2.(wave))
  fig
end