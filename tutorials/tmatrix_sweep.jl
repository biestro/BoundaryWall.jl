using CairoMakie
using LinearAlgebra
using Meshes
using CircularArrays

import BoundaryWall as BWM
using StaticArrays

begin # definitions
  HBAR        = 1.0
  MASS        = HBAR/2
  SIGMA       = (2*MASS/HBAR^2)*(1/4*im)
  N           = 500
  NDOM        = 150
  ϕ           = 135
  waveNumber  = sqrt(2.064)
  waveVector  = waveNumber*SVector(cosd(ϕ), sind(ϕ)) # parabolic billiard eigenstate
end

y, x,ym, xm, distance_matrix, arc_lengths = BWM.createConfocalBilliard(2.0, 3.0, N)

mats = []

waveRange = LinRange(sqrt(2.064), sqrt(2.405), 100);

for wavenum in waveRange

  M=BWM.calcMmatrix(
    wavenum,
    SIGMA,
    distance_matrix,
    arc_lengths,
    2,
    SVector.(xm, ym),
    CircularArray(SVector.(x, y)[1:end]),
    N
  )
  push!(mats, inv(M))
end


heatmap(log.(abs.(mats[end])))