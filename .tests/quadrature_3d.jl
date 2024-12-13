using HCubature

import BoundaryWall as BWM

_rm = SVector.(xm, ym)
_r = CircularArray(SVector.(x,y)[1:end])
#calcGreenFun(_type, θ, waveVector, norm(_rm[i] - r(_r[_j], _r[_j+1], t)))
_j = 1

_band = 1

MIJ = zeros(ComplexF64,N,N)

using QuadGK
quadgk(t->BWM.calcGreenFun(:yz, atand(waveVector[2], waveVector[1]), waveVector, norm(_rm[i] - BWM.r(_r[_j], _r[_j+1], t))), 0,0.5, 1)
hquadrature(t->BWM.calcGreenFun(:yz, atand(waveVector[2], waveVector[1]), waveVector, norm(_rm[i] - BWM.r(_r[_j], _r[_j+1], t))), 0,1)

to_graph(t) = BWM.calcGreenFun(:yz, atand(waveVector[2], waveVector[1]), waveVector, norm(_rm[i] - BWM.r(_r[_j], _r[_j+1], t)))

T = LinRange(0.0, 1.0, 1000)
lines(T, abs2.(to_graph.(T)))

for i in 1:N
  println("outer for $i")
  for (j,_j) in zip(1:N, BWM.segmentIterator(_r, N))
    if abs(i-j) < _band
    println("outer for $_j")
    MIJ[i,j] = first(hquadrature(t -> BWM.calcGreenFun(:yz, atand(waveVector[2], waveVector[1]), waveVector, norm(_rm[i] - BWM.r(_r[_j], _r[_j+1], t))), 0, 1, maxevals=10000))
    end
  end
end

