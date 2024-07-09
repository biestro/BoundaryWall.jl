
# obtain symbolic coefficients

# normal coeffs
using OffsetArrays # for simplicity

function myfact(x::Int)
  if x < 21
    return factorial(x)
  else
    return factorial(big(x))
  end
end

Base.:/(_a::OffsetVector, _b::Vector) = _a.parent ./ _b

@inline function evenCoeff(n::Union{Int64, BigInt}, a::Float64)
  _coeffs = [1.0,  a]
  if n < length(_coeffs)
    return _coeffs[1:n+1]
  end
  
  _coeff_array = OffsetArray(zeros(Float64,n+1), 0:n)

  _coeff_array[0:1] .= _coeffs
  @inbounds for _n in eachindex(_coeff_array)[2:n]
    _coeff_array[_n] = a * _coeff_array[_n-1] * (_n-1)*(_n-1.5)*_coeff_array[_n-2]
  end

  return _coeff_array
end

@inline function evenCoeff(n::Union{Int64, BigInt}, a::BigFloat)
  _coeffs = [1.0,  a]
  if n < length(_coeffs)
    return _coeffs[1:n+1]
  end
  
  _coeff_array = OffsetArray(zeros(BigFloat,n+1), 0:n)

  _coeff_array[0:1] .= _coeffs
  @inbounds for _n in eachindex(_coeff_array)[2:n]
    _coeff_array[_n] = a * _coeff_array[_n-1] - (_n-1)*(_n-1.5)*_coeff_array[_n-2]
  end

  return _coeff_array
end

@inline function oddCoeff(n::Union{Int64, BigInt}, a::Float64)
  _coeffs = [1.0,  a]
  if n < length(_coeffs)
    return _coeffs[1:n+1]
  end
  
  _coeff_array = OffsetArray(zeros(Float64,n+1), 0:n)

  _coeff_array[0:1] .= _coeffs
  @inbounds for _n in eachindex(_coeff_array)[2:n]
    _coeff_array[_n] = a * _coeff_array[_n-1] - (_n-1)*(_n-0.5)*_coeff_array[_n-2]
  end

  return _coeff_array
end

@inline function oddCoeff(n::Union{Int64, BigInt}, a::BigFloat)
  _coeffs = [1.0,  a]
  if n < length(_coeffs)
    return _coeffs[1:n+1]
  end
  
  _coeff_array = OffsetArray(zeros(BigFloat,n+1), 0:n)

  _coeff_array[0:1] .= _coeffs
  @inbounds for _n in eachindex(_coeff_array)[2:n]
    _coeff_array[_n] = a * _coeff_array[_n-1] - (_n-1)*(_n-0.5)*_coeff_array[_n-2]
  end

  return _coeff_array
end

function evenSeries(n::Union{Int64, BigInt}, a::Union{BigFloat, Float64})
  _even_coeff  = evenCoeff(n, a)
  _fact_coeff = [myfact(2*_n) for _n in 0:n]

  return OffsetArray(_even_coeff / _fact_coeff, 0:n)
end

function oddSeries(n::Union{Int64, BigInt}, a::Union{BigFloat, Float64})
  _odd_coeff   = oddCoeff(n, a)
  _fact_coeff = [myfact(2*_n+1) for _n in 0:n]

  return OffsetArray(_odd_coeff / _fact_coeff, 0:n)
end

x = collect(LinRange(0,10,200))
a = 5.0
N = 200

even_coeffs_pos = evenSeries(N, BigFloat(a))
even_coeffs_neg = evenSeries(N, BigFloat(-a))

function pcfEven(_coeffs::OffsetVector, x::Float64)

  return sum([_coeffs[n_] * x^(2*n_) for n_ in eachindex(_coeffs)])
end

using BoundaryWall: par2cart, cart2par

even_coeffs_pos[40:end] .= 0.0

y = sum([even_coeffs_pos[n_] * x.^(2*n_) for n_ in eachindex(even_coeffs_pos[1:21])])

lines(x, y)

