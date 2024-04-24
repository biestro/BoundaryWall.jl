using GLMakie
using SpecialFunctions

function simps(y::Vector, x::Union{Vector,LinRange})
  n = length(y)-1
  n % 2 == 0 || error("`y` length (number of intervals) must be odd")
  length(x)-1 == n || error("`x` and `y` length must be equal")
  h = (x[end]-x[1])/n
  s = sum(slice(y,1:2:n) + 4slice(y,2:2:n) + slice(y,3:2:n+1))
  return h/3 * s
end

G0(r::Vector, rp::Vector) = hankelh1(1,norm(r - rp)) / norm(r - rp)


include("../src/GeometryUtils.jl")

circ = GeometryUtils.createCircle(1.0, LinRange(-pi,pi, 100), (0.0, 0.0))
x  =getindex.(circ, 1)
y  =getindex.(circ, 2)
xm =getindex.(circ, 3)
ym =getindex.(circ, 4)
ds =getindex.(circ, 5)

using QuadGK

curve(t::Float64) = [cos(t), sin(t)]

first(quadgk(t->G0([1.0,0.01], curve(t)), 0,0.001*pi))

ds

boundary_point    = [1.0, 1.0]
observation_point = [1.1, 1.0]

boundary_start = []


G0(x1, x2)

nx= 20
x = LinRange(0.01,5,nx)
y = hankelh1.(0,x)




let
fig = Figure()
ax  = Axis(fig[1,1])
lines!(ax, x,abs2.(y))
fig
end