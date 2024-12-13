using StaticArrays
const SV = SVector{2,Float64}

abstract type AbstractParticle end

const EDGECOLOR = :black
mutable struct Particle <: AbstractParticle
  pos::SV
  vel::SV
end

# constructor
Particle(x0, y0, φ0) = Particle(SV(x0, y0), SV(cos(φ0), sin(φ0)))

abstract type Obstacle end

struct Wall <: Obstacle
    sp::SV
    ep::SV
    normal::SV
end

struct Disk <: Obstacle
    c::SV
    r::Float64
end

const Billiard = NTuple{N, Obstacle} where N

using LinearAlgebra: dot, normalize

"""
    collision(p::AbstractParticle, o::Obstacle) → t, cp
Find the collision (if any) between given particle and obstacle.
Return the time until collision and the estimated collision point `cp`.
"""
@inline function collision(p::Particle, w::Wall)
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision()
    else
        t = dot(w.sp - p.pos, n)/denom
        return t, p.pos + t * p.vel
    end
end

normalvec(w::Wall, pos) = w.normal



@inline function collision(p::Particle, d::Disk)
  dotp = dot(p.vel, normalvec(d, p.pos))
  dotp ≥ 0.0 && return nocollision()

  dc = p.pos - d.c
  B = dot(p.vel, dc)           #pointing towards circle center: B < 0
  C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
  Δ = B*B - C

  Δ ≤ 0.0 && return nocollision()
  sqrtD = sqrt(Δ)
  # Closest point:
  t = -B - sqrtD
  return t, p.pos + t * p.vel
end

normalvec(d::Disk, pos) = normalize(pos - d.c)

@inline nocollision() = (Inf, SV(0.0, 0.0))



function next_collision(p::AbstractParticle, bd)
  j, ct, cp = 0, Inf, SV(0.0, 0.0)
  for i in eachindex(bd)
      t, c = collision(p, bd[i])
      if t < ct
          j = i
          ct = t
          cp = c
      end
  end
  return j, ct, cp
end

propagate!(p::Particle, pos, t) = (p.pos = pos)


function resolvecollision!(p::AbstractParticle, o::Obstacle)
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
end

"""
    bounce!(p, bd)
Evolve the particle for one collision (in-place).
"""
@inline function bounce!(p::AbstractParticle, bd)
    i::Int, tmin::Float64, cp::SV = next_collision(p, bd)
    if tmin != Inf
        propagate!(p, cp, tmin)
        resolvecollision!(p, bd[i])
    end
    return i, tmin, p.pos, p.vel
end



"""
    timeseries!(p::AbstractParticle, bd, n) -> xt, yt, t
Evolve the particle in the billiard `bd` for `n` collisions
and return the position timeseries `xt, yt` along with time vector `t`.
"""
function timeseries!(p::AbstractParticle, bd, n::Int)

    t = [0.0]; xt = [p.pos[1]]; yt = [p.pos[2]]; c = 0

    while c < n
        prevpos = p.pos; prevvel = p.vel
        i, ct = bounce!(p, bd)
        xs, ys = extrapolate(p, prevpos, prevvel, ct)
        push!(t, ct)
        append!(xt, xs)
        append!(yt, ys)
        c += 1
    end

    return xt, yt, t
end

extrapolate(p::Particle, prevpos, prevvel, ct) = p.pos

x, y, r = 1.0, 1.0, 0.3
sp = [0.0,y]; ep = [0.0, 0.0]; n = [x,0.0]
leftw = Wall(sp, ep, n)
sp = [x,0.0]; ep = [x, y]; n = [-x,0.0]
rightw = Wall(sp, ep, n)
sp = [x,y]; ep = [0.0, y]; n = [0.0,-y]
topw = Wall(sp, ep, n)
sp = [0.0,0.0]; ep = [x, 0.0]; n = [0.0,y]
botw = Wall(sp, ep, n)
disk = Disk([x/2, y/2], r)

box = (botw, rightw, topw, leftw)
bd = (botw, rightw, topw, leftw, disk)
using DynamicalBilliards


el = DynamicalBilliards.Ellipse(SV(0.0, 0.0), 1.0, 1.0)
using DynamicalBilliards
bd = Obstacle[]



p = Particle(0.1, 0.1, 2π*rand())

xt, yt, t = timeseries!(p, bd, 10)

using GLMakie


@inline function Makie.plot!(_ax::Axis, _d::Disk)
  poly!(_ax, Circle(Point2f(_d.c...), _d.r), color=(EDGECOLOR, 0.5), strokewidth=2, strokecolor=EDGECOLOR)
end

@inline function Makie.plot!(_ax::Axis, _w::Wall)
  lines!(_ax, first.([_w.sp, _w.ep]), last.([_w.sp, _w.ep]), color=EDGECOLOR)
end

@inline function Makie.plot!(_ax::Axis, _bd::Billiard)
  #[lines!(_ax,first.([_w.sp, _w.ep]), last.([_w.sp, _w.ep]), color=EDGECOLOR) for _w in _bd]
  for _o ∈ _bd; plot!(_ax, _o); end
end

let 
  fig = Figure()  
  ax = Axis(fig[1,1])
  lines!(ax, xt, yt, linewidth=5)
  #[lines!(ax,first.([_w.sp, _w.ep]), last.([_w.sp, _w.ep]), color=:black) for _w in box]
  plot!(ax, bd)
  ax.aspect=DataAspect()
  fig
end


# own billiard
struct Semicircle{T<:AbstractFloat} <: Circular{T} # <: Obstacle{T}
  c::SVector{2,T} # this MUST be a static vector
  r::T
  facedir::SVector{2,T} # this MUST be a static vector
  name::String # this is an OPTIONAL field
end
