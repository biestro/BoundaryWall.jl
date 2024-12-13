using DynamicalBilliards
using GLMakie
const SV = SVector{2}

bd = Obstacle{Float64}[]

@inline function Makie.plot!(_ax::Axis, _d::Disk)
  poly!(_ax, Circle(Point2f(_d.c...), _d.r), color=(EDGECOLOR, 0.5), strokewidth=2, strokecolor=EDGECOLOR)
end

@inline function Makie.plot!(_ax::Axis, _e::Ellipse)
  _th = range(-pi,pi,200)
  _x = _e.a*cos.(_th)
  _y = _e.b*sin.(_th)
  #poly!(_ax, Circle(Point2f(_d.c...), _d.r), color=(EDGECOLOR, 0.5), strokewidth=2, strokecolor=EDGECOLOR)
  lines!(_ax, _x,_y, color=EDGECOLOR)
end

@inline function Makie.plot!(_ax::Axis, _w::Wall)
  println("in fun")
  lines!(_ax, first.([_w.sp, _w.ep]), last.([_w.sp, _w.ep]), color=EDGECOLOR)
end

@inline function Makie.plot!(_ax::Axis, _bd::Billiard)
  #[lines!(_ax,first.([_w.sp, _w.ep]), last.([_w.sp, _w.ep]), color=EDGECOLOR) for _w in _bd]
  [plot!(_ax, _o) for _o ∈ _bd]
end


const EDGECOLOR=:black

hexagon_vertex = (r) -> [ [r*cos(2π*i/6), r*sin(2π*i/6)] for i in 1:6]
hexver = hexagon_vertex(2.0)

for i in eachindex(hexver)
  starting = hexver[i]
  ending = hexver[mod1(i+1, length(hexver))]
  w = ending - starting
  normal = [-w[2], w[1]]
  wall = InfiniteWall(starting, ending, normal, "wall $i")
  push!(bd, wall)
end

# billiard = Billiard(bd)
billiard = Billiard(Ellipse(SVector(0.0, 0.0), 1.5, 1.0,false))
x0 = 0.0; y0=0.0;

φ0 = rand()*2pi
p = Particle(x0, y0, φ0)

#ct, poss, vels = evolve(p, billiard, 100)
xt, yt, vxt, vyt, t = DynamicalBilliards.timeseries(p, billiard, 100)


let 
  fig = Figure()  
  ax = Axis(fig[1,1])
  # lines!(ax, xt, yt, linewidth=5)
  lines!(ax, xt, yt, color=1,colormap=(:tab10,0.5), colorrange=(1,10))
  #[lines!(ax,first.([_w.sp, _w.ep]), last.([_w.sp, _w.ep]), color=:black) for _w in box]
  # plot!(ax, bd)
  #[plot!(ax, _w) for _w in bd]
  plot!(ax, billiard)
  
  ax.aspect=DataAspect()
  fig
end


# parabolic billiard

