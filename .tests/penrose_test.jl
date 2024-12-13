using Random
using Printf

compose(f::Function, n::Int) = reduce( ∘, fill(f,n))


function drawpenrose(num_compositions::Int64)
  lindenmayer_rules = Dict(
      "A" => "",
      "M" => "OA++PA----NA[-OA----MA]++", 
      "N" => "+OA--PA[---MA--NA]+",
      "O" => "-MA++NA[+++OA++PA]-", 
      "P" => "--OA++++MA[+PA++++NA]--NA")

  rul(x) = lindenmayer_rules[x]

  # penrose = replace(replace(replace(replace("[N]++[N]++[N]++[N]++[N]",
      # r"[AMNOP]" => rul), r"[AMNOP]" => rul), r"[AMNOP]" => rul), r"[AMNOP]" => rul)
  axiom = "[N]++[N]++[N]++[N]++[N]"
  axiom_len = 5
  penrose = compose(k->replace(k, r"[AMNOP]" => rul), num_compositions)(axiom)

  x, y, theta, r, svglines, stack = 0, 0, π / 5, 20.0, String[], Vector{Real}[]

  x_vec = Float64[]
  y_vec = Float64[]
  for c in split(penrose, "")
      if c == "A"
          xx, yy = x + r * cos(theta), y + r * sin(theta)
          line = @sprintf("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' style='stroke:rgb(255,165,0)'/>\n", x, y, xx, yy)
          x, y = xx, yy
          push!(svglines, line)
          push!(x_vec, x); push!(y_vec, y)
      elseif c == "+"
          theta += π / axiom_len
      elseif c == "-"
          theta -= π / axiom_len
      elseif c == "["
          push!(stack, [x, y, theta])
      elseif c == "]"
          x, y, theta = pop!(stack)
          push!(x_vec, x); push!(y_vec, y)

      end
  end

  # svg = join(unique(svglines), "\n")
  # fp = open("penrose_tiling.svg", "w")
  # write(fp, """<svg xmlns="http://www.w3.org/2000/svg" height="350" width="350"> <rect height="100%" """ *
  #           """width="100%" style="fill:black" />""" * "\n$svg</svg>")
  # close(fp)
  return x_vec, y_vec
end

n = 4
X,Y = drawpenrose(5)
distance = hypot.(X,Y)

using GLMakie; GLMakie.activate!(inline=false, float=true)

let 
  fig = Figure(theme=theme_dark())
  ax=Axis(fig[1,1])
  scatterlines!(ax,X,Y, color=sin.(distance))
  ax.aspect=DataAspect()
  fig
end