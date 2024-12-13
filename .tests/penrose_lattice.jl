using LinearAlgebra
using StaticArrays: SVector, MVector
using GLMakie; GLMakie.activate!(inline=false, float=true)
using GLMakie
import BoundaryWall as BWM

compose(f::Function, n::Int) = reduce( ∘, fill(f,n))

function penrose_centers(num_compositions::Int64, _rad::Float64)
  # from Rosetta Code
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

  x, y, theta, r, svglines, stack = 0, 0, π / 5, _rad, String[], Vector{Real}[]

  x_vec = Float64[]
  y_vec = Float64[]
  for c in split(penrose, "")
      if c == "A"
          xx, yy = x + r * cos(theta), y + r * sin(theta)
          # line = @sprintf("<line x1='%.1f' y1='%.1f' x2='%.1f' y2='%.1f' style='stroke:rgb(255,165,0)'/>\n", x, y, xx, yy)
          x, y = xx, yy
          # push!(svglines, line)
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
  return SVector.(x_vec, y_vec)
end

begin # CONSTANTS
  N    = 20
  HBAR = 1.0
  MASS = 0.5
  HBAR = 1.0
  SIGMA= (2*MASS/HBAR^2)*(1/4*im)
  NDOM = 50
  zero = 13.3237
  R    = 0.1
  θ    = LinRange(0, 2pi, N+1)
  TH   = 180                  # incident angle
  # KVEC = zero/R * SVector(cosd(TH), sind(TH))
  KVEC = SVector(cosd(TH), sind(TH))
  POTENTIAL_STRENGTH = -100.0
  BANDED = 3


  # circle 1
end

begin
  # WINDOW = 5.0
  # GAMMA_MAX = 30.0
    # show model
    n_penrose = 3
  STEP = 10.0R
  N_CIRCLES = (10,7)
  RANGES = [-(N_CIRCLES[1]-1)*STEP/2:STEP:(N_CIRCLES[1]-1)*STEP/2,-(N_CIRCLES[2]-1)*STEP/2:STEP:(N_CIRCLES[2]-1)*STEP/2]
  N_STEPS =length(RANGES)
  CENTERS = penrose_centers(n_penrose, STEP)
  
  CENTERS = @. MVector(round(first(CENTERS), digits=5),round(last(CENTERS), digits=5))
  sort!(CENTERS, by = norm)
  replace!.(CENTERS,-0.0=>0.0)
  
  unique!(CENTERS)
  center_wave = popat!(CENTERS, 1)
 


  CIRCLES = [BWM.createCircle(R, θ, SVector(cen)) for cen in CENTERS]
  x = vcat(getindex.(CIRCLES, 1)...)
  y = vcat(getindex.(CIRCLES, 2)...)
  xm= vcat(getindex.(CIRCLES, 3)...)
  ym= vcat(getindex.(CIRCLES, 4)...)
  ds= vcat(getindex.(CIRCLES, 5)...)
  rij = BWM.calcDistances(xm,ym)
  fig = Figure(theme=merge(BWM.theme, theme_dark()))
  ax = Axis(fig[1,1])


 
  for cen in CENTERS
    circ = BWM.createCircle(R, θ, SVector(cen))
  # [scatter!(ax, BWM.createCircle(R, θ, SVector(cen))) for cen in CENTERS]
  lines!(ax, circ[1], circ[2], linestyle=:solid)
  end
  text!(ax, CENTERS; text=string.(eachindex(CENTERS)))
  # arrows!(ax,[-12.0],[0.0], [-cosd(TH)],[-sind(TH)], align=:origin)
  ax.aspect=DataAspect()
  # save("docs/src/assets/photonic_diagram.png", fig)
  # save("docs/src/assets/photonic_diagram.svg", fig)

  fig
end




# domain 
x0, xf = (-5.0,5.0)
y0, yf = (-5.0,5.0)
xdom = LinRange(x0, xf, NDOM)
ydom = LinRange(y0, yf, NDOM)
COORDS = [(_x,_y) for _x in xdom, _y in ydom]
XDOM, YDOM = first.(COORDS)[:], last.(COORDS)[:]

# simulation
r_fun(r::SVector{2, Float64}) = SVector(r[1] - center_wave[1], r[2] - center_wave[2])
wave = BWM.boundaryWallWave(3.0 * KVEC, (k,r)->BWM.besselWave(norm(k), hypot(r_fun(r)...),2), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, BANDED, Inf)

# plotting
begin
fig = Figure(theme=theme_dark())
ax = Axis(fig[1,1])
heatmap!(ax, xdom, ydom, abs2.(reshape(wave, length(xdom), length(ydom))), colormap=BWM.cmap_inv)
for cen in CENTERS
  circ = BWM.createCircle(R, θ, SVector(cen))
# [scatter!(ax, BWM.createCircle(R, θ, SVector(cen))) for cen in CENTERS]
lines!(ax, circ[1], circ[2], linestyle=:solid, color=:gray50)
end
ax.aspect=DataAspect()
fig
end
