using LinearAlgebra
using StaticArrays: SVector
using SpecialFunctions
# using GLMakie; GLMakie.activate!(inline=false, float=true)
using CairoMakie
import BoundaryWall as BWM

begin # CONSTANTS
  N    = 9                         # number of points per element
  HBAR = 1.0
  MASS = 0.5
  HBAR = 1.0
  SIGMA= (2*MASS/HBAR^2)*(1/4*im)
  NDOM = 100
  zero = 13.3237
  R    = 1.0
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
  STEP = 2.0R + R/2 # diameter  + constant
  N_CIRCLES = (10,7)
  RANGES = [-(N_CIRCLES[1]-1)*STEP/2:STEP:(N_CIRCLES[1]-1)*STEP/2,-(N_CIRCLES[2]-1)*STEP/2:STEP:(N_CIRCLES[2]-1)*STEP/2]
  N_STEPS =length(RANGES)
  CENTERS = vec([(i, j) for i in RANGES[1], j in RANGES[2]])

  # CENTERS  = vec([(i, j) for i in RANGES, j in RANGES])
  # STRENGTH = ones(length(RANGES), length(RANGES));
  # reverse!(STRENGTH, dims=1)
  #STRENGTH=STRENGTH[:] .* 10.0e10
  # STRENGTH=STRENGTH[:]
  # INDICES  = sort([N_STEPS-5,
  #                 2N_STEPS-5,
  #                 3N_STEPS-5,
  #                 4N_STEPS-5,
  #                 5N_STEPS-5,
  #                 6N_STEPS-5, 6N_STEPS-4,6N_STEPS-3,6N_STEPS-6,6N_STEPS-7,
  #                 7N_STEPS-3, 7N_STEPS-7,
  #                 8N_STEPS-3, 8N_STEPS-7,
  #                 9N_STEPS-3, 9N_STEPS-7,
  #                10N_STEPS-3,10N_STEPS-7,
  #                11N_STEPS-3,11N_STEPS-7])


  # STRENGTH[INDICES] .= rand(length(INDICES))
  INDICES = sort([31,32,33,43,44,45,46,47,48,38,39, 40, 23, 24, 25, 26, 27, 28])
  
  deleteat!(CENTERS, INDICES)

  # POTENTIAL_STRENGTH = repeat(STRENGTH, inner=N)
  CIRCLES = [BWM.createCircle(R, θ, SVector(cen)) for cen in CENTERS]
  x = vcat(getindex.(CIRCLES, 1)...)
  y = vcat(getindex.(CIRCLES, 2)...)
  xm= vcat(getindex.(CIRCLES, 3)...)
  ym= vcat(getindex.(CIRCLES, 4)...)
  ds= vcat(getindex.(CIRCLES, 5)...)
  rij = BWM.calcDistances(xm,ym)
  fig = Figure(theme=BWM.theme)
  ax = Axis(fig[1,1])

  for cen in CENTERS
    circ = BWM.createCircle(R, θ, SVector(cen))
  # [scatter!(ax, BWM.createCircle(R, θ, SVector(cen))) for cen in CENTERS]
  lines!(ax, circ[1], circ[2], linestyle=:solid)
  end
  text!(ax, CENTERS; text=string.(eachindex(CENTERS)))
  arrows!(ax,[-12.0],[0.0], [-cosd(TH)],[-sind(TH)], align=:origin)
  ax.aspect=DataAspect()
  # save("docs/src/assets/photonic_diagram.png", fig)
  save("docs/src/assets/photonic_diagram.svg", fig)
  fig
end



# domain
x0, xf = (-15.,15.)
y0, yf = (-10.,10.)
xdom = LinRange(x0, xf, NDOM)
ydom = LinRange(y0, yf, NDOM)
COORDS = [(_x,_y) for _x in xdom, _y in ydom]
XDOM, YDOM = first.(COORDS)[:], last.(COORDS)[:]

@inline function wave_function(_freq::Float64, _kvec::SVector{2, Float64}, _width::Float64)
  return abs2.(reshape(
    BWM.boundaryWallWave(_freq * _kvec, (k,r)->BWM.gaussianWave(k,r - SVector(-10.25, 0.0), _width; abstol=1e-3), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, BANDED, POTENTIAL_STRENGTH),
    NDOM, NDOM
    )
  )
end
FREQS = [1.35, 1.55, 1.65]

waves = [wave_function(f, KVEC, 4.0) for f in FREQS]
idx_sensor = NDOM - 3
  
begin
  fig = Figure(theme=BWM.theme, size=(600,600))
  gl = fig[1,1] = GridLayout()
  COLOR_RANGE = (0,maximum([maximum(w) for w in waves]))

  ax_l = Makie.Axis[]
  ax_r = Makie.Axis[]
  for (i,w) in enumerate(waves)
    ax1,hm = heatmap(gl[i,1], xdom, ydom, w, colormap=BWM.cmap, colorrange=COLOR_RANGE)
    ax2,_=lines(gl[i,2], w[idx_sensor,:], ydom, color=:blue)
    append!(ax_l, [ax1])
    append!(ax_r, [ax2])
    ax1.xautolimitmargin = (0.05f0, 0.05f0)
    ax1.yautolimitmargin = (0.05f0, 0.05f0)

    ax1.aspect=DataAspect()
    xlims!(ax2, -10, 100.0)
   
    # draw circles
    for cen in CENTERS; 
      circ = BWM.createCircle(R, θ, SVector(cen));  
      lines!(ax1, circ[1], circ[2], linestyle=:solid)
    end

    # sensor
    lines!(ax1, [xdom[idx_sensor], xdom[idx_sensor]], [ydom[1], ydom[end]], color=:blue, linestyle=:dash)
  end
  colgap!(gl, 10)

  [hidexdecorations!(ax, ticks=false, minorticks=false) for ax in ax_l[1:2]]
  [hidexdecorations!(ax, ticks=false, minorticks=false) for ax in ax_r[1:2]]
  [hideydecorations!(ax, ticks=false, minorticks=false) for ax in ax_r]
  # [ax.yaxisposition = :right for ax in ax_r]
  [ax.yticksmirrored = true for ax in ax_r]
  ax_l[1].title="Density"
  ax_r[1].title="Sensor readout"
  ax_l[3].xlabel="x"
  ax_r[3].xlabel="ΨΨ*"
  # colgap!(fig.layout, Relative(0.02))
  # rowsize!(fig.layout,1, Aspect(1, 1.))
  colsize!(gl,1, Aspect(1, 1.8))
  Colorbar(fig[1,2], label="|Ψ|²", labelrotation=0, colorrange=COLOR_RANGE, colormap=BWM.cmap, ticks=(collect(COLOR_RANGE), ["min", "max"]))

  # ax2 = Axis(fig[1,2], tellwidth=true)
  # heatmap!(ax, XDOM, YDOM, abs2.(reshape))
  # scatter!(ax, x,y,color=:black, markersize=5)

  save("crystal.pdf", fig, px_per_unit=3)
  fig
end