"""
This tutorial adds interactivity to the multiplexer
"""


using LinearAlgebra
using StaticArrays: SVector, SizedMatrix
using SpecialFunctions
using GLMakie; GLMakie.activate!(inline=false, float=true)
import BoundaryWall as BWM

begin # CONSTANTS
  N    = 11                         # number of points per element
  HBAR = 1.0
  MASS = 0.5
  HBAR = 1.0
  SIGMA= (2*MASS/HBAR^2)*(1/4*im)
  NDOM = 50
  zero = 13.3237
  R    = 1.0
  θ    = LinRange(0, 2pi, N+1)
  TH   = 180                  # incident angle
  # KVEC = zero/R * SVector(cosd(TH), sind(TH))
  KVEC = SVector(cosd(TH), sind(TH))
  FREQ = 1.0
  # circle 1
end

begin
  # WINDOW = 5.0
  # GAMMA_MAX = 30.0
    # show model
  STEP = 2.25R + R/2 # diameter  + constant
  N_CIRCLES = 9
  RANGES = -(N_CIRCLES-1)*STEP/2:STEP:(N_CIRCLES-1)*STEP/2
  N_STEPS =length(RANGES)
  CENTERS = vec([(i, j) for i in RANGES, j in RANGES])

  CENTERS  = vec([(i, j) for i in RANGES, j in RANGES])
  STRENGTH = ones(length(RANGES), length(RANGES));
  reverse!(STRENGTH, dims=1)
  #STRENGTH=STRENGTH[:] .* 10.0e10
  STRENGTH=STRENGTH[:]

  INDICES = sort([64,65,66,67,68,69,60,51,40,41,42,43,44,45,31,22,13,12,11,10])
  
  deleteat!(CENTERS, INDICES)

  POTENTIAL_STRENGTH = repeat(STRENGTH, inner=N)
  CIRCLES = [BWM.createCircle(R, θ, SVector(cen)) for cen in CENTERS]
  x = vcat(getindex.(CIRCLES, 1)...)
  y = vcat(getindex.(CIRCLES, 2)...)
  xm= vcat(getindex.(CIRCLES, 3)...)
  ym= vcat(getindex.(CIRCLES, 4)...)
  ds= vcat(getindex.(CIRCLES, 5)...)
  rij = BWM.calcDistances(xm,ym)
  fig = Figure(theme=BWM.theme)
  ax = Axis(fig[1,1])
  # scatter!(ax,xm, ym, markersize=5)
  for cen in CENTERS
    circ = BWM.createCircle(R, θ, SVector(cen))
  # [scatter!(ax, BWM.createCircle(R, θ, SVector(cen))) for cen in CENTERS]
  lines!(ax, circ[1], circ[2], linestyle=:solid)
  end
  # text!(ax, CENTERS[INDICES]; text=string.(INDICES))

  text!(ax, CENTERS; text=string.(eachindex(CENTERS)))

  # text!(ax, CENTERS; text=string.(STRENGTH))
  inset_ax = Axis(fig[1, 1],
                  width=Relative(0.1),
                  height=Relative(0.1),
                  halign=0.0,
                  valign=1.0,
                  backgroundcolor=:snow2,
                  # xautolimitmargin = (0.2, 0.2),
                  # yautolimitmargin = (0.2, 0.2),
  )
  # inset_ax.aspect = DataAspect()
  inset_ax.limits=(-1,1,-1,1)
  hidedecorations!(inset_ax)
  # hidedecorations!(inset_ax)
  translate!(inset_ax.scene, 0, 0, 10)
  # translate!(inset_ax.elements[:background], 0, 0, 9)
  arrows!(inset_ax,[cosd(TH)/2],[sind(TH)/2], [-cosd(TH)],[-sind(TH)], align=:origin)
  

  ax.aspect=DataAspect()
  inset_ax.aspect=(DataAspect())
  fig
end



# domain
x0, xf = (-2+floor(minimum(x)),2+ceil(maximum(x)))
y0, yf = (-2+floor(minimum(y)),2+ceil(maximum(y)))
xdom = LinRange(x0, xf, NDOM)
ydom = LinRange(y0, yf, NDOM)
COORDS = [(x,y) for x  in xdom, y in ydom]
XDOM, YDOM = first.(COORDS)[:], last.(COORDS)[:]

begin

  @inline function double_beam(k,r, w, a1::Int,a2::Int)
    # coeffs
    return a1*BWM.gaussianWave(k,r - SVector(-14.0, 10.5), w; abstol=1e-3) + a2*BWM.gaussianWave(k,r - SVector(-14.0, -10.5), w; abstol=1e-3)
  end

  @inline function wave_function(_freq::Float64, _kvec::SVector{2, Float64}, _width::Float64, _coeff::SVector{2, Int64})
    return abs2.(
        reshape(
      BWM.boundaryWallWave(_freq * _kvec, (k,r) -> double_beam(k,r, _width, _coeff...), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, POTENTIAL_STRENGTH),
      NDOM, NDOM
        )
      )
  end

  
  POTENTIAL_STRENGTH = Inf
  DELTA_FREQ         = 0.05
  DELTA_WIDTH        = 0.1
  DELTA_ANGLE        = 15 # degrees
  TH   = 180               
  MAX_FREQ           = 2.0

  frequency_obs = Observable(1.0)
  width_obs     = Observable(3.0)
  toggle_obs    = Observable(SVector(0,1))
  frequency_str = Observable("Freq: $(round(frequency_obs.val, digits=3))")
  theta_obs = Observable(TH)
  angle_obs = Observable(SVector(cosd(theta_obs[]), sind(theta_obs[])))
  angle_str = Observable("Angle: $(round(theta_obs.val, digits=3))")

  # arrow     = Observable([[Point2f(1,1)], [Point2f(2,2)]])
  # point1y   = Observable(sind(theta_obs[])/2)
  # KVEC = zero/R * SVector(cosd(TH), sind(TH))
  # KVEC = SVector(cosd(TH), sind(TH))
  banded = 3
  fig = Figure(theme=BWM.theme)

  wave = Observable(
          (
            # BWM.boundaryWallWave(frequency_obs[] * KVEC, (k,r)->BWM.gaussianWave(k,r - SVector(-10.0, 0.0), width_obs[]; abstol=1e-3), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, -3.0);
            wave_function(frequency_obs[], angle_obs[], width_obs[], toggle_obs[])
          )
        )
  wave_end = Observable(wave[][end,:])


  ax = Axis(fig[1,1])
  ax.aspect=DataAspect()
  # viz!(ax, MESH, color=wave, colormap=:linear_kbgyw_5_98_c62_n256) #:linear_protanopic_deuteranopic_kbw_5_98_c40_n256
  # viz!(ax, MESH, color=wave, colormap=:linear_protanopic_deuteranopic_kbw_5_98_c40_n256)

  heatmap!(ax, xdom, ydom, wave, colormap=BWM.cmap)
  # ax2 = Axis(fig[1,2], tellwidth=true)
  # heatmap!(ax, XDOM, YDOM, abs2.(reshape))
  scatter!(ax, x,y,color=:black, markersize=5)


  fig[2, 1] = buttongrid = GridLayout(tellwidth = false)
  buttons = buttongrid[1, 1:4] = [Button(fig,label="+ f"), Button(fig, label="- f"),Button(fig,label="Toggle A"), Button(fig, label="Toggle B") ]


  gl_labels = fig[0,1] = GridLayout()
  # gl_output = fig[1,2] = GridLayout(tellheight=false)
  ax2 = Axis(fig[1,2], autolimitaspect=1)
  lines!(ax2,  wave_end, ydom)
  ax2.yticklabelsvisible=false
  # ax2.aspect=DataAspect()
  # ax2.aspect=(2,1)
  colsize!(fig.layout, 1, Aspect(2, 10.0))
  # ax_arrow,_ = arrows(gl_labels[1,1], arrow, align=:origin)
  # ax_arrow.aspect=DataAspect()
  # ax_arrow.limits=(-1,1,-1,1)

  # hidedecorations!(ax_arrow)


  Label(gl_labels[1,1], frequency_str, font = :italic, tellheight=true, tellwidth=false)
  Label(gl_labels[1,2], angle_str, font = :italic, tellheight=true, tellwidth=false)

  on(buttons[1].clicks) do _
    if abs(frequency_obs[]) <= MAX_FREQ
    frequency_obs[] += DELTA_FREQ
    end
    frequency_str[] = "Freq: $(round(frequency_obs.val, digits=3))"
    wave[] = (
      wave_function(frequency_obs[], angle_obs[], width_obs[], toggle_obs[])
    )
    wave_end[] = wave[][end,:]
    # end
    # return "frequency_obs[]"
  end

  # Label(gl_labels[2,1], "titleposition = :top\norientation = :vertical\nnbanks = 2", font = :italic, tellheight=false)

  on(buttons[2].clicks) do _
    if frequency_obs[] >= 0.0 # MIN FREQ
    frequency_obs[] -=  DELTA_FREQ
    end
    frequency_str[] = "Freq: $(round(frequency_obs.val, digits=3))"

    wave[] = (
            wave_function(frequency_obs[], angle_obs[], width_obs[], toggle_obs[])
    )
    wave_end[] = wave[][end,:]

    # ax2.autolimitaspect()

  end

  # TOGGLES A
  on(buttons[3].clicks) do _

    toggle_obs[] = SVector(1 - toggle_obs[][1], toggle_obs[][2])

    wave[] = (
            wave_function(frequency_obs[], angle_obs[], width_obs[], toggle_obs[])
    )
    wave_end[] = wave[][end,:]

  end

  # TOGGLES B
  on(buttons[4].clicks) do _

    toggle_obs[] = SVector(toggle_obs[][1], 1 - toggle_obs[][2])
    
    wave[] = (
            wave_function(frequency_obs[], angle_obs[], width_obs[], toggle_obs[])
    )
    # println(width_obs[])
    wave_end[] = wave[][end,:]

  end


  # sliderobservables = [s.value for s in sg.sliders]
  # _ = lift(sliderobservables...) do slvalues...
  #   _freq, _theta = slvalues
  #   wave[] = real(
  #           BWM.boundaryWallWave(_freq * KVEC, (k,r)->BWM.gaussianWave(k,r - SVector(-10.0, 0.0), 3.0; abstol=1e-3), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, banded, -3.0);
  #   )
  # end

  fig
end
