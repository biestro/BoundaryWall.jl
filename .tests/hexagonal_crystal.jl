using LinearAlgebra
using StaticArrays: SVector
using SpecialFunctions
# using GLMakie; GLMakie.activate!(inline=false, float=true)
using CairoMakie
import BoundaryWall as BWM
using ThreadsX

# TODO: create optical path length for destructive interference

begin # CONSTANTS
  N    = 13                         # number of points per element
  HBAR = 1.0
  MASS = 0.5
  HBAR = 1.0
  SIGMA= (2*MASS/HBAR^2)*(1/4*im)
  NDOM = 30
  zero = 13.3237
  R    = 0.5
  θ    = LinRange(0, 2pi, N+1)
  TH   = 180                  # incident angle
  # KVEC = zero/R * SVector(cosd(TH), sind(TH))
  KVEC = SVector(cosd(TH), sind(TH))
  POTENTIAL_STRENGTH = Inf
  BANDED = 3


  # circle 1
end

begin
  # WINDOW = 5.0
  # GAMMA_MAX = 30.0
    # show model
  STEP = 5R + R/2 # diameter  + constant
  N_CIRCLES = (15,15)
  centers = BWM.buildGrid(BWM.HexagonalGrid(SVector(0.0, 0.0),STEP), N_CIRCLES...)
  x0 = maximum(first.(centers))/2
  y0 = maximum(last.(centers))/2
  CENTERS = @. SVector(first(centers) - x0, last(centers) - y0)


  # STRENGTH[INDICES] .= rand(length(INDICES))
  INDICES = sort([
    collect(136:150)...,
    collect(151:165)...,
    collect(31:37)...,
    collect(16:21)...,
    51,52,67,68,82,83,98,99,113,114,129,130,
    ])
  
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
  # arrows!(ax,[-12.0],[0.0], [-cosd(TH)],[-sind(TH)], align=:origin)
  ax.aspect=DataAspect()
  # save("docs/src/assets/photonic_diagram.png", fig)
  fig
end



# domain
x0, xf = (-25.,25.)
y0, yf = (-20.,20.)
xdom = [xf]#LinRange(x0, xf, NDOM)
ydom = LinRange(y0, yf, NDOM)
COORDS = [(_x,_y) for _x in xdom, _y in ydom]
XDOM, YDOM = first.(COORDS)[:], last.(COORDS)[:]

@inline function double_beam(k,r, w, a1::Int,a2::Int)
  # coeffs
  return a1*BWM.gaussianWave(k,r - SVector(-20.0, 3.5), w; abstol=1e-3) + a2*BWM.gaussianWave(k,r - SVector(-20.0, -13.0), w; abstol=1e-3)
end

@inline function wave_function(_freq::Float64, _kvec::SVector{2, Float64}, _width::Float64, _coeff::SVector{2, Int64})
  return BWM.boundaryWallWave(_freq * _kvec, (k,r) -> double_beam(k,r, _width, _coeff...), x, y, xm, ym, XDOM, YDOM, SIGMA, ds, rij, length(ds), N, BANDED, POTENTIAL_STRENGTH)
    
end


begin
wave = wave_function(0.6, KVEC, 2.4, SVector(1,1))

lines(abs2.(wave))
end

heatmap(xdom, ydom,abs2.(wave))

inputs = SVector.([1,0,1], [0,1,1])
waves = ThreadsX.map(gate -> wave_function(0.89, KVEC, 1.4, gate), inputs)
COLOR_RANGE = (0.0, maximum([maximum(abs2.(w)) for w in waves]))

# COLOR_RANGE = (0,5)

# see the multiplexer path
begin
idx_sensor = NDOM - 3
labels=["(1,0)", "(0,1)", "(1,1)"]
fig = Figure(theme=BWM.theme, size=(600, 800).*2, fontsize=20)
gl = fig[1,1] = GridLayout()
thresh = 0.5
for (j,w) in enumerate(waves)
    
  xp, yp, u, v = BWM.gradient(xdom, ydom, w)
  wave_reshaped = reshape(w, length(xdom), length(ydom))

  k = 2
  XP = xp[1:k:end, 1:k:end]
  YP = yp[1:k:end, 1:k:end]
  U  = u[1:k:end, 1:k:end]
  V  = v[1:k:end, 1:k:end]
  C  = abs2.(wave_reshaped)[1:k:end, 1:k:end]

  # C = angle.(wave)[1:k:end, 1:k:end]
  # ax = Axis(fig[1,1])
  ax,_=arrows(gl[j,1], Makie.Point2f.(XP, YP)[:],Makie.Point2f.(U,V)[:], normalize=true,arrowsize = 5, align=:center,lengthscale=1.0,color=C[:], colormap=BWM.cmap, )
  # heatmap!(ax, xdom, ydom, angle.(wave_reshaped), colormap=:twilight)
  for cen in CENTERS; 
    circ = BWM.createCircle(R, θ, SVector(cen));  
    lines!(ax, circ[1], circ[2], linestyle=:solid)
  end
  ax.aspect=DataAspect()
  # text!(ax, -15.0, 9.0; text=L"\omega=%$(FREQS[j])")
  Label(gl[j,0], labels[j], tellheight=false)
  hidedecorations!(ax, ticks=false, minorticks=false)
  ax,_=heatmap(gl[j,2], xdom, ydom, abs2.(w), colormap=BWM.cmap, colorrange=COLOR_RANGE)
  for cen in CENTERS; 
    circ = BWM.createCircle(R, θ, SVector(cen));  
    lines!(ax, circ[1], circ[2], linestyle=:solid)
  end
  ax.aspect=DataAspect()
  hidedecorations!(ax, ticks=false, minorticks=false)
  lines!(ax, [xdom[idx_sensor], xdom[idx_sensor]], [ydom[1], ydom[end]], linestyle=:dash, color=:blue)
  ax,_=lines(gl[j,3],abs2.(w)[idx_sensor,:], ydom,color=:blue)
  xlims!(ax, -0.1, 1.0)
  hidedecorations!(ax, ticks=false, minorticks=false)
  condition_met = maximum(abs2.(w[idx_sensor,:])) > thresh
  lines!(ax, [thresh, thresh],[ydom[1], ydom[end]], color=:red)
  Label(gl[j,4], "$(Int(condition_met))")
  # ax.aspect=DataAspect()
end
# Label(gl[0,0], "Inputs")
Colorbar(gl[4,1:3], vertical=false, colorrange=COLOR_RANGE, colormap=BWM.cmap, ticks=(collect(COLOR_RANGE), ["min", "max"]), tickalign=0, flipaxis=false, label="ΨΨ*")

# colsize!(gl,3, Aspect(2, 0.5))
# [rowsize!(gl,i,Auto(0.1)) for i in eachindex(gl.rowsizes)]
Label(gl[0,0], "In", tellwidth=false)
Label(gl[0,1], "Flow", tellwidth=false)
Label(gl[0,2], "Density", tellwidth=false)
Label(gl[0,3], "Sensor", tellwidth=false)
Label(gl[0,4], "Out", tellwidth=false)
colsize!(gl,3,Fixed(200))
[rowsize!(gl,i,Fixed(350)) for i in 1:3]

# rowsize!(gl,1,Fixed(100))
# rowsize!(gl,2,Auto(0.5))
# rowsize!(gl,3,Auto(0.1))
resize_to_layout!(fig)

# save("logic_gate.svg", fig)
fig

end