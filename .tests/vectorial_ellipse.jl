using GLMakie, StaticArrays, LinearAlgebra
using Meshes

include("../src/GeometryUtils.jl")
include("../src/BoundaryWall.jl")

R = 1.0
T = LinRange(-pi, pi, 1001)
ellipse = GeometryUtils.createEllipse(R, R, T, pi/2, (0.0,0.0))
x  = getindex(ellipse, 1)
y  = getindex(ellipse, 2)
xm = getindex(ellipse, 3)
ym = getindex(ellipse, 4)
ds = getindex(ellipse, 5)
rij= GeometryUtils.calcDistances(xm,ym)

function polarized_field(k::SVector{2, Float64}, x::Float64, y::Float64, ϕ::Float64)
  return exp(im*ϕ)*exp(-im * (k[1] * x + k[2] * y))
end

Ex0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarized_field(k, x, y, 0.0)
Ey0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarized_field(k, x, y, 0.0)
Ez0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarized_field(k, x, y, 0.0)

incident_waves = SVector(Ex0, Ey0, Ez0)

begin
GC.gc()
Nx, Ny = 150,150
xdom = LinRange(floor(minimum(x))-1,ceil(maximum(x))+1, Nx)
ydom = LinRange(floor(minimum(y))-1,ceil(maximum(y))+1, Ny)
# xdom = [centers[end][2]+2R, centers[end÷2][2], centers[1][2]-2R]
# ydom = [-1.0, 21.0]

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coords.(vertices(GRID))

XDOM, YDOM = first.(COORDS), last.(COORDS)

wave_vector = 10.1735*SVector(cosd(45),sind(45))

wave = BoundaryWall.boundaryWallVec(wave_vector, 
                                    xm, 
                                    ym, 
                                    XDOM, 
                                    YDOM, 
                                    -0.25im,
                                    rij,
                                    ds,
                                    length(ds),
                                    diagind(size(rij)...),
                                    incident_waves)

end

# begin 
#   fig = Figure()
#   # rs_h = IntervalSlider(fig[2, 1], range = LinRange(-100, 100, 1000),
#   #   startvalues = (-10.0, 10.0))
#   cr = Slider(fig[2,1], range=range(0, 100, 100), startvalue=5)

#   ax = Axis(fig[1,1])
#   hm = lift(cr.value) do int
#     # empty!(ax)
#   heatmap!(ax, xdom, ydom, abs2.(reshape(wave[1], Nx, Ny)), interpolate=true,colormap=:turbo, colorrange=(0,int))
#   lines!(ax, x,y, color=:black)
  
#   end
#   ax.aspect=DataAspect()
#   fig
  
# end
total_field = mapreduce(x->abs2.(x), +, wave)


# stokes = BoundaryWall.calcStokes.(wave[1], wave[2])

let
  GC.gc()
  # ψ = copy(real(reshape(getindex(wave,3), Nx, Ny)))
  ψ = map(i->reshape(abs.(i), Nx,Ny), wave)
  # crange = (-1,1) .* (max(maximum(ψ), -minimum(ψ)))
  
  fig = Figure(backgroundcolor=:white, fontsize=14)
  ga  = fig[1,1]=GridLayout()
  # ax = [Axis(ga[1,1], title=L"\mathrm{Re}[E_x]"),
  #       Axis(ga[1,2], title=L"\mathrm{Re}[E_y]"),
  #       Axis(ga[2,1], title=L"\mathrm{Re}[E_z ]"),
  #       Axis(ga[2,2], title=L"2\mathrm{Re}[E_xE_y^*]")]
  ax = [Axis(ga[1,1], title=L"|E_x|^2"),
        Axis(ga[1,2], title=L"|E_y|^2"),
        Axis(ga[2,1], title=L"|E_z|^2"),
        Axis(ga[2,2], title=L"|E|^2")]
  
  # hm=heatmap!(ax, xdom, ydom,ψ ,colormap=:balance, interpolate=false)
  
  contour!(ax[1], xdom, ydom, ψ[1], color=:black, levels=range(0.0, 3.0, 10))
  # contour!(ax[1], xdom, ydom, ψ[1], color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)
  
  contour!(ax[2], xdom, ydom, ψ[2], color=:black, levels=range(0.0, 3.0, 10))
  # contour!(ax[2], xdom, ydom, ψ[2], color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)
  
  contour!(ax[3], xdom, ydom, ψ[3], color=:black, levels=range(0.0, 2.0, 10))
  # contour!(ax[3], xdom, ydom, ψ[3], color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)
  
  
  contour!(ax[4], xdom, ydom, reshape(total_field, Nx, Ny), color=:black,levels=range(0,5,10))
  # contour!(ax[4], xdom, ydom, reshape(first.(stokes), Nx, Ny), color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)
  for ax in ax
  ax.backgroundcolor=:transparent
  ax.xticksmirrored=true
  ax.yticksmirrored=true
  ax.xtickalign=1
  ax.ytickalign=1
  ax.xgridvisible=false; ax.ygridvisible=false
  ax.xautolimitmargin=(0.025f0, 0.025f0)
  ax.yautolimitmargin=(0.05f0, 0.05f0)
  ax.xticks=xdom[1]:2:xdom[end]
  ax.yticks=ydom[1]:2:ydom[end]
  ax.xminorticksvisible=true
  ax.yminorticksvisible=true
  ax.xminortickalign=1
  ax.yminortickalign=1
  ax.xminorticks=IntervalsBetween(2)
  ax.yminorticks=IntervalsBetween(2)
  hidedecorations!(ax, ticks=false, minorticks=false)
  end
  # [lines!(ax, getindex(circ, 1), getindex(circ,2), color=Makie.wong_colors()[1], linewidth=2.0) for circ in circles, ax in ax]
  [lines!(ax, x,y) for ax in ax]
  [ax.aspect=DataAspect() for ax in ax]

  fig
end