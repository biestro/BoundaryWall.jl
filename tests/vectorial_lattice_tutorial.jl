
# Tutorial for generating lattices with small 


using GLMakie, StaticArrays, LinearAlgebra
using Meshes
using StatsBase: sample
using Distributions: Uniform

include("../src/GeometryUtils.jl")

begin
nx, ny  = (10, 10)
origin  = (0.0, 0.0)
n_holes = nx*ny ÷ 3
grid    = GeometryUtils.HexagonalGrid(origin, 1.0)
R       = grid.a / 4.5
N       = 10
TH      = LinRange(-pi,pi,N)

centers = GeometryUtils.buildGrid(grid, nx, ny)
indices_to_remove = sample(LinearIndices((1:nx, 1:ny)), n_holes, replace=false)
Δr = maximum(centers) ./ 2
map!(c->c .- Δr, centers, centers)
# deleteat!(centers, sort(indices_to_remove))

# circles = [GeometryUtils.createEllipse(R, R/2, TH,rand()*2pi, c) for c in centers[:]]
circles = [GeometryUtils.createCircle(R,TH, c) for c in centers[:]]

x  = vcat(getindex.(circles, 1)...)
y  = vcat(getindex.(circles, 2)...)
xm = vcat(getindex.(circles, 3)...)
ym = vcat(getindex.(circles, 4)...)
ds = vcat(getindex.(circles, 5)...)
rij= GeometryUtils.calcDistances(xm,ym)
end


let
  f,ax = scatter(centers)
  lines!(ax, x,y)
  ax.aspect=DataAspect()
  f
end

include("../src/BoundaryWall.jl")

function polarized_field(k::SVector{2, Float64}, x::Float64, y::Float64, ϕ::Float64)
  return exp(im*ϕ)*exp(-im * (k[1] * x + k[2] * y))
end

Ex0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarized_field(k, x, y, 0.0)
Ey0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarized_field(k, x, y, pi/2)
Ez0(k::SVector{2, Float64}, x::Float64, y::Float64) = polarized_field(k, x, y, 0.0)

incident_waves = SVector(Ex0, Ey0, Ez0)


begin
GC.gc()
Nx, Ny = 10,10
xdom = LinRange(floor(minimum(x))-1,ceil(maximum(x))+1, Nx)
ydom = LinRange(floor(minimum(y))-1,ceil(maximum(y))+1, Ny)
# xdom = [centers[end][2]+2R, centers[end÷2][2], centers[1][2]-2R]
# ydom = [-1.0, 21.0]

GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coordinates.(vertices(GRID))

XDOM, YDOM = first.(COORDS), last.(COORDS)

wave_vector = 5.0*SVector(-1.0,0.0)
# wave = BoundaryWall.boundaryWallWave(wave_vector, xm,ym,XDOM, YDOM, -0.25im,rij,ds,length(ds), diagind(size(rij)...),20.0)
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
                                    incident_waves, 10.0)



end


# heatmap(xdom, ydom, abs2.(reshape(getindex(wave,2),Nx,Ny)), colorrange=(0,3),colormap=:balance, interpolate=true)

stokes = BoundaryWall.calcStokes.(wave[1], wave[2])
GLMakie.activate!()
let
GC.gc()
# ψ = copy(real(reshape(getindex(wave,3), Nx, Ny)))
ψ = map(i->reshape(abs2.(i), Nx,Ny), wave)
# crange = (-1,1) .* (max(maximum(ψ), -minimum(ψ)))

fig = Figure(backgroundcolor=:white, fontsize=14)
ga  = fig[1,1]=GridLayout()
ax = [Axis(ga[1,1], title=L"\mathrm{Re}[E_x]"),
      Axis(ga[1,2], title=L"\mathrm{Re}[E_y]"),
      Axis(ga[2,1], title=L"\mathrm{Re}[E_z ]"),
      Axis(ga[2,2], title=L"2\mathrm{Re}[E_xE_y^*]")]

# hm=heatmap!(ax, xdom, ydom,ψ ,colormap=:balance, interpolate=false)

contour!(ax[1], xdom, ydom, ψ[1], color=:black, levels=range(0.1, 1.0, 6))
contour!(ax[1], xdom, ydom, ψ[1], color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)

contour!(ax[2], xdom, ydom, ψ[2], color=:black, levels=range(0.1, 3.0, 6))
contour!(ax[2], xdom, ydom, ψ[2], color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)

contour!(ax[3], xdom, ydom, ψ[3], color=:black, levels=range(0.1, 3.0, 6))
contour!(ax[3], xdom, ydom, ψ[3], color=:black, levels=range(-3.0,-0.1, 6), linestyle=:dot)


contour!(ax[4], xdom, ydom, reshape(getindex.(stokes,3), Nx, Ny), color=:black,levels=range(0,2,5))
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
# colsize!(ga, 1, Aspect(1, 1))

# xlims!(ax, xdom[1]-1, xdom[end]+1)
# ylims!(ax, ydom[1]-1, ydom[end]+1)
# scatter!(ax, centers, color=:black)
[lines!(ax, getindex(circ, 1), getindex(circ,2), color=Makie.wong_colors()[1], linewidth=2.0) for circ in circles, ax in ax]
# [scatter!(ax, centers, color=Makie.wong_colors()[1], markersize=5) for ax in ax]
# Colorbar(ga[1,2], colorrange=(0,1), colormap=Reverse(:grayC))
[ax.aspect=DataAspect() for ax in ax]

# cb = Colorbar(ga[1,2], hm)
# rowsize!(fig.layout, 1, Aspect(1, 1))
# rowsize!(fig.layout, 1, ax.scene.px_area[].widths[2])

# cb.height=Relative()
# rowsize!(ga[1,1], 1, Aspect(1,1))
# resize_to_layout!(fig)
# save("vector_lattice_circular.eps", fig)
fig
end

