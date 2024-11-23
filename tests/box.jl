# This file tries to find the T-matrix for multiple circles, given they all have
# the same individual T-matrix. Preliminary tests btw, nothing serious

using GLMakie, LinearAlgebra, StaticArrays
include("../src/BoundaryWall.jl")
include("../src/BWM_geometry.jl")

using Distributions: Uniform

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
TH   = 180                       # incident angle
# KVEC = zero/R * SVector(cosd(TH), sind(TH))
KVEC = 2.5 * SVector(cosd(TH), sind(TH))
# circle 1
end

begin
  WINDOW = 5.0
  GAMMA_MAX = 30.0
    # show model
  STEP = 2R + R/2 # diameter  + constant
  N_CIRCLES = 11
  RANGES = -(N_CIRCLES-1)*STEP/2:STEP:(N_CIRCLES-1)*STEP/2
  N_STEPS =length(RANGES)
  CENTERS = vec([(i, j) for i in RANGES, j in RANGES])

  CENTERS  = vec([(i, j) for i in RANGES, j in RANGES])
  # STRENGTH = ones(length(RANGES), length(RANGES));
  STRENGTH = rand(Uniform(0.0, GAMMA_MAX),length(RANGES), length(RANGES)) ;
  # STRENGTH[N_CIRCLES÷2+1, 1:end] .= 0.0
  # STRENGTH[4:end,N_CIRCLES÷2+1] .= 0.0
  # STRENGTH[diagind(length(RANGES), length(RANGES))] .= 0.5
  # STRENGTH[1:4, N_CIRCLES ÷ 2+2] .= 0.0
  # STRENGTH[2, N_CIRCLES ÷ 2+2] = 0
  # STRENGTH[3, N_CIRCLES ÷ 2+2] = 0
  # STRENGTH[4, N_CIRCLES ÷ 2+2] = 0
  reverse!(STRENGTH, dims=1)
  #STRENGTH=STRENGTH[:] .* 10.0e10
  STRENGTH=STRENGTH[:]
  INDICES  = sort([N_STEPS-5,
                  2N_STEPS-5,
                  3N_STEPS-5,
                  4N_STEPS-5,
                  5N_STEPS-5,
                  6N_STEPS-5, 6N_STEPS-4,6N_STEPS-3,6N_STEPS-6,6N_STEPS-7,
                  7N_STEPS-3, 7N_STEPS-7,
                  8N_STEPS-3, 8N_STEPS-7,
                  9N_STEPS-3, 9N_STEPS-7,
                 10N_STEPS-3,10N_STEPS-7,
                 11N_STEPS-3,11N_STEPS-7])

  # STRENGTH[INDICES] .= rand(length(INDICES))

  POTENTIAL_STRENGTH = repeat(STRENGTH, inner=N)
  CIRCLES = [BoundaryWall.createCircle(R, θ, SVector(cen)) for cen in CENTERS]
  x = vcat(getindex.(CIRCLES, 1)...)
  y = vcat(getindex.(CIRCLES, 2)...)
  xm= vcat(getindex.(CIRCLES, 3)...)
  ym= vcat(getindex.(CIRCLES, 4)...)
  ds= vcat(getindex.(CIRCLES, 5)...)
  rij = BoundaryWall.calcDistances(xm,ym)
  fig, ax = scatter(xm, ym, color=POTENTIAL_STRENGTH, markersize=5, colormap=:viridis)
  # text!(ax, CENTERS; text=string.(eachindex(CENTERS)))
  # text!(ax, CENTERS; text=string.(STRENGTH))
  quiver!(ax,[-3.0],[5.0], [-cosd(TH)],[-sind(TH)])
  ax.aspect=DataAspect()


  # creating box
  nbox = 100
  # NBOX = 229 # odd
  # # top box

  # L_END = RANGES[end]+2R
  # L_BEG = RANGES[1]-2R

  # x_top_box = [ones(nbox)*(L_BEG);LinRange(L_BEG,L_END, nbox);ones(nbox)*(RANGES[end]+2R)]
  # y_top_box = [LinRange(WINDOW/2, RANGES[end]+2R, nbox);ones(nbox) *(RANGES[end]+2R) ;LinRange(RANGES[end]+2R,WINDOW/2, nbox)]

  # # bottom right box
  # x_bot_right_box = [ones(nbox)*(L_END);LinRange(RANGES[end]+2R,WINDOW/2, nbox)]
  # y_bot_right_box = [LinRange(-WINDOW/2, RANGES[1]-2R,nbox);ones(nbox)*(RANGES[1]-2R)]

  # x_bot_left_box = -reverse(x_bot_right_box)
  # y_bot_left_box = reverse(y_bot_right_box)

  # XTOP, YTOP = GeometryUtils.divideResonatorCurve(x_top_box,y_top_box, NBOX)
  # # XTOP, YTOP = GeometryUtils.divideResonatorCurve(x_top_box,y_top_box, NBOX)
  # CORNER_R =  tuple.(reverse(XTOP[XTOP .> WINDOW/2]), -reverse(YTOP[XTOP .> WINDOW/2]))
  # CORNER_L =  tuple.(reverse(XTOP[XTOP .< -WINDOW/2]), -reverse(YTOP[XTOP .< -WINDOW/2]))

  # TOP_MID      = GeometryUtils.calcMidpoints(XTOP, YTOP)
  # DS_TOP       = GeometryUtils.calcArcLength(XTOP, YTOP)
  # CORNER_R_MID = GeometryUtils.calcMidpoints(first.(CORNER_R), last.(CORNER_R))
  # CORNER_R_DS  = GeometryUtils.calcArcLength(first.(CORNER_R), last.(CORNER_R))
  # CORNER_L_MID = GeometryUtils.calcMidpoints(first.(CORNER_L), last.(CORNER_L))
  # CORNER_L_DS  = GeometryUtils.calcArcLength(first.(CORNER_L), last.(CORNER_L))

  XBOT = ones(nbox)*(L_BEG)
  XTOP = ones(nbox)*(L_BEG);
  
  YBOT   = collect(LinRange(RANGES[1]-2R,-WINDOW/2,nbox));
  YTOP   = collect(LinRange(WINDOW/2, RANGES[end]+2R, nbox));

  TOP_MID  = GeometryUtils.calcMidpoints(XTOP, YTOP)
  BOT_MID  = GeometryUtils.calcMidpoints(XBOT, YBOT)


  DS_TOP = GeometryUtils.calcArcLength(XTOP, YTOP)
  DS = [DS_TOP; DS_TOP]

  BOX = [tuple.(XBOT, YBOT); tuple.(XTOP, YTOP)]
  MID = [tuple.(BOT_MID...); tuple.(TOP_MID...)]
  RIJ = GeometryUtils.calcDistances(first.(MID), last.(MID))


  # scatter!(BOX,color=eachindex(x_top_box))
  scatter!(MID,color=eachindex(x_top_box))

  # scatter!(CORNER_R, color=eachindex(x_top_box))
  # scatter!(CORNER_L, color=eachindex(x_top_box))

  fig
end

using Meshes

Nx, Ny = 150,150
xdom = LinRange(-17.0, 17.0, Nx)
ydom = LinRange(-17.0, 17.0, Ny)
GRID = RectilinearGrid(xdom, ydom)
# MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = coords.(vertices(GRID))

XDOM, YDOM = first.(COORDS), last.(COORDS)

wave = BoundaryWall.boundaryWallWave(KVEC, first.(MID), last.(MID), XDOM,YDOM, SIGMA, RIJ, DS, length(DS),diagind(size(RIJ)...))

let fig = Figure()
ax = Axis(fig[1,1])
viz!(ax,GRID, color=abs2.(wave), colorscheme=:dense)
ax.aspect=DataAspect()
scatter!(ax, MID, color=:black)
fig
end