using DynamicalBilliards
using GLMakie
const SV = SVector{2}

@inline function cart2par(x::Float64,y::Float64)
  # convert cartesian to parabolic coordinates
  # convention to have the positive x axis as the 
  # branch cut
  z = -x + 1im * y; # convert to imaginary unit
  z = sqrt(2*z);
  _xi = imag(z);
  _eta = real(z);
  return _xi, _eta
end 

@inline function par2cart(_xi::Float64, _eta::Float64)
  X = @. 1/2 * (_xi^2 - _eta^2)
  Y = @. _xi * _eta
  return X,Y
end

N = 100
XI0 = 3.0; ETA0 = 2.0;
xi = LinRange(0, XI0, N)
eta = LinRange(0, ETA0, N)

z1 = [par2cart(_r...) for _r in SVector.(xi,ETA0)]
z2 = [par2cart(_r...) for _r in SVector.(XI0,eta)]

# transform operation
function dilate(_x::Float64, _y::Float64, _dil::SVector{2,Float64},_c::SVector{2,Float64})
  # dilates coordinates by _dil amount
  # _c is the center
  _x_dilated = (_x - _c[1]) * _dil[1] + _c[1]
  _y_dilated = (_y - _c[2]) * _dil[2] + _c[2]

  return _x_dilated, _y_dilated
end

dilation_fwd = SVector(1.2, 1.0)
dilation_bwd = SVector(0.8, 1.0)

z1_dilated_fwd = [dilate(_r[1], _r[2], dilation_fwd , SVector(first(par2cart(XI0,ETA0)), 0.0)) for _r in z1]
z1_dilated_bwd = [dilate(_r[1], _r[2], dilation_bwd , SVector(first(par2cart(XI0,ETA0)), 0.0)) for _r in z1]

original_parabola = [z1; reverse(z2)]
dilated_parabola_fwd  = [z1_dilated_fwd; reverse(z2)]
dilated_parabola_bwd  = [z1_dilated_bwd; reverse(z2)]

begin
  _y = last.(original_parabola)
  _y = -_y
  _x = reverse(first.(original_parabola))
  confocal_parabolic = [tuple.(_x, _y);original_parabola]

  _y = last.(dilated_parabola_fwd)
  _y = -_y
  _x = reverse(first.(dilated_parabola_fwd))
  chaotic_parabolic_fwd = [tuple.(_x, _y); dilated_parabola_fwd]

  _y = last.(dilated_parabola_bwd)
  _y = -_y
  _x = reverse(first.(dilated_parabola_bwd))
  chaotic_parabolic_bwd = [tuple.(_x, _y); dilated_parabola_bwd]

end




# boundary wall method
import BoundaryWall as BWM

# chaotic forwards
begin
N = 400
xc_fwd, yc_fwd   = first.(chaotic_parabolic_fwd), last.(chaotic_parabolic_fwd)
xc_fwd, yc_fwd   = BWM.divideCurve(xc_fwd,yc_fwd, N)
xmc_fwd, ymc_fwd = BWM.calcMidpoints(xc_fwd,yc_fwd)
dsc_fwd      = BWM.calcArcLength(xc_fwd,yc_fwd)
rij_fwd     = BWM.calcDistances(xmc_fwd, ymc_fwd)
end

# chaotic backwards
begin
N = 400
xc_bwd, yc_bwd   = first.(chaotic_parabolic_bwd), last.(chaotic_parabolic_bwd)
xc_bwd, yc_bwd   = BWM.divideCurve(xc_bwd,yc_bwd, N)
xmc_bwd, ymc_bwd = BWM.calcMidpoints(xc_bwd,yc_bwd)
dsc_bwd      = BWM.calcArcLength(xc_bwd,yc_bwd)
rij_bwd     = BWM.calcDistances(xmc_bwd, ymc_bwd)
end

# confocal
begin
N = 400
x, y   = first.(confocal_parabolic), last.(confocal_parabolic)
x, y   = BWM.divideCurve(x,y, N)
xm, ym = BWM.calcMidpoints(x,y)
ds     = BWM.calcArcLength(x,y)
rij    = BWM.calcDistances(xm, ym)
end



let fig = Figure(fontsize=20)
  ga = fig[1,1] = GridLayout()
  ax = Axis(ga[1,1], xlabel=L"y", ylabel=L"x", ylabelrotation=0, yreversed=true)
  # ax,_=lines(fig[1,1],original_parabola, color=eachindex(original_parabola))
  ln1 = lines!(ax, reverse.(confocal_parabolic), color=:black,  linewidth=2.0, label=L"D_x\;>1")
  ln2 = lines!(ax, reverse.(chaotic_parabolic_fwd), color=:red, linewidth=2.0, linestyle=:dash, label=L"D_x\;=1")
  ln3 = lines!(ax, reverse.(chaotic_parabolic_bwd), color=:red, linewidth=2.0, linestyle=:dot, label=L"D_x\;<1")
  hidedecorations!(ax, label=false)
  Legend(ga[1,2], ax)
  # colgap!(ga, -200.0)
  ax.aspect=DataAspect()
  fig
end


begin using Meshes
  nx,ny = 100,100
  x0, xf = (minimum(xc_fwd)-1.0, maximum(xc_fwd)+1.0)
  y0, yf = (minimum(yc_fwd)-1.0, maximum(yc_fwd)+1.0)
  xdom = LinRange(x0, xf, nx)
  ydom = LinRange(y0, yf, ny)
  GRID = RectilinearGrid(xdom, ydom)
  MESH = SimpleMesh(vertices(GRID), GRID.topology)
  COORDS = SVector.(coords.(vertices(MESH)))
  XDOM = (first.(COORDS))
  YDOM = (last.(COORDS))
end



using CircularArrays, ProgressBars, Peaks
rmid   = SVector.(xm, ym)
rpos   = CircularArray(SVector.(x, y))
rmid_c_fwd = SVector.(xmc_fwd,ymc_fwd)
rpos_c_fwd = CircularArray(SVector.(xc_fwd,yc_fwd))
rmid_c_bwd = SVector.(xmc_bwd,ymc_bwd)
rpos_c_bwd = CircularArray(SVector.(xc_bwd,yc_bwd))

resonanceFunction(_k::Float64,rij::Matrix{Float64},ds::Vector{Float64}, _band::Int64, rmid::Vector,rpos::CircularVector{SVector{2, Float64}, Vector{SVector{2, Float64}}}) = sum(abs.(inv(BWM.calcMmatrix(_k, SIGMA, rij,ds, _band, rmid,rpos, N))))

# chaotic
begin # first sweep
  wave_numbers = LinRange(0.05, 1.5, 300)
  SIGMA = -0.25im
  
  Tij_confocal = zeros(length(wave_numbers))
  Tij_chaotic_fwd = zeros(length(wave_numbers))
  Tij_chaotic_bwd = zeros(length(wave_numbers))
  
  # calculate spectrum
  for (_i,_k) in ProgressBar(enumerate(wave_numbers))
    Tij_confocal[_i]    = resonanceFunction(_k,rij,ds, 5, rmid,rpos)
    Tij_chaotic_fwd[_i] = resonanceFunction(_k,rij_fwd,dsc_fwd, 5, rmid_c_fwd,rpos_c_fwd)
    Tij_chaotic_bwd[_i] = resonanceFunction(_k,rij_bwd,dsc_bwd, 5, rmid_c_bwd,rpos_c_bwd)
  end
end

begin
  peaks_confocal, vals_confocal = findmaxima(Tij_confocal)
  peaks_chaotic_fwd, vals_chaotic_fwd = findmaxima(Tij_chaotic_fwd)
  peaks_chaotic_bwd, vals_chaotic_bwd = findmaxima(Tij_chaotic_bwd)

  resonances_confocal = wave_numbers[peaks_confocal]
  resonances_chaotic_fwd = wave_numbers[peaks_chaotic_fwd]
  resonances_chaotic_bwd = wave_numbers[peaks_chaotic_bwd]
  # window = minimum(diff(resonances))/2 # minimum window between resonances
end

GC.gc()


let 
  COLORMAP = :dense
  fig = Figure(fontsize=20)
  ax = [Axis(fig[1,1]),Axis(fig[2,1]),Axis(fig[3,1])]
  heatmap!(ax[1], wave_numbers, [0],hcat(Tij_chaotic_fwd...)', colorscale=log10, colormap=COLORMAP)
  heatmap!(ax[2], wave_numbers, [0],hcat(Tij_confocal...)'   , colorscale=log10, colormap=COLORMAP)
  heatmap!(ax[3], wave_numbers, [0],hcat(Tij_chaotic_bwd...)', colorscale=log10, colormap=COLORMAP)
  hidexdecorations!(ax[2])
  hidexdecorations!(ax[1])
  [hideydecorations!(ax) for ax in ax]
  ax[3].xlabel= L"k^2"
  [text!(ax, 1.01, 0.20, text=_t) for (ax,_t) in zip(ax, ["a)", "b)", "c)"])]
  fig
end


# refine resonances
# TODO: update _resonances when two are found, refineResonances!
function refineResonances(_resonances::Vector,_rij::Matrix{Float64},_ds::Vector{Float64},_rmid::Vector,_rpos::CircularVector{SVector{2, Float64}, Vector{SVector{2, Float64}}})

  GC.gc()
  subsearch_density = 55
  
  _tij_subsearch = zeros(subsearch_density)

  _b = 5 # integration band

  MAX_ITER = 20
  window = minimum(diff(_resonances))/2 # minimum window between _resonances
  _resonances_subsearch = Float64[]
  ϵ = 1e-5
  println("Beginning refinement process")
  for (_i,_r) in enumerate(_resonances)
  println("Current step: $_i / $(length(_resonances))")
  _window = window # original specturm window
  _r_prev = _r     # current resonance
  _delta_r = 1.0
    for _iter in 1:MAX_ITER  
    
      _k_subsearch = LinRange(_r - _window, _r + _window, subsearch_density)
      for (_i,_k) in ProgressBar(enumerate(_k_subsearch))
        _tij_subsearch[_i] = resonanceFunction(_k,_rij,_ds, _b, _rmid, _rpos)
      end
      _peaks, _ = findmaxima(_tij_subsearch)
      
      # try to find a peak
      try 
        _k_subsearch[_peaks][1]
      catch e
        @warn ("NO PEAK FOUND!! Reverting to previous iteration's peak....")
        push!(_resonances_subsearch, _r_prev)
        break
      end

      length(_k_subsearch[_peaks]) > 1 ? @warn("More than one peak found around $_r_prev") : nothing;
      _r_new = _k_subsearch[_peaks][1]

      # update 
      _delta_r = abs(_r_prev - _r_new)
      println("Δk = $_delta_r")
      if _delta_r< ϵ
        println("TOLERANCE REACHED! Finishing process...")
        push!(_resonances_subsearch,_r_new)
        break
      end
      _window = _window/2 # halve the window (reduce the search)
      _r_prev = _r_new
      if _iter % MAX_ITER == 0
        println("MAX ITERATIONS REACHED! Finishing process...")
        push!(_resonances_subsearch, _r_new)
      end
    end
  end
  return _resonances_subsearch
end


# refine resonances around a ball
_refined_confocal = refineResonances(resonances_confocal, rij, ds, rmid, rpos)
_refined_fwd = refineResonances(resonances_chaotic_fwd, rij_fwd, dsc_fwd, rmid_c_fwd, rpos_c_fwd)
_refined_bwd = refineResonances(resonances_chaotic_bwd, rij_bwd, dsc_bwd, rmid_c_bwd, rpos_c_bwd)

# get wavenumber to evaluate wave
_wavenum = getindex.([_refined_bwd,_refined_confocal, _refined_fwd ], 5)

# get waves
wave_confocal = BWM.boundaryWallWave(_wavenum[2]*SVector(cosd(45), sind(45)), BWM.planeWave, 
        x,y, xm, ym,XDOM, YDOM, -0.25im, ds, rij, length(ds), length(ds),1, Inf)

wave_chaotic_fwd = BWM.boundaryWallWave(_wavenum[3]*SVector(cosd(45), sind(45)), BWM.planeWave, 
        xc_fwd,yc_fwd, xmc_fwd, ymc_fwd,XDOM, YDOM, -0.25im, dsc_fwd, rij_fwd, length(dsc_fwd), length(dsc_fwd), 1, Inf)

wave_chaotic_bwd = BWM.boundaryWallWave(_wavenum[1]*SVector(cosd(45), sind(45)), BWM.planeWave, 
        xc_bwd,yc_bwd, xmc_bwd, ymc_bwd,XDOM, YDOM, -0.25im, dsc_bwd, rij_bwd, length(dsc_bwd), length(dsc_bwd), 1, Inf)

let fig = Figure()
  c_max = 0.5
  ax = [Axis(fig[1,1]), Axis(fig[1,2]), Axis(fig[1,3]),Axis(fig[2,1]), Axis(fig[2,2]), Axis(fig[2,3])]
  heatmap!(ax[1], ydom, xdom, reshape(abs2.(wave_chaotic_bwd), nx, ny)',interpolate=true, colormap=:dense, colorrange=(0,c_max))
  heatmap!(ax[2], ydom, xdom, reshape(abs2.(wave_confocal), nx, ny)',   interpolate=true, colormap=:dense, colorrange=(0,c_max))
  heatmap!(ax[3], ydom, xdom, reshape(abs2.(wave_chaotic_fwd), nx, ny)',interpolate=true, colormap=:dense, colorrange=(0,c_max))

  heatmap!(ax[4], ydom, xdom, reshape(angle.(wave_chaotic_bwd), nx, ny)', colormap=:twilight)
  heatmap!(ax[5], ydom, xdom, reshape(angle.(wave_confocal), nx, ny)',    colormap=:twilight)
  heatmap!(ax[6], ydom, xdom, reshape(angle.(wave_chaotic_fwd), nx, ny)', colormap=:twilight)
  [ax.aspect=DataAspect() for ax in ax]
  [ax.yreversed=true for ax in ax]
  [hidedecorations!(ax) for ax in ax]
  [lines!(ax, last.(_r.data), first.(_r.data), color=:black) for (ax,_r) in zip(ax, [rpos_c_bwd,rpos,rpos_c_fwd,rpos_c_bwd,rpos,rpos_c_fwd])]
  # [lines!(_ax, first.(_r), last.(_r), color=:black) for (_ax,_x,_y) in zip(ax, [rpos_c_bwd,rpos,rpos_c_fwd,rpos_c_bwd,rpos,rpos_c_fwd])]
  fig
end

T = [abs.(inv(BWM.calcMmatrix(_k,-0.25im,_rij,_ds, 5,_rmid,_rpos, N))) for (_k,_rij,_ds,_rmid,_rpos) in 
     zip(_wavenum, [rij_bwd, rij, rij_fwd], [dsc_bwd, ds, dsc_fwd], [rmid_c_bwd, rmid, rmid_c_fwd], [rpos_c_bwd, rpos, rpos_c_fwd])]

let fig = Figure()
  ax = [Axis(fig[1,1]), Axis(fig[1,2]), Axis(fig[1,3])]
  [heatmap!(ax, t, colormap=:grayC, colorscale=log10) for (ax,t) in zip(ax,T)]
  [ax.aspect=DataAspect() for ax in ax]
  [ax.yreversed=true for ax in ax]
  # lines!(ax, xc_fwd,yc_fwd, color=:black)
  fig
end
