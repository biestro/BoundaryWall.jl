# TODO: if process finds two peaks, continue with the first one, but insert the 
# second one in the original resonance vector

using BoundaryWall
using StaticArrays
using CircularArrays
using GLMakie
using Peaks
using ProgressBars

# geometry
N = 500
θ = LinRange(0,2pi,N)
R = 1.0
x,y,xm,ym,ds = BoundaryWall.createCircle(R, θ, SVector(0.0, 0.0))
rij = BoundaryWall.calcDistances(xm,ym)

rmid = SVector.(xm,ym)
rpos = CircularArray(SVector.(x,y))
# calculate eigenstates numerically

# initial scan
begin # first step
  wave_numbers = LinRange(0.001, 20.0, 100)
  SIGMA = -0.25im

  Tij = zeros(length(wave_numbers))

  # calculate spectrum
  for (_i,_k) in ProgressBar(enumerate(wave_numbers))
    Tij[_i] = sum(abs.(inv(BoundaryWall.calcMmatrix(_k, SIGMA,rij,ds, 1, rmid,rpos, N))))
  end

  peaks, vals = findmaxima(Tij)

  resonances = wave_numbers[peaks]
  window = minimum(diff(resonances))/2 # minimum window between resonances
end

lines(wave_numbers[1:end-1],diff(Tij))


resonanceFunction(_k::Float64,_band::Int64) = sum(abs.(inv(BoundaryWall.calcMmatrix(_k, SIGMA,rij,ds, _band, rmid,rpos, N))))


# second step
begin
  GC.gc()
  subsearch_density = 25
  
  _tij_subsearch = zeros(subsearch_density)

  _b = 5 # integration band

  MAX_ITER = 20
  window = minimum(diff(resonances))/2 # minimum window between resonances
  resonances_subsearch = Float64[]
  ϵ = 1e-5
  for (_i,_r) in enumerate(resonances)
  println("Current step: $_i / $(length(resonances))")
  _window = window # original specturm window
  _r_prev = _r     # current resonance
  _delta_r = 1.0
    for _iter in 1:MAX_ITER  
    
      _k_subsearch = LinRange(_r - _window, _r + _window, subsearch_density)
      for (_i,_k) in ProgressBar(enumerate(_k_subsearch))
        _tij_subsearch[_i] = resonanceFunction(_k,_b)
      end
      _peaks, _ = findmaxima(_tij_subsearch)
      
      # try to find a peak
      try 
        _k_subsearch[_peaks][1]
      catch e
        @warn ("NO PEAK FOUND!! Reverting to previous iteration's peak....")
        push!(resonances_subsearch, _r_prev)
        break
      end

      length(_k_subsearch[_peaks]) > 1 ? @warn("More than one peak found around $_r_prev") : nothing;
      _r_new = _k_subsearch[_peaks][1]

      # update 
      _delta_r = abs(_r_prev - _r_new)
      println("Δk = $_delta_r")
      if _delta_r< ϵ
        println("TOLERANCE REACHED! Finishing process...")
        push!(resonances_subsearch,_r_new)
        break
      end
      _window = _window/2 # halve the window (reduce the search)
      _r_prev = _r_new
      if _iter % MAX_ITER == 0
        println("MAX ITERATIONS REACHED! Finishing process...")
        push!(resonances_subsearch, _r_new)
      end
    end
  end
end

scatter(abs.(resonances - resonances_subsearch))

begin using Meshes
nx,ny = 100,100
x0, xf = (-2.0R, 2.0R)
y0, yf = (-2.0R, 2.0R)
xdom = LinRange(x0, xf, nx)
ydom = LinRange(y0, yf, ny)
GRID = RectilinearGrid(xdom, ydom)
MESH = SimpleMesh(vertices(GRID), GRID.topology)
COORDS = SVector.(coords.(vertices(MESH)))
XDOM = (first.(COORDS))
YDOM = (last.(COORDS))
end
wave = BoundaryWall.boundaryWallWave(resonances_subsearch[end]*SVector(1.0,0.0), planeWave, x,y,xm,ym,XDOM, YDOM, SIGMA, ds, rij, length(ds), N, 2, Inf)

f,ax=heatmap(xdom, ydom, reshape(angle.(wave), nx, ny), colormap=:twilight);ax.aspect=DataAspect();f
f,ax=heatmap(xdom, ydom, reshape(abs.(wave), nx, ny), colormap=:dense);ax.aspect=DataAspect();f