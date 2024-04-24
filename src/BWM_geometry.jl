"""
Developed by Alberto RB for use in Billiard problems.
If you use this, thank me at least!!

Geometry utils for creating billiards, meshes, and more.
"""

# module BWM_geometry
import AbaqusReader, Meshes
using LinearAlgebra
using StaticArrays: SVector
using HOHQMesh, Meshes
# using Makie

# export createCircle, createEllipse, createConfocalBilliard, calcArcLength, calcDistances, calcMidpoints, divideCurve, divideResonatorCurve

"""
Redefinition, don't remember what for, but do not delete it
"""
Base.getindex(v::Meshes.Point2, i::Int64) = v.coords[i]


"""
  createCircle(r, THETA, SHIFT)

Creates a circle
"""
function createCircle(r::Float64, θ::LinRange, SHIFT::SVector{2, Float64})
  x,y    = r * cos.(θ) .+ SHIFT[1], r * sin.(θ) .+ SHIFT[2]
  xm, ym = BWM_geometry.calcMidpoints(x,y)
  ds     = BWM_geometry.calcArcLength(x,y)
  return x,y,xm,ym,ds
end

"""
  createEllipse(r, THETA, SHIFT)

Creates an ellipse
"""
function createEllipse(r::Float64, μ::Float64, θ::LinRange, ϕ::Float64, SHIFT::SVector{2, Float64})
  x,y    = r * cos.(θ), μ * sin.(θ)
  u,v    = @. x * cos(ϕ) - y*sin(ϕ), x*sin(ϕ) + y*cos(ϕ) # rotate
  u = u .+ SHIFT[1]
  v = v .+ SHIFT[2]
  u,v    = BWM_geometry.divideResonatorCurve(u,v,length(u))
  xm, ym = BWM_geometry.calcMidpoints(u,v)
  ds     = BWM_geometry.calcArcLength(u,v)
  return u,v,xm,ym,ds
end




"""
  calcDistances(a, b)

Computes distance matrix for |r - r'|

# Arguments
- `a::Vector`: r
- `b::Vector`: r'
"""
function calcDistances(a::Vector, b::Vector)
  return @. hypot(a - a', b - b')
end

"""
  confocalBilliard(sigma0, tau0, numSegments)

Constructs a confocal parabolic billiard formed by two parabolic contours given
by σ0 and τ0.

# Arguments
- `sigma0::Float64`: parabolic coordinate 1
- `tau0::Float64`: parabolic coordinate 2
- `numSegments::Int`: number of segment discretization 
"""
function confocalBilliard(sigma0::Float64, tau0::Float64, numSegments::Int)
  # return parabolic confocal billiard
  sigma = sigma0;
  tau = tau0 * (4/numSegments) * (0:numSegments/4); 

  x = @. (tau^2 - sigma^2)/2;
  y = sigma*tau;

  sigma = sigma0 * (4/numSegments) * (0:numSegments/4);
  tau = tau0;

  y1 = sigma*tau;
  x1 = @. (tau^2 - sigma^2)/2;

  x = [x; reverse(x1)[2:end]]
  y = [y; reverse(y1)[2:end]]

  x = [x;reverse(x)[2:end]]
  y = [y; -y[2:end]]

return -x,y  # flip false midpoint
end


"""
Builds an assymetric parabolic resonator with input channels. Needed for the
midpoints in the BWM, returned in `BWM_geometry.createResonator`.

# Arguments
- `XI0::Float64`: parabolic coordinate 1
- `ETA0::Float64`: parabolic coordinate 2
- `dXI::Float64`: upper channel/arm opening
- `dETA::Float64`: lower channel/arm opening
- `channelLength::Float64`: channel/arm length (both channels)
- `numSegments::Int`: number of segment discretization 

# Returns
- `x`: x coord
- `y`: y coord
- `xm`: x midpoint
- `ym`: y midpoint
- `rij_mid`: distance matrix |r - r'|
- `ds`: arc lengths
"""
function buildAsymmetricResonator(XI0::Float64, ETA0::Float64; dXI::Float64=0.5, dETA::Float64=0.5,channelLength::Float64=0.5, nSegments::Int=100)
  xi = LinRange(0, XI0, nSegments)
  eta = ETA0
  x,y = par2cart(xi, eta)
  # x = @. 1/2 * (xi^2 - eta^2)
  # y = @. xi*eta

  @assert (dXI < XI0) ? true : throw(AssertionError("dXi cannot be greater than ETA0 "))

  # eta = LinRange(0, (ETA0-dEta), nSegments)
  # xi  = XI0;
  eta = ETA0;
  xi = LinRange(-(XI0-dXI),XI0+channelLength, nSegments)
  x1 = @. 1/2 * (xi^2 - eta^2)
  # xi = LinRange(-(XI0-), XI0, nSegments)
  y1 = @. -xi*eta

  # resonator wall at 90 deg
  eta = LinRange(-(ETA0+channelLength),-ETA0, nSegments)
  xi  = -(XI0-dXI);

  x2 = @. 1/2 * (xi^2 - eta^2)
  y2 = @. xi*eta

  x1 = [x2;x1]
  y1 = [y2;y1]
  
  # x1 = [reverse(x1); x1[2:end]]
  # y1 = [reverse(y1);-y1[2:end]]
  # x1 = reverse(x1);
  # y1 = reverse(y1)
  x1d, y1d = divideResonatorCurve(collect(x1),collect(y1),nSegments)
  x1m, y1m = calcMidpoints(x1d,y1d)
  

  # # second curve
  eta = LinRange(-(ETA0-dETA), ETA0+channelLength, nSegments)
  xi  = XI0;
  x = @. 1/2 * (xi^2 - eta^2)
  y = @. xi*eta

  eta = (ETA0-dETA)
  xi = LinRange(-(XI0+channelLength), -XI0,nSegments)
  x2,y2 = par2cart(xi,eta)

  x = [x2; x]
  y = [y2; y]
  # x = [reverse(x); x]
  # y = [reverse(y);-y]

  x,y = divideResonatorCurve(collect(x),collect(y),nSegments)
  xm, ym = calcMidpoints(x,y)

  ds = diag( (@. hypot(x - x', y - y')),1);
  ds1 = diag( (@. hypot(x1d - x1d', y1d - y1d')),1);
  ds = [ds; ds1]

  x = [x; x1d]
  y = [y; y1d]

  
  # return (x,y), (xm,ym), (x1,y1), (x1m, y1m)
  xm = [xm; x1m]
  ym = [ym; y1m]


  rij_mid = @. hypot(xm - xm', ym - ym');
  return x,y,xm, ym,rij_mid,ds

end

"""
  function buildTopChannel(XI0, ETA0, dξ, cL, N)

Builds the top channel in an assymetric parabolic resonator. Used for the mesh
discretization in `BWM_geometry.createResonator`.

# Arguments
- `XI0::Float64`: parabolic coordinate 1 (billiard's body)
- `ETA0::Float64`: parabolic coordinate 2 (billiard's body)
- `dξ::Float64`: upper channel/arm opening
- `cL::Float64`: channel/arm length (both channels)
- `N::Int`: number of segment discretization 
"""
function buildTopChannel(XI0::Float64, ETA0::Float64, dξ::Float64, cL::Float64, N::Int)
  @assert (abs(dξ) < abs(XI0)) ? true : throw(AssertionError("dXi cannot be greater than ETA0 "))
   # resonator wall at 90 deg
  eta = LinRange(-(ETA0+cL),-ETA0, N)
  xi  = -(XI0-dξ);

  x2, y2 =par2cart(xi,eta)

  x_channel_top_1 = x2
  y_channel_top_1 = y2
  # rest of wall, not used for midpoints, only for mesh
  eta = LinRange(-ETA0,-(ETA0+cL), N)
  xi  = -XI0;
  x2, y2 =par2cart(xi,eta)
  x_channel_top_2 = x2
  y_channel_top_2 = y2
  return tuple.(x_channel_top_1, y_channel_top_1),tuple.(x_channel_top_2, y_channel_top_2)

end

"""
  function buildTopChannel(XI0, ETA0, dη, cL, N)

Builds the bottom channel in an assymetric parabolic resonator. Used for the mesh
discretization in `BWM_geometry.createResonator`.

# Arguments
- `XI0::Float64`: parabolic coordinate 1 (billiard's body)
- `ETA0::Float64`: parabolic coordinate 2 (billiard's body)
- `dη::Float64`: upper channel/arm opening
- `cL::Float64`: channel/arm length (both channels)
- `N::Int`: number of segment discretization 
"""
function buildBotChannel(XI0::Float64, ETA0::Float64, dη::Float64, cL::Float64, N::Int)
  # variable xi
  @assert (abs(dη) < abs(ETA0)) ? true : throw(AssertionError("dXi cannot be greater than ETA0 "))
   # resonator wall at 90 deg
  eta = (ETA0 - dη)
  xi  = LinRange(-(XI0+cL),-XI0, N);

  x2, y2 = par2cart(xi, eta)

  x_channel_bot_1 = x2
  y_channel_bot_1 = y2
  # rest of wall, not used for midpoints, only for mesh
  eta = -ETA0
  xi  = LinRange(XI0,(XI0+cL), N);
  
  x2, y2 = par2cart(xi, eta)

  x_channel_bot_2 = x2
  y_channel_bot_2 = y2
  return tuple.(x_channel_bot_2, y_channel_bot_2),tuple.(x_channel_bot_1, y_channel_bot_1)

end

"""
  par2cart(_xi, _eta)

Converts parabolic-cylindrical to cartesian coordinates.
# Arguments
- `_xi`: PC1
- `_eta`: PC2
"""
function par2cart(_xi, _eta)
  X = @. 1/2 * (_xi^2 - _eta^2)
  Y = @. _xi * _eta
  return X,Y
end

"""
  par2cart(_xi, _eta)

Converts cartesian to parabolic-cylindrical coordinates.
# Arguments
- `x::Float64`
- `y::Float64`
"""
function cart2par(x::Float64,y::Float64)
  # convert cartesian to parabolic coordinates
  # convention to have the positive x axis as the 
  # branch cut
  z = -x + 1im * y; # convert to imaginary unit
  z = sqrt(2*z);
  _xi = imag(z);
  _eta = real(z);
  return _xi, _eta
end 

"""
  par2cart(_xi, _eta)

Converts cartesian to parabolic-cylindrical coordinates (vector form).
# Arguments
- `x::Vector`
- `y::Vector`
"""
function cart2par(x::Vector,y::Vector)
  # convert cartesian to parabolic coordinates
  # convention to have the positive x axis as the 
  # branch cut
  z = @. (-x + 1im * y); # convert to imaginary unit
  z = sqrt.(2*z);
  _xi = imag(z);
  _eta = real(z);
  return _xi, _eta
end 


"""
  buildAsymmetricResonator(XI0, ETA0; dXI=0.5, dETA=0.5,channelLength=0.5, nSegments=100)

Constructs the wall for an assymetric parabolic resonator. Returns all variables
needed for the BWM (regarding geometry).

# Arguments
- `XI0::Float64`: Parabolic wall 1
- `ETA0::Float64`: Parabolic wall 2
- `dXI::Float64`: Upper channel opening
- `dETA::Float64`: Lower channel opening
- `channelLength::Float64`: channel length
- `nSegments::Int`: number of discretization points

# Returns
- `x`: x coord
- `y`: y coord
- `xm`: x midpoint
- `ym`: y midpoint
- `rij_mid`: distance matrix |r - r'|
- `ds`: arc lengths
"""
function buildAsymmetricResonator(XI0::Float64, ETA0::Float64; dXI::Float64=0.5, dETA::Float64=0.5,channelLength::Float64=0.5, nSegments::Int=100)
  xi = LinRange(0, XI0, nSegments)
  eta = ETA0
  x,y = par2cart(xi, eta)
  # x = @. 1/2 * (xi^2 - eta^2)
  # y = @. xi*eta

  dXi = 0.2
  @assert (dXI < XI0) ? true : throw(AssertionError("dEta cannot be greater than ETA0 "))

  # eta = LinRange(0, (ETA0-dEta), nSegments)
  # xi  = XI0;
  eta = ETA0;
  xi = LinRange(-(XI0-dXI),XI0+channelLength, nSegments)
  x1 = @. 1/2 * (xi^2 - eta^2)
  # xi = LinRange(-(XI0-), XI0, nSegments)
  y1 = @. -xi*eta

  # resonator wall at 90 deg
  eta = LinRange(-(ETA0+channelLength),-ETA0, nSegments)
  xi  = -(XI0-dXI);

  x2 = @. 1/2 * (xi^2 - eta^2)
  y2 = @. xi*eta

  
  x1 = [x2;x1]
  y1 = [y2;y1]

  
  # x1 = [reverse(x1); x1[2:end]]
  # y1 = [reverse(y1);-y1[2:end]]
  # x1 = reverse(x1);
  # y1 = reverse(y1)
  x1d, y1d = divideResonatorCurve(collect(x1),collect(y1),nSegments)
  x1m, y1m = calcMidpoints(x1d,y1d)
  

  # # second curve
  eta = LinRange(-(ETA0-dETA), ETA0+channelLength, nSegments)
  xi  = XI0;
  x = @. 1/2 * (xi^2 - eta^2)
  y = @. xi*eta

  eta = (ETA0-dETA)
  xi = LinRange(-(XI0+channelLength), -XI0,nSegments)
  x2,y2 = par2cart(xi,eta)

  x = [x2; x]
  y = [y2; y]
  # x = [reverse(x); x]
  # y = [reverse(y);-y]

  x,y = divideResonatorCurve(collect(x),collect(y),nSegments)
  xm, ym = calcMidpoints(x,y)

  ds = diag( (@. hypot(x - x', y - y')),1);
  ds1 = diag( (@. hypot(x1d - x1d', y1d - y1d')),1);
  ds = [ds; ds1]

  x = [x; x1d]
  y = [y; y1d]

  
  # return (x,y), (xm,ym), (x1,y1), (x1m, y1m)
  xm = [xm; x1m]
  ym = [ym; y1m]


  rij_mid = @. hypot(xm - xm', ym - ym');
  return x,y,xm, ym,rij_mid,ds

end

"""
  divideCurve(x, y, numSegments)

Divides a curve into N equal segments (couldn't figure out how to interpolate in
`julia`).

# Arguments
- `x::Vector`: x array
- `y::Vector`: y array
- `numSegments::Int`: number of discretization points

# Returns
- points (x,y)
"""
function divideCurve(x::Vector{Float64}, y::Vector{Float64}, numSegments::Int)
  # calculate arc length T (dt=each segment)
  dx = diag(x .- x',1) 
  dy = diag(y .- y',1)
  dt = @. sqrt(dx^2+dy^2)
  dt = [0.0; dt]

  cumLengths = cumsum(dt)
  totalLength = cumLengths[end]
  segmentLength = totalLength / numSegments
  
  #xDivided = zeros(numSegments+1)
  #yDivided = zeros(numSegments+1)
  xDivided = zeros(numSegments)
  yDivided = zeros(numSegments)
  xDivided[1] = x[1]
  yDivided[1] = y[1]
  
  #for segment in 1:numSegments-1
  for segment in 1:numSegments-1
      targetLength = segmentLength * segment
      
      # Find the index of the point corresponding to the target length
      index = searchsortedfirst(cumLengths, targetLength)  # Returns the index of the first value in a greater than or equal to x
      
      # Interpolate to find the (x, y) coordinates at the target length
      fraction = (targetLength - cumLengths[index - 1]) / (cumLengths[index] - cumLengths[index - 1])
      xSegment = x[index - 1] + fraction * (x[index] - x[index - 1])
      ySegment = y[index - 1] + fraction * (y[index] - y[index - 1])
      
      xDivided[segment+1] = xSegment
      yDivided[segment+1] = ySegment
  end
  # return 
  return [xDivided; xDivided[1]], [yDivided; yDivided[1]]
end

"""
  divideCurve(x, y, numSegments)

Divides a curve (for use in an open resonator geometry) into N equal segments. 
For use with open curves.

# Arguments
- `x::Vector`: x array
- `y::Vector`: y array
- `numSegments::Int`: number of discretization points

# Returns
- points (x,y)
"""
function divideResonatorCurve(x::Vector{Float64}, y::Vector{Float64}, numSegments::Int)
  # calculate arc length T (dt=each segment)
  dx = diag(x .- x',1) 
  dy = diag(y .- y',1)
  dt = @. sqrt(dx^2+dy^2)
  dt = [0.0; dt]

  cumLengths = cumsum(dt)
  totalLength = cumLengths[end]
  segmentLength = totalLength / numSegments
  
  #xDivided = zeros(numSegments+1)
  #yDivided = zeros(numSegments+1)
  xDivided = zeros(numSegments)
  yDivided = zeros(numSegments)
  xDivided[1] = x[1]
  yDivided[1] = y[1]
  
  #for segment in 1:numSegments-1
  for segment in 1:numSegments-1
      targetLength = segmentLength * segment
      
      # Find the index of the point corresponding to the target length
      index = searchsortedfirst(cumLengths, targetLength)  # Returns the index of the first value in a greater than or equal to x
      
      # Interpolate to find the (x, y) coordinates at the target length
      fraction = (targetLength - cumLengths[index - 1]) / (cumLengths[index] - cumLengths[index - 1])
      xSegment = x[index - 1] + fraction * (x[index] - x[index - 1])
      ySegment = y[index - 1] + fraction * (y[index] - y[index - 1])
      
      xDivided[segment+1] = xSegment
      yDivided[segment+1] = ySegment
  end
  # return 
  return [xDivided; x[end]], [yDivided; y[end]]
end

"""
  divideCurve(x, y, ds)

Divides a curve (for use in an open resonator geometry) into N equal segments of
length `ds`. Use this when you want to match the arc length of another object.
For use with open boundaries.

# Arguments
- `x::Vector`: x array
- `y::Vector`: y array
- `ds::Float64`: length of line segments

# Returns
- points (x,y)
"""
function divideResonatorCurve(x::Vector{Float64}, y::Vector{Float64}, ds::Float64)
  # calculate arc length T (dt=each segment)
  dx = diag(x .- x', 1) 
  dy = diag(y .- y', 1)
  dt = @. sqrt(dx^2 + dy^2)
  dt = [0.0; dt]

  cumLengths = cumsum(dt)
  totalLength = cumLengths[end]
  numSegments = ceil(Int,totalLength / ds)  # Calculate the number of segments based on arc length ds
  segmentLength = totalLength / numSegments
  println(numSegments)
  
  xDivided = zeros(numSegments)
  yDivided = zeros(numSegments)
  xDivided[1] = x[1]
  yDivided[1] = y[1]
  
  for segment in 1:numSegments-1
      targetLength = segmentLength * segment
      
      index = searchsortedfirst(cumLengths, targetLength)
      fraction = (targetLength - cumLengths[index - 1]) / (cumLengths[index] - cumLengths[index - 1])
      xSegment = x[index - 1] + fraction * (x[index] - x[index - 1])
      ySegment = y[index - 1] + fraction * (y[index] - y[index - 1])
      
      xDivided[segment+1] = xSegment
      yDivided[segment+1] = ySegment
  end
  
  return [xDivided; x[end]], [yDivided; y[end]]
end

"""
  calcMidpoints(x, y)

Computes the midpoints of a set of (x,y) points.

# Arguments
- `x::Vector`: x array
- `y::Vector`: y array

# Returns
- midpoints (x, y)
"""
function calcMidpoints(x::Vector,y::Vector)
  # calculates calcMidpoints of two sets of points
  xmid = diag((x .+ x')/2,1);
  ymid = diag((y .+ y')/2,1);

  return xmid, ymid
end

function singleParabola(sigma0::Float64, tau0::Float64, numSegments::Int)
  # return parabolic confocal billiard
  sigma = sigma0;
  tau = tau0 * (2/numSegments) * (-numSegments/2:numSegments/2); 

  x = @. (tau^2 - sigma^2)/2;
  y = sigma*tau;

return x,y  # flip false midpoint
end



function buildBaseResonator(XI0::Float64, ETA0::Float64; XIC=XI0/4, ETAC=ETA0*2,nSegments=52)
  # XI0  = 1.0
  # ETA0 = XI0
  
  # XIC  = XI0 / 4
  # ETAC = ETA0 * 2
  
  # dXI  = XI0/ 6
  # dETA = dXI
  # channelLength = 1.0
  # nSegments = 52
  # channelLength=ETAC-ETA0
  eta = LinRange(ETAC,ETA0, nSegments)
  xi  = XIC
  x1  = @. 1/2 * (xi^2 - eta^2)
  y1  = @. xi*eta
  
  eta = ETA0
  xi  = LinRange(XIC, XI0, nSegments)
  x2  = @. 1/2 * (xi^2 - eta^2)
  y2  = @. xi*eta
  
  xs = [x1; x2]; ys = [y1; y2];
  
  ds = ConfocalParabolic.calcArcLength(xs, ys)
  
  # ds = diag( (@. hypot(xs - xs', ys - ys')),1);
  xs = [xs; -reverse(xs)]; ys = [ys;reverse(ys)]
  xs, ys = ConfocalParabolic.divideResonatorCurve(xs, ys, nSegments-1)
  xm, ym = ConfocalParabolic.calcMidpoints(xs, ys)
  ds  = ConfocalParabolic.calcArcLength(xs, ys)
  ds  = [ds; reverse(ds)]
  xs = [xs; reverse(xs)];  ys = [ys;reverse(-ys)]
  xm = [xm; reverse(xm)]; ym = [ym; reverse(-ym)]
  # ds  = [ds; ConfocalParabolic.calcArcLength(reverse(xs), reverse(-ys))]
  
  # xb, yb = ConfocalParabolic.createConfocalBilliard(XI0, ETA0, 300)
  
  rij = @. hypot(xm - xm', ym - ym');
  
  return xs,ys,xm,ym,rij,ds
end

function abaqusToMeshes(filename::String)
  # converts abaqus file (exported from HOHQMesh)
  # into readable Mesh by Meshes.jl

  _abaqus_mesh = AbaqusReader.abaqus_read_mesh(filename)
  _nodes = _abaqus_mesh["nodes"]
  _elems = _abaqus_mesh["elements"]
  
  _points = [tuple(_nodes[i][1:2]...) for i in 1:_nodes.count] # get set of points
  _connec = Meshes.connect.([tuple(_elems[i]...) for i in 1:_elems.count]) # get faces ?  
  return _points, _connec
end


function createConfocalBilliard(XI0::Float64, ETA0::Float64, n_segments::Int)
x,y = confocalBilliard(XI0, ETA0,n_segments)
x,y = divideCurve(x,y,n_segments)
xm, ym = calcMidpoints(x,y);
rij_mid =  @. hypot(xm - xm', ym - ym');
ds = diag( (@. hypot(x-x',y-y')),1 )
return x,y,xm,ym,rij_mid,ds
end 

function createResonator(n_segments::Int, XI0::Float64, ETA0::Float64, upper_channel_opening::Float64,lower_channel_opening::Float64, channel_length::Float64, dens::SVector)
  # n_segments             = 1200
  # XI0                    = 3.0
  # ETA0                   = 2.0
  # upper_channel_opening  = 0.4
  # lower_channel_opening  = 0.4
  # channel_length         = 1.0

  # 
  x,y,xm,ym,rijMid,ds = ConfocalParabolic.buildAsymmetricResonator(XI0, ETA0;
  dXI=upper_channel_opening,
  dETA=lower_channel_opening,
  channelLength=channel_length, 
  nSegments=n_segments)

  channel_top_left, channel_top_right = ConfocalParabolic.buildTopChannel(XI0, ETA0, upper_channel_opening, channel_length, n_segments)
  channel_bot_left ,channel_bot_right = ConfocalParabolic.buildBotChannel(XI0, ETA0, lower_channel_opening, channel_length, n_segments)
  Nx = dens[1]; Ny = dens[2]; N = [Nx, Ny, 0]


  CHANNEL = :top
  channel_upper = newProject("CHANNEL_TOP","out")

  bndry_points_top_left  = [LinRange(0.0,1.0,n_segments) first.(channel_top_left)  last.(channel_top_left) zeros(n_segments)]
  bndry_points_top_right = [LinRange(0.0,1.0,n_segments) first.(channel_top_right) last.(channel_top_right) zeros(n_segments)]
  spline_top_left = newSplineCurve("left",  n_segments, bndry_points_top_left)
  spline_top_right = newSplineCurve("right",  n_segments, bndry_points_top_right)

  setPolynomialOrder!(channel_upper, 2)
  setMeshFileFormat!(channel_upper, "ABAQUS")

  top_line_top = newEndPointsLineCurve("top",[last(channel_top_right)...,0.0],[first(channel_top_left)..., 0.0])
  top_line_bot = newEndPointsLineCurve("bot",[last(channel_top_left)...,0.0],[first(channel_top_right)..., 0.0])

  top_line_left = newEndPointsLineCurve("left",[first(channel_top_left)...,0.0],[last(channel_top_left)..., 0.0])
  top_line_right = newEndPointsLineCurve("right",[first(channel_top_right)...,0.0],[last(channel_top_right)..., 0.0])
  addCurveToOuterBoundary!(channel_upper, spline_top_left)
  addCurveToOuterBoundary!(channel_upper, top_line_bot)
  addCurveToOuterBoundary!(channel_upper, top_line_top)
  addCurveToOuterBoundary!(channel_upper, spline_top_right)
  
  bounds = [channel_top_right[end][2], channel_top_left[1][1],channel_top_left[end][2], channel_top_right[1][1]]; 
  addBackgroundGrid!(channel_upper, bounds, N)

  generate_mesh(channel_upper)

  CHANNEL=:bot
  channel_bottom = newProject("CHANNEL_BOT","out")
  setPolynomialOrder!(channel_bottom, 4)
  setMeshFileFormat!(channel_bottom, "ABAQUS")

  bndry_points_bot_left  = [LinRange(0.0,1.0,n_segments) first.(channel_bot_left)  last.(channel_bot_left) zeros(n_segments)]
  bndry_points_bot_right = [LinRange(0.0,1.0,n_segments) first.(channel_bot_right) last.(channel_bot_right) zeros(n_segments)]
  spline_bot_left = newSplineCurve("bot_left",  n_segments, bndry_points_bot_left)
  spline_bot_right = newSplineCurve("bot_right",  n_segments, bndry_points_bot_right)


  bot_line_top = newEndPointsLineCurve("bot_top",[last(channel_bot_right)...,0.0],[first(channel_bot_left)..., 0.0])
  bot_line_bot = newEndPointsLineCurve("bot_bot",[last(channel_bot_left)...,0.0],[first(channel_bot_right)..., 0.0])

  bot_line_left = newEndPointsLineCurve("bot_left",[first(channel_bot_left)...,0.0],[last(channel_bot_left)..., 0.0])
  bot_line_right = newEndPointsLineCurve("bot_right",[first(channel_bot_right)...,0.0],[last(channel_bot_right)..., 0.0])
  addCurveToOuterBoundary!(channel_bottom, spline_bot_left)
  addCurveToOuterBoundary!(channel_bottom, spline_bot_right)
  addCurveToOuterBoundary!(channel_bottom, bot_line_bot)
  addCurveToOuterBoundary!(channel_bottom, bot_line_top)
  

  # boudns y,x,y,x
  channel_bot_left ,channel_bot_right
  bounds = [channel_bot_right[end][2], channel_bot_left[1][1],channel_bot_left[end][2], channel_bot_right[1][1]]; 
  # Nx = 100; Ny = 100; N = [Nx, Ny, 0]
  addBackgroundGrid!(channel_bottom, bounds, N)
  generate_mesh(channel_bottom)

  # read mesh
  abaqus_top_p, abaqus_top_c = ConfocalParabolic.abaqusToMeshes(joinpath(@__DIR__, "out/CHANNEL_TOP.inp"))
  meshes_top = Meshes.SimpleMesh(abaqus_top_p, abaqus_top_c)
  abaqus_bot_p, abaqus_bot_c = ConfocalParabolic.abaqusToMeshes(joinpath(@__DIR__, "out/CHANNEL_BOT.inp"))
  meshes_bot = Meshes.SimpleMesh(abaqus_bot_p, abaqus_bot_c)

  return x,y,xm,ym,rijMid,ds, meshes_top, meshes_bot
end

"""
  confocalBilliardMesh(XI0, ETA0, n_curve)

Calculates a set of (x,y) points inside a billiard defined within XI0 and ETA0.

# Arguments
- `XI0::Float64`: 
- `ETA0::Float64`: 
- `n_curve::Int`: number of discretization points
"""
function confocalBilliardMesh(XI0::Float64, ETA0::Float64, n_curve::Int, Nx::Int, Ny::Int)
  # XI0 = 3.0 
  # ETA0 = 2.0
  # n_curve = 20
  x1,y1 = singleParabola(XI0, ETA0, n_curve)
  x2,y2 = singleParabola(ETA0, XI0, n_curve)

  bndry_points_right = collect.(tuple.(collect(LinRange(0,1,n_curve+1)),-x1,y1, 0.0))
  bndry_points_right=mapreduce(permutedims, vcat, bndry_points_right)
  spline_right = newSplineCurve("right", n_curve+1, bndry_points_right)
  
  bndry_points_left = collect.(tuple.(collect(LinRange(1,0,n_curve+1)),x2,y2, 0.0))
  bndry_points_left=mapreduce(permutedims, vcat, bndry_points_left)
  bndry_points_left=reverse(bndry_points_left, dims=1)
  spline_left = newSplineCurve("left", n_curve+1, bndry_points_left)

  bndry_points = [bndry_points_right; bndry_points_left]

  confocal = newProject("billiard","out")
  setPolynomialOrder!(confocal, 2)
  setMeshFileFormat!(confocal, "ABAQUS")
  # points = [(0.0, 0.0),(1.0, 0.0), (1.0, 1.0), (0.5, 1.0)]
  # boudns y,x,y,x
  # minx = min(minimum(x1), minimum(x2))
  # maxx = max(maximum(x1),maximum(x2))
  # miny = min(minimum(y1), minimum(y2))
  # maxy = max(maximum(y1),maximum(y2))
  bounds = [1.5, -1.0, -1.5, 1.0]
  # Nx = 50; Ny = 70
  N = [Nx, Ny, 0]
  addCurveToOuterBoundary!(confocal, spline_right)
  addCurveToOuterBoundary!(confocal, spline_left)
  # HOHQMesh.addRefinementRegionPoints!(confocal, spline_left)
  rline = newRefinementLine("ref", "smooth", [0.0,0.0,0.0], [0.0, 1.0,0.0],0.05, 2.0)
  # addRefinementRegion!(confocal, rline)

  addBackgroundGrid!(confocal, bounds, N)

  corner_top = newRefinementCenter("corner_top", "smooth", [bndry_points_right[end,2], bndry_points_right[end,3], 0.0], 0.05, 2.0)
  corner_bot = newRefinementCenter("corner_bot", "smooth", [bndry_points_right[1,2], bndry_points_right[1,3], 0.0], 0.05, 2.0)
  addRefinementRegion!(confocal, corner_top)
  addRefinementRegion!(confocal, corner_bot)

  setSmoothingStatus!(confocal,"ON")
  # fig
  generate_mesh(confocal)
  # abaqus_top_p, abaqus_top_c = abaqusToMeshes(joinpath(@__DIR__, "out/billiard.inp"))
  abaqus_top_p, abaqus_top_c = abaqusToMeshes(joinpath(@__DIR__, "../out/billiard.inp"))
  meshes_top = Meshes.SimpleMesh(abaqus_top_p, abaqus_top_c)
  return meshes_top
end

"""
Calculates arc length for `a` and `b` vectors.
"""
function calcArcLength(a::Vector, b::Vector)
  return diag( calcDistances(a,b),1);
end

function buildBaseTop(XI0::Float64, ETA0::Float64; XIC=XI0/4, ETAC=ETA0*2,nSegments=52)
  # XI0  = 1.0
  # ETA0 = XI0
  
  # XIC  = XI0 / 4
  # ETAC = ETA0 * 2
  
  # dXI  = XI0/ 6
  # dETA = dXI
  # channelLength = 1.0
  # nSegments = 52
  # channelLength=ETAC-ETA0
  eta = LinRange(ETAC,ETA0, nSegments)
  xi  = XIC
  x1  = @. 1/2 * (xi^2 - eta^2)
  y1  = @. xi*eta
  
  eta = ETA0
  xi  = LinRange(XIC, XI0, nSegments)
  x2  = @. 1/2 * (xi^2 - eta^2)
  y2  = @. xi*eta
  
  xs = [x1; x2]; ys = [y1; y2];
  
  ds = ConfocalParabolic.calcArcLength(xs, ys)
  
  # ds = diag( (@. hypot(xs - xs', ys - ys')),1);
  xs = [xs; -reverse(xs)]; ys = [ys;reverse(ys)]
  xs, ys = ConfocalParabolic.divideResonatorCurve(xs, ys, nSegments-1)
  xm, ym = ConfocalParabolic.calcMidpoints(xs, ys)
  ds  = ConfocalParabolic.calcArcLength(xs, ys)
  # ds  = [ds; reverse(ds)]
  # xs = [xs; reverse(xs)];  ys = [ys;reverse(-ys)]
  # xm = [xm; reverse(xm)]; ym = [ym; reverse(-ym)]
  # ds  = [ds; ConfocalParabolic.calcArcLength(reverse(xs), reverse(-ys))]
  
  # xb, yb = ConfocalParabolic.createConfocalBilliard(XI0, ETA0, 300)
  
  rij = @. hypot(xm - xm', ym - ym');
  
  return xs,ys,xm,ym,rij,ds
end

function createSymmetricResonator(n_segments::Int, XI0::Float64, ETA0::Float64, XIC::Float64,ETAC::Float64, dens::SVector)
  # n_segments             = 1200
  # XI0                    = 3.0
  # ETA0                   = 2.0
  # upper_channel_opening  = 0.4
  # lower_channel_opening  = 0.4
  # channel_length         = 1.0

  # 
  x,y,xm, ym, rij,ds = buildBaseResonator(XI0, ETA0; XIC=XIC, ETAC=ETAC,nSegments=n_segments)
  x_top,y_top,_,_,_,_ = buildBaseTop(XI0, ETA0; XIC=XIC, ETAC=ETAC,nSegments=n_segments)
  x_right, y_right = [x_top[end]; x_top[end]], [y_top[end]; -y_top[end]]
  x_bot, y_bot = reverse(x_top), -reverse(y_top);
  x_left, y_left = [x_top[1]; x_top[1]], [-y_top[1]; y_top[1]]

  Nx, Ny = dens; 
  NMESH = [Nx, Ny, 0];

  symmetric_res = newProject("SYM","out")
  setPolynomialOrder!(symmetric_res, 2) 
  setMeshFileFormat!(symmetric_res, "ABAQUS")

  bndry_points_top  = [LinRange(0.0,1.0,n_segments) first.(x_top)  last.(y_top) zeros(n_segments)]
  bndry_points_bot  = [LinRange(0.0,1.0,n_segments) first.(x_bot)  last.(y_bot) zeros(n_segments)]

  spline_top = newSplineCurve("top",  n_segments, bndry_points_top)
  spline_bot = newSplineCurve("bot",  n_segments, bndry_points_bot)

  line_top = newEndPointsLineCurve("top", [x_top[1],y_top[1],0.0],[x_top[end],y_top[end], 0.0])
  line_bot = newEndPointsLineCurve("bot", [x_bot[1],y_bot[1],0.0],[x_bot[end],y_bot[end], 0.0])
  # setPolynomialOrder!(symmetric_res, 2)
  # setMeshFileFormat!(symmetric_res, "ABAQUS")

  line_left  = newEndPointsLineCurve("left", [x_left[1],y_left[1],0.0],[x_left[end],y_left[end], 0.0])
  line_right  = newEndPointsLineCurve("right", [x_right[1],y_right[1],0.0],[x_right[end],y_right[end], 0.0])

  # all_lines = [bndry_points_top[1,2:end],line_top.

  addCurveToOuterBoundary!(symmetric_res,line_top )
  addCurveToOuterBoundary!(symmetric_res,line_right )
  addCurveToOuterBoundary!(symmetric_res,line_bot )
  addCurveToOuterBoundary!(symmetric_res,line_left )
  # bounds = [maximum(y_top) * 1.5, x_top[1] * 1.5, x_top[end] * 1.5, minimum(y_bot) * 1.5]; 
  bounds = [5.5, -5.5, -5.5, 5.5]
  addBackgroundGrid!(symmetric_res, bounds, NMESH)

  # generate_mesh(symmetric_res)
  HOHQMesh.plotProject!(symmetric_res, MODEL+GRID)
  abaqus_p, abaqus_c = ConfocalParabolic.abaqusToMeshes(joinpath(@__DIR__, "out/SYM.inp"))
  sym_mesh= Meshes.SimpleMesh(abaqus_p, abaqus_c)

  return x,y,xm,ym,rij,ds, sym_mesh
end

# lattices

# abstract type BravaisLattice end;

# struct ObliqueGrid <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
#   b::Float64
#   ϕ::Float64
# end

# struct RectangularGrid <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
#   b::Float64
# end

# struct CenteredRectangular <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
#   b::Float64
# end

# struct SquareGrid <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
# end

# struct HexagonalGrid <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
# end

# struct HoneyLattice <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
# end

# struct TriangularGrid <: BravaisLattice
#   o::SVector{2, Float64}
#   a::Float64
#   b::Float64
#   s::Float64 # shift
# end

# """
#   Builds a grid for certain lattice types centered at g.o (origin)
# """
# function buildGrid(g::SquareGrid, nx::Int, ny::Int)
#   return vec([SVector(g.a * i - g.o[1], g.a * j - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
# end

# function buildGrid(g::RectangularGrid, nx::Int, ny::Int)
#   return vec([SVector(g.a * i - g.o[1], g.b * j - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
# end

# function buildGrid(g::HexagonalGrid, nx::Int, ny::Int)
#   return vec([SVector(g.a * i - g.o[1] + g.a/2 * (j%2), (g.a*j)* sqrt(3)/2 - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
# end

# function buildGrid(g::CenteredRectangular, nx::Int, ny::Int)
#   return vec([SVector(g.a * i - g.o[1] + g.a/2 * (j%2), (g.a*j)* g.b/2 - g.o[2]) for i in 0:nx-1, j in 0:ny-1])
# end

# function buildGrid(g::TriangularGrid, _n::Int)
#   _points = SVector{2, Float64}[]

#   for i in 0:_n
#     for j in 0:i
#       push!(_points, SVector(g.o[1] + i*g.a - g.s*j, g.o[2] + g.b*j)) #sqrt(3)/2*
#     end
#   end 
#   return _points
# end

# function buildGrid(g::HoneyLattice, _n::Int)
#   _centers = buildGrid(TriangularGrid(g.o, g.a, g.a*sqrt(3)/2, g.a/2), _n)
#   [push!(_centers, SVector(c[1]+g.a/2, c[2]+g.a/3)) for c in _centers[1:end-_n-1]]

#   return _centers
# end




function circularLatticeMesh(centers_3d::Vector{Vector{Float64}}, 
                             R::Float64, 
                             x::Vector{Float64}, 
                             y::Vector{Float64}, 
                             n_grid::Tuple{Int64,Int64},
                              )
  scattering_project = newProject("lattice", "out")
  setMeshFileFormat!(scattering_project, "ABAQUS")
  
  circ = newCircularArcCurve.("rod",   
                               centers_3d, 
                               R, 
                               0.0, 360.0) 
  
  
  
  xf = maximum(x) + 2R ; x0 = minimum(x) - 2R; 
  y0 = minimum(y) - 2R; yf = maximum(y) + 2R
  # n_grid = (100,90)
  bounds = [yf, x0, y0, xf]
  N_grid = [n_grid[1], n_grid[2], 0]
  addBackgroundGrid!(scattering_project, bounds, N_grid)
  [addCurveToInnerBoundary!(scattering_project, c, "circle$i") for (i,c) in enumerate(circ)]
  generate_mesh(scattering_project)
  
  plotProject!(scattering_project, 5)
  
  _points, _connec = BWM_geometry.abaqusToMeshes("./out/lattice.inp")
  _mesh = Meshes.SimpleMesh(_points, _connec)
  return _mesh
  end
# end

