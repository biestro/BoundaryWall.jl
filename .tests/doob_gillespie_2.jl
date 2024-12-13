using GLMakie
using StatsBase
using StaticArrays
using Distributions


"""
  Transition Operator for SEIR

Only when N is a function of the current population

# Parameters
- `_STATE::Int64`: current state 
"""
function seirTransitionOperator(_STATE::SVector{4, Int64}, _MATRIX::SMatrix{4,4,Float64})
  # S*I*_beta/N, E*_sigma, I*_gamma
  _N = sum(_STATE)
  return (_STATE[1]*_STATE[3]*_MATRIX[2,3]/_N, _STATE[2]*_MATRIX[1,2], _STATE[3]*_MATRIX[3,4])

end

"""
  Gillespie

Runs the Doob Gillespie algorithm
# Parameters
- `_MAX_ITER::Int64`: max number of iterations
- `_INIT::SVector{N, Int64}`: initial state
- `_PARAMS::SMatrix{N,N, Float64}`: stochastic matrix

"""
function Gillespie(_MAX_ITER::Int64, _INIT::SVector, _PARAMS::SMatrix)

  # initialize vectors
  _TIME = [0.0]
  _STATE = [SVector(_INIT...)]
  
  _N = sum(_INIT)
   
  

  for _ in 1:_MAX_ITER
    # S*I*_beta/N, E*_sigma, I*_gamma
    current_state = last(_STATE)
    propensities = seirTransitionOperator(current_state, _PARAMS)
    # println(propensities)
    sum_propensity = cumsum(propensities)
  
    if sum_propensity[end] <= 0
      println("finished...")
      break
    end
  
     # geometric time jump (time to next reaction)
    delta_time = rand(Distributions.Exponential(1/sum_propensity[end]))
    push!(_TIME, last(_TIME)+delta_time)
  
    # uniform sample for 
    rand_num = rand(Distributions.Uniform(0,1))
  
    if (rand_num * sum_propensity[end]) < sum_propensity[1]
      push!(_STATE, current_state .+ SVector(-1,1,0,0))  # susceptible becomes exposed
    # elseif (rand_num * propensity) < sum(propensities[1:2])
    elseif (rand_num * sum_propensity[end]) < sum_propensity[2]
      push!(_STATE, current_state .+ SVector(0,-1,1,0))  # exposed becomes infected
    elseif (rand_num * sum_propensity[end]) < sum_propensity[3]
      push!(_STATE, current_state .+ SVector(0,0,-1,1))  # infected recovers
    # elseif (rand_num * propensity) < sum(propensities[1:4])
    #   push!(_STATE, current_state .+ SVector(0,0,0,-1))  # lost immunity
    # elseif (rand_num * propensity) < propensities[4]
    #   push!(_STATE, current_state .+ SVector(1,0,0,-1))
    else
      push!(_STATE, current_state)
    end
  end
  return _TIME,_STATE
  # return time, state
end


begin
S₀ = 1_000
E₀ = 10
I₀ = 10
R₀ = 0

σ, β, γ, μ = 0.01, 0.3, 0.05, 0.001

n_states = 50

init              = SA[S₀,E₀,I₀,R₀]
transition_matrix = SA[0.0 σ   0.0 0.0;
                       0.0 0.0 β   0.0;
                       0.0 0.0 0.0 γ  ;
                       0.0 0.0 0.0 0.0 ]


RUNS = [Gillespie(4000,init, transition_matrix) for _ in 1:n_states]
time = first.(RUNS)
state = last.(RUNS)

using Interpolations

n_interp = 100
interp_data = zeros(SVector{4, Float64},n_interp)
minimum_time, min_id = findmin(map(t -> t[end], time))
interp_time = LinRange(0.0, minimum_time, n_interp)

total_people = sum(init)

# number_of_people = zeros(length(time), length(init))
susceptible_hist = [[], [], [], []]

for i in eachindex(time)
itp = linear_interpolation(time[i], state[i])
interp_data .+= itp.(interp_time)

# number_of_people[i,:] = itp(half_time)
push!(susceptible_hist[1], getindex.(itp.(interp_time), 1))
push!(susceptible_hist[2], getindex.(itp.(interp_time), 2))
push!(susceptible_hist[3], getindex.(itp.(interp_time), 3))
push!(susceptible_hist[4], getindex.(itp.(interp_time), 4))
end
interp_data ./= length(time) # mean

fig = Figure()
ax = Axis(fig[1,1])
for (t,s) in zip(time, state)
lines!(ax,t, getindex.(s,1), color=(Makie.wong_colors()[1], 0.25), linewidth=0.5)
lines!(ax,t, getindex.(s,2), color=(Makie.wong_colors()[2], 0.25), linewidth=0.5)
lines!(ax,t, getindex.(s,3), color=(Makie.wong_colors()[3], 0.25), linewidth=0.5)
lines!(ax,t, getindex.(s,4), color=(Makie.wong_colors()[4], 0.25), linewidth=0.5)
end

g = ["S", "E", "I", "R"]
for i in eachindex(interp_data[1])
lines!(ax, interp_time, getindex.(interp_data, i), linewidth=3, label=g[i])
end
axislegend(ax, position=:rc)

total_interp = sum(reduce(vcat, transpose.(interp_data)), dims=2)

lines!(ax, interp_time,    vec(total_interp) )
xlims!(ax, interp_time[1], interp_time[end]  )

# ax.aspect=DataAspect()
fig

end

lines(getindex.(interp_data,2), getindex.(interp_data,3))
