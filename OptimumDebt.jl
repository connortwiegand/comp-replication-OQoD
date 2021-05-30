
using Parameters # mutable structure 
# parameterization for our benchmark
@with_kw mutable struct BenchmarkParameters
    g::Float64 = 0.0185 # growth rate
    γ::Float64 = 0.217  # ratio of government purchases to GDP
    χ::Float64 = 0.082  # ratio of government tranfers to GDP
    b::Float64 = 2/3    # average of the US federal plus US state debt divided by gross domestic product over the post-war period
    θ::Float64 = 0.3   # labor’s share of income equal to 70%
    k::Float64 = 2.5    # capital output ratio
    δ::Float64 = 0.075 #Depreciation Rate
    ρ::Float64 = 0.6 #  parameter values for the earnings process are in the range of estimates in the literature
    σ::Float64 = 0.3# (See Aiyagari (1994a) for a discussion of these estimates
    η::Float64 = 0.328
    μ::Float64 = 1.5
    β::Float64 = 0.991 #Discount Factor
end;



# graphs to replicate on page 16/23
#=
Welfare Gain
Interest Rate
Tax Rate
Aggregate Hours
=#




