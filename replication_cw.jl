using Pkg,  LinearAlgebra, Parameters#, Roots, Gadfly, BasisMatrices, Optim, Plots, DataFrames, FastGaussQuadrature

@with_kw mutable struct RParams
    g::Float64  = 0.0185 #
    χ::Float64 = 0.082 #
    δ::Float64 = 0.075 #
    ρ::Float64 = 0.6 # 
    μ::Float64 = 1.5 #
    γ::Float64 = 0.217 #
    θ::Float64 = 0.3 #
    η::Float64 = 0.328 #
    σ::Float64 = 0.3 #
    β::Float64 = 0.991 #
end

