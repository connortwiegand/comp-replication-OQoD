using Pkg,  LinearAlgebra, Parameters, Roots#, Gadfly, BasisMatrices, Optim, Plots, DataFrames, FastGaussQuadrature

@with_kw mutable struct RParams
    g::Float64  = 0.0185 # growth rate
    χ::Float64 = 0.082 # ratio of government tranfers to GDP
    δ::Float64 = 0.075 # depreciation
    ρ::Float64 = 0.6 # 
    μ::Float64 = 1.5 #
    γ::Float64 = 0.217 # ratio of government purchases to GDP
    θ::Float64 = 0.3 # captial’s share of income
    η::Float64 = 0.328 #
    σ::Float64 = 0.3 #
    β::Float64 = 0.991 # discount factor
    k::Float64 = 2.5 # capital output ratio
    b::Float64 = 2/3 # debt:GDP ratio
end

"""
Section 3.3, from the technical appendix
'#' Indicates regular comments
'##' indicates to-dos 

N, r_l, and r_u are guessed
θ, δ, and b are parameters
tol is a tolerance level

Need methods for computing α and H to do this part
"""

function try_r(N, r_guess, para::RParams) 
    @unpack g,χ,γ,b,θ = para
    ## Determine τ_y from government budget constraint
    τ_y = (γ + χ + r_guess*b - g*b)/(1 + r_guess*b - δ*θ/(r+δ))

    #"at" stands for after tax
    ir_at = (1 - τ_y) * r_guess
    wr_at = (1 - θ) / N

    ## Use these two items to compute the finite element approximation for α

    ## Also compute a finite element approximation for H 
        ##-- TA p. 198, compute (13) with finite element method with linear basis functions & residual from (9)
        ##-- Possibly from Davids notes
    ## Use this to calculate E[ã_t] -- TA p. 199, in between (15) and (16)
    return E_at
end

function iterate_r(N, r_l, r_u, θ, δ, b, tol = 1e-6)
    r_guess = 0.5 * (r_l + r_u)
    E_at = try_r(N, r_guess, θ)

    while norm(Ea_t - ((θ / (r_guess + δ)) + b)) > tol
        if Ea_t < (θ / (r_guess + δ)) + b
            r_l = r_guess 
        elseif Ea_t > (θ / (r_guess + δ)) + b
            r_u = r_guess 
        end
        r_guess = 0.5 * (r_l + r_u)
        E_at = try_r(N, r_guess, θ)
    end

    return r_guess
end

## Now we need to check N
    ## Use Newton-Raphson, according to TA
