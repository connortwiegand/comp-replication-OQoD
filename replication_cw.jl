using Pkg,  LinearAlgebra, Parameters, Roots#, Gadfly, BasisMatrices, Optim, Plots, DataFrames, FastGaussQuadrature

@with_kw mutable struct HHModel
    #Production Params
    g::Float64  = 0.0185 # growth rate -- rate of technical progress
    δ::Float64 = 0.075 # depreciation
    θ::Float64 = 0.3 # captial’s share of income
    #Utility Params
    μ::Float64 = 1 # risk aversion parameter -- 1.5 in Aiyagari and McGrattan, we will use log preferences
    η::Float64 = 0.328 # related to labor elasticity -- consumption's share of utility
    β::Float64 = 0.991 # discount factor (0.998 in transfored economy)
    #Economy Params
    χ::Float64 = 0.082 # transfers -- ratio of government tranfers to GDP 
    γ::Float64 = 0.217 # government purchases -- ratio of government purchases to GDP 
    k::Float64 = 2.5 # capital output ratio
    b::Float64 = 2/3 # debt:GDP ratio
    #Labor Productivity is log AR(1) -- log(e_t) ~ AR(1) 
    ρ_ϵ::Float64 = 0.6 # AR(1) coefficient
    σ_ϵ::Float64 = 0.3 # AR(1) s.d.
    Nϵ::Int64 = 7 # number of states? **Unsure about this**
    ϵ::Vector{Float64} = zeros(0)
    Π::Matrix{Float64} = zeros(0,0)
    
    
    #Asset Constraints
    a̲::Float64 = 0. # borrowing constraint
    a̅::Float64 = 600. # uppber bound on asset accumulation **Unsure about this, consider changing**
    Na::Int64 = 100 #number of grid points for splines

    #Prices **Unsure about these**
    r̄::Float64 = .01
    w̄::Float64 = 1.

    #Solution -- may need to change this
    k::Int = 2 #type of interpolation
    Vf::Vector{Interpoland} = Interpoland[]
    cf::Vector{Interpoland} = Interpoland[]
    lf::Vector{Interpoland} = Interpoland[] # Unsure
end



"""
    U(HH::HHModel,c,l) -- Updated

    Takes in consumption, leisure, and η, spits out utility value
"""
function U(HH::HHModel, c, l)
    η = HH.η
    return η * log(c) + (1 - η) * log(l)
end


"""
    setupgrids_shocks!(HH::HHModel, curv=1.7) -- Started updating

Set up non-linear grids for interpolation

Notes: 
    *Old utility maximization constraint was 
        𝑐+𝑎′=(1+𝑟¯)𝑎+𝑤¯𝜖_𝑠 
    New constraint is 
        𝑐+(1+g)𝑎′≤(1+r̄)𝑎 + w̄*e*(1-l)+χ
    I think ϵ is David's notation, e(t) is Aiyagari's notation
"""
function setupgrids_shocks!(HH::HHModel, curv=1.7)
    @unpack a̲,a̅,Na,ρ_ϵ,σ_ϵ,Nϵ,k,r̄,w̄,β = HH
    #Compute grid on A
    agrid = (a̅-a̲).*LinRange(0,1,Na).^curv .+ a̲

    #Store markov chain
    mc = rouwenhorst(Nϵ,ρ_ϵ,σ_ϵ)
    HH.Π = Π = mc.p
    HH.ϵ = exp.(mc.state_values)

    #First guess of interpolation functions
    abasis = Basis(SplineParams(agrid,0,k))
    a = nodes(abasis)[1]

    Vf = HH.Vf = Vector{Interpoland}(undef,Nϵ)
    cf = HH.cf = Vector{Interpoland}(undef,Nϵ)
    for s in 1:Nϵ
        c = @. r̄*a + w̄*HH.ϵ[s] #* ## To Do
        l =  ## To Do
        V = U(HH,c,l)./(1-β)

        Vf[s]= Interpoland(abasis,V)
        cf[s]= Interpoland(abasis,c)
        lf[s]= Interpoland(abasis,l)
    end
end;



"""
    iterate_endogenousgrid(HH,a′grid,cf′)

Iterates on Euler equation using endogenous grid method
"""
function iterate_endogenousgrid(HH,a′grid,cf′)
    @unpack γ,ϵ,β,Nϵ,Π,r̄,w̄,a̲= HH
    c′ = zeros(length(a′grid),Nϵ)
    for s in 1:Nϵ
        c′[:,s]= cf′[s](a′grid)
    end

    EERHS = β*(1+r̄)*(c′).^(-γ)*Π' #RHS of Euler Equation
    c = EERHS.^(-1/γ)

    #compute implies assets
    a = ((c .+ a′grid) .- w̄ .*ϵ')./(1+r̄)

    cf = Vector{Interpoland}(undef,Nϵ)
    for s in 1:Nϵ
        if a[1,s]> a̲
            c̲ = r̄*a̲ + w̄*ϵ[s]
            cf[s]= Interpoland(Basis(SplineParams([a̲; a[:,s]],0,1)),[c̲;c[:,s]])
        else
            cf[s]= Interpoland(Basis(SplineParams(a[:,s],0,1)),c[:,s])
        end
    end
    return cf
end;


"""
    solveHHproblem_eg!(HH)

Solves the HH problem using the endogeneous grid method
"""
function solveHHproblem_eg!(HH,verbose=false)
    a′grid = nodes(HH.cf[1].basis)[1]#Get nodes for interpolation
    
    cf′ = iterate_endogenousgrid(HH,a′grid,HH.cf)
    diff = 1.
    while diff  > 1e-8
        HH.cf = iterate_endogenousgrid(HH,a′grid,cf′)
        diff = maximum(norm(cf′[s](a′grid)-HH.cf[s](a′grid),Inf) for s in 1:HH.Nϵ) 
        if verbose
            println(diff)
        end
        cf′ = HH.cf
    end
end
solveHHproblem_eg!(HH)
setupgrids_shocks!(HH)
@time solveHHproblem_eg!(HH);
