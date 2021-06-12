using Pkg,BasisMatrices,LinearAlgebra,Parameters,Roots,Optim,QuantEcon,DataFrames,Gadfly,SparseArrays,Arpack

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
    i_type::Int = 2 #type of interpolation -- e.g. 2 = quadratic, **look out for this, used to be called k**
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
    return η*log.(c) .+ (1-η)*log.(l)
end


"""
    setupgrids_shocks!(HH::HHModel, curv=1.7) -- Good to go I think

Set up non-linear grids for interpolation

Notes: 
    *Old utility maximization constraint was 
        𝑐+𝑎′=(1+𝑟¯)𝑎+𝑤¯𝜖_𝑠 
    New constraint is 
        𝑐+(1+g)𝑎′≤(1+r̄)𝑎 + w̄*e*(1-l)+χ
    I think ϵ is David's notation, e(t) is Aiyagari's notation
"""
function setupgrids_shocks!(HH::HHModel, curv=1.7)
    @unpack a̲,a̅,Na,ρ_ϵ,σ_ϵ,Nϵ,i_type,r̄,w̄,β,η,γ,g,b = HH
    #Compute grid on A
    agrid = (a̅-a̲).*LinRange(0,1,Na).^curv .+ a̲

    τ = γ + (r̄-g)b 

    #Store markov chain
    mc = rouwenhorst(Nϵ,ρ_ϵ,σ_ϵ)
    HH.Π = Π = mc.p
    HH.ϵ = exp.(mc.state_values)

    #First guess of interpolation functions
    abasis = Basis(SplineParams(agrid,0,i_type))
    a = nodes(abasis)[1]

    Vf = HH.Vf = Vector{Interpoland}(undef,Nϵ)
    cf = HH.cf = Vector{Interpoland}(undef,Nϵ)
    #lf = HH.lf = Vector{Interpoland}(undef,Nϵ)
    for s in 1:Nϵ
        # NOTE: These are guesses! So they can be off
            # - We will see later if they need to be adjusted, but don't get too caught up on them
        #c = @. r̄*a + w̄*HH.ϵ[s] + χ
        #l = @. (η/(1-η)) * w̄ * HH.ϵ[s] 
        c = @. max(1/(η)*( (r̄-g)*a + w̄*HH.ϵ[s] + (-τ)), 0.001)# I interpreted these as the budget contraient where a = a' and we plug in static FOC
        l = @. max(((1-η)*c)/(w̄*HH.ϵ[s]* η), 0.001)# static FOC
        V = U(HH,c,l)./(1-β)

        Vf[s]= Interpoland(abasis,V)
        cf[s]= Interpoland(abasis,c)
        #It doesn't currently run with lf[s]. Unsure if this needed anyay, but maybe try setting i_type = 3
        ###lf[s]= Interpoland(abasis,l) 
    end
end;


HH = HHModel()
setupgrids_shocks!(HH)

"""
    iterate_endogenousgrid(HH,a′grid,cf′) -- Started updating

Iterates on Euler equation using endogenous grid method
"""
function iterate_endogenousgrid(HH,a′grid,cf′)
    @unpack γ,ϵ,β,Nϵ,Π,r̄,w̄,a̲,g,η,χ,b = HH
    c′ = zeros(length(a′grid),Nϵ)
    for s in 1:Nϵ
        c′[:,s]= cf′[s](a′grid)
    end

    EERHS = β*(1+r̄)/(1+g)*η*(c′).^(-1)*Π' #RHS of Euler Equation
    c = η * EERHS.^(-1)

    #Compute leisure from consumption (using FOCs)
        #NOTE: I have no idea what needs dots and what doesn't
    l = ((1-η)/η) * (c*(w̄*ϵ).^(-1))

    # Compute taxes
    τ = γ + (r̄-g)b 

    #compute implied assets from BC
    a = ((c .+ w̄.*ϵ'.*l .+ ((1+g) * a′grid)) .- w̄ .* ϵ' .+ τ) ./ (1+r̄) # note: that is ϵ', not ϵ′ -- dumb

    cf = Vector{Interpoland}(undef,Nϵ)
    for s in 1:Nϵ
        if a[1,s]> a̲
            c̲ = η * ( (r̄-g)*a̲ + w̄*ϵ[s] - τ )
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

@with_kw mutable struct AiyagariModel
    HH::HHModel = HHModel()

    #Production Parameters
    α::Float64 = 0.3
    δ::Float64 = 0.025
    Θ̄::Float64 = 1.

    #Moments to match/prices
    W̄::Float64 = 1.
    R̄::Float64 = 1.01
    K2Y::Float64 = 10.2 #capital to output ratio
    N̄::Float64 = 1.

    #Distribution Parameters
    Ia::Int = 1000 #Number of gridpoints for distribution
    z̄::Matrix{Float64} = zeros(0,0) #Gridpoints for the state variables
    ω̄::Vector{Float64} = zeros(0) #Fraction of agents at each grid level
    H::SparseMatrixCSC{Float64,Int64} = spzeros(Ia,Ia) #Transition matrix
end;