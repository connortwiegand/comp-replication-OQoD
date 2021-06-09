# Aiyagari Models
using BasisMatrices,LinearAlgebra,Parameters,Optim,QuantEcon,DataFrames,Gadfly,SparseArrays,Arpack
using Roots

#=
 1. added leisure to utility function   U(HH::HHModel,c)
 2. budget constraint                   setupgrids_shocks!(HH::HHModel, curv=1.7)
 3. Value function with  growth         setupgrids_shocks!(HH::HHModel, curv=1.7)
 =#

# structure and utility
#------------------------------------------------------------------------------------------------

    @with_kw mutable struct HHModel
        #Preference Parameters
        γ::Float64 = 1. #Risk aversion
        η::Float64 = 0.328 # consumption elasticity?            Evan: June 9
        β::Float64 = 0.985 #Discount Rate

        #Prices
        r̄::Float64 = .01
        w̄::Float64 = 1.

        #Asset Grid Parameters
        a̲::Float64 = 0. #Borrowing Constraint
        a̅::Float64 = 400. #Upper Bound on assets
        Na::Int64 = 100.

        #Income Process
        ρ_ϵ::Float64 = 0.9923
        σ_ϵ::Float64 = 0.0983
        Nϵ::Int64 = 7
        ϵ::Vector{Float64} = zeros(0)
        Π::Matrix{Float64} = zeros(0,0)

        #Solution
        k::Int = 2 #type of interpolation
        Vf::Vector{Interpoland} = Interpoland[]
        cf::Vector{Interpoland} = Interpoland[]
        lf::Vector{Interpoland} = Interpoland[]

        #Extra
        EΦ′::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
        Φ::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)

        #additional parameters without a Home
        χ::Float64 = 0.9923
        g::Float64 = 0.9923
    end


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
    end

    """
        U(HH::HHModel,c)
    """
    function U(HH::HHModel,c,l)
        γ = HH.γ
        η = HH.η
        if γ == 1
            return log.(c.^(η) .* l.^(1-η))
        else
            return ((c.^(η) .* l.^(1-η)).^(1-γ))./(1-γ)
        end
    end
#------------------------------------------------------------------------------------------------



# shocks
#------------------------------------------------------------------------------------------------
    """
        setupgrids_shocks!(HH::HHModel, curv=1.7)
    Set up non-linear grids for interpolation
    """
    function setupgrids_shocks!(HH::HHModel, curv=1.7)
        @unpack a̲,a̅,Na,ρ_ϵ,σ_ϵ,Nϵ,k,r̄,w̄,β,χ,γ,η,g = HH
        
        #Compute grid on A
        # grid over l ? 
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
        lf = HH.cf = Vector{Interpoland}(undef,Nϵ)
        for s in 1:Nϵ
            #------------------------------------------------------------------------------------------------------------------------ 
            l = 0.5 # fixed l until we get that figured out. 
            c = @. r̄*a + (1-l)*w̄*HH.ϵ[s] + χ            # Evan June 6
                #= 
                Im a bit confused about this, what about saving for next period?
                This is a budget constraint right hand side
            =#
            if γ == 1
                V = U(HH,c,l)./(1-β)
            else
                V = U(HH,c,l)./(1-(β(1+g)^[η(1-γ)] ))
            end
            # V = U(HH,c,l)./(1-β)
                #= 
                infinite sum  ∑(β(1+g)^[η(1-μ)] )^t  replaces ∑ β^t
                setting μ to 1 gave us a log utility function. It also cancels growth out of our value funciton
            =#
            #------------------------------------------------------------------------------------------------------------------------ 

            Vf[s]= Interpoland(abasis,V)
            cf[s]= Interpoland(abasis,c)
            lf[s]= Interpoland(abasis,l)
        end

        #Expectations of 1st derivative of Basis functions
        HH.EΦ′ = kron(Π,BasisMatrix(abasis,Direct(),nodes(abasis)[1],[1]).vals[1])
        HH.Φ = kron(Matrix{Float64}(I,Nϵ,Nϵ),BasisMatrix(abasis,Direct()).vals[1])
    end


    """
        setupgrids_shocks!(AM::AiyagariModel)
    Setup the grids and shocks for the aiyagari model
    """
    function setupgrids_shocks!(AM::AiyagariModel,curv=2.)
        @unpack HH,Ia,N̄= AM
        @unpack a̲,a̅,Nϵ = HH
        setupgrids_shocks!(HH)
        #Normalize so that average labor supply is 1
        πstat = real(eigs(HH.Π',nev=1)[2])
        πstat ./= sum(πstat)
        HH.ϵ = HH.ϵ./dot(πstat,HH.ϵ)*N̄
        #Grid for distribution
        agrid = (a̅-a̲).*LinRange(0,1,Ia).^curv .+ a̲
        AM.z̄ = hcat(kron(ones(Nϵ),agrid),kron(1:Nϵ,ones(Ia)))
        AM.ω̄ = ones(Ia*Nϵ)/(Ia*Nϵ)
    end

#------------------------------------------------------------------------------------------------------------------------ 






#                                          Endogenous Grid
#------------------------------------------------------------------------------------------------

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
    end


    """
    solveHHproblem_eg!(HH)
    Solves the HH problem using the endogeneous grid method
    """
    function solveHHproblem_eg!(HH)
    a′grid = nodes(HH.cf[1].basis)[1]#Get nodes for interpolation

    cf′ = iterate_endogenousgrid(HH,a′grid,HH.cf)
    diff = 1.
    while diff  > 1e-8
        HH.cf = iterate_endogenousgrid(HH,a′grid,cf′)
        diff = maximum(norm(cf′[s](a′grid)-HH.cf[s](a′grid),Inf) for s in 1:HH.Nϵ) 
        cf′ = HH.cf
    end
    end
#------------------------------------------------------------------------------------------------

# junk
#------------------------------------------------------------------------------------------------
    # only implemented for HH control below not for Aiyagari
    """
    ee_errors(HH,agrid) 
    Check's the Euler equation errors of a solution on a grid of assets
    """
    function ee_errors(HH,agrid)
    @unpack cf,γ,β,ϵ,Π,Nϵ,r̄,w̄,a̲ = HH

    EE = zeros(length(agrid),Nϵ)
    for s in 1:Nϵ
        c = cf[s](agrid)
        a′ = @. (1+r̄)*agrid + w̄*ϵ[s]- c
        c′ = hcat([cf[s′](a′) for s′ in 1:Nϵ]...)
        ĉ  = (β.*(1+r̄).*c′.^(-γ)*Π[s,:]).^(-1/γ)
        EE[:,s]= log10.(abs.(ĉ.-c)./ĉ)
        EE[abs.(a′.-a̲).<1e-4,s].= -6 #record a small error when on borrowing constraints
    end

    df = DataFrame(EE,["ϵ_$s" for s in 1:Nϵ])
    df.a = agrid
    return df,maximum(EE)
    end


#------------------------------------------------------------------------------------------------



# Main course

# find stationary distribution and calibratesteadystate
#------------------------------------------------------------------------------------------------
    """
        find_stationarydistribution!(AM::AiyagariModel,V)
    Computes the stationary distribution 
    """
    function find_stationarydistribution!(AM::AiyagariModel)
        @unpack Ia,z̄,HH,W̄,R̄ = AM
        @unpack ϵ,Π,Nϵ,cf,a̲,a̅ = HH

        a = z̄[1:Ia,1] #grids are all the same for all shocks
        c = hcat([cf[s](a) for s in 1:Nϵ]...) #consumption policy IaxNϵ
        a′ = R̄.*a .+ W̄.*ϵ' .- c #create a IaxNϵ grid for the policy rules
        
        #make sure we don't go beyond bounds.  Shouldn't bind if bmax is correct
        a′ = max.(min.(a′,a̅),a̲)
        
        Qa = BasisMatrix(Basis(SplineParams(a,0,1)),Direct(),reshape(a′,Ia*Nϵ)).vals[1]
        Q = spzeros(Ia*Nϵ,Ia*Nϵ)
        for s in 1:Nϵ
            Q[1+(s-1)*Ia:s*Ia,:] = kron(reshape(Π[s,:],1,:),Qa[1+(s-1)*Ia:s*Ia,:]) 
        end
        
        AM.H = Q'
        AM.ω̄ .= real(eigs(AM.H;nev=1)[2])[:]
        AM.ω̄ ./= sum(AM.ω̄) 
    end

    function calibratesteadystate!(AM)
        @unpack Θ̄,α,N̄,K2Y,R̄ = AM
        AM.HH.r̄ = R̄ - 1
        Y2K = 1/K2Y
        AM.δ = α*Y2K + 1 - R̄ #matching capital to output ratio and interest rate gives depreciation rate
        K2N = (Y2K/Θ̄)^(1/(α-1)) #relationship between capital to output and capital to labor
        K̄ = K2N*N̄
        AM.W̄ = AM.HH.w̄ = (1-α)*Θ̄*K2N^α
    
        setupgrids_shocks!(AM)
        function βres(β)
            AM.HH.β=β
            # solvebellman!(AM.HH)      # collocation
            solveHHproblem_eg!(AM.HH)   # Endogenopus grid
            find_stationarydistribution!(AM)
            return dot(AM.ω̄,AM.z̄[:,1]) -K̄
        end
    
        Q̄ = 1/R̄
        fzero(βres,Q̄^2,Q̄^1.2)
    end
    
#------------------------------------------------------------------------------------------------

# interest determination
#------------------------------------------------------------------------------------------------    
    function capital_supply_demand(AM,R)
        @unpack Θ̄,α,N̄,δ = AM
        AM.R̄ = R
        AM.HH.r̄ = R-1

        Y2K = (R-1+δ)/α
        K2N = (Y2K/Θ̄)^(1/(α-1))
        AM.W̄ = AM.HH.w̄ = (1-α)*Θ̄*K2N^α
        KD = K2N * AM.N̄ 

        # solvebellman!(AM.HH)      # collocation
        solveHHproblem_eg!(AM.HH)   # Endogenopus grid
        find_stationarydistribution!(AM)
        KS = dot(AM.ω̄,AM.z̄[:,1])

        return [KS,KD]
    end



# 


# run functions
#------------------------------------------------------------------------------------------------  
    AM = AiyagariModel()
    AM.HH.β = 0.99
    setupgrids_shocks!(AM)




    # Figure 1: capital steady state: 
    Rgrid = LinRange(1.,1.007,10)
    KSKD = hcat([capital_supply_demand(AM,R) for R in Rgrid]...)
    plot(layer(y=Rgrid,x=KSKD[1,:],color=["Capital Supplied"],Geom.line),
        layer(y=Rgrid,x=KSKD[2,:],color=["Capital Demanded"],Geom.line),
        Guide.ylabel("Gross Interest Rate"), Guide.xlabel("Capital"))


    # 



    # this is a bit out of the blue

    AM.R̄ = 1.01 #target a quarterly interest rate of 1%


    # one last function

    AM = AiyagariModel()

    calibratesteadystate!(AM)

    plot(x=AM.z̄[AM.z̄[:,1].>1.,1],y=AM.ω̄[AM.z̄[:,1].>1.],Geom.bar)
#------------------------------------------------------------------------------------------------  