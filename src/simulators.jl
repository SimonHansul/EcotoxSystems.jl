# simulators.jl
# these functions are central to easily run simulations from parameters
# central functions are ODE_simulator and IBM_simulator 

@enum ReturnType odesol dataframe

"""
    ODE_simulator(
        p::ComponentVector; 
        alg = Tsit5(),
        saveat = 1,
        reltol = 1e-6,
        model = default_ODE!,
        statevars_init = initialize_statevars,
        gen_ind_params = generate_individual_params,
        param_links::Union{Nothing,NamedTuple} = nothing,  
        callbacks = DEBODE_callbacks,
        returntype::ReturnType = dataframe,
        kwargs...
    )

Run the model as ODE system. 
This function is essentially a wrapper around `OrdinaryDiffEq.solve` with additional input and output processing.

args:

- `p`: Parameters, given as component vector. 
At least two components are given: `glb` for global parameters, 
`spc` for species-level parameters. `glb` as to contain an entry `t_max` for the maximum simulation time. 
`spc` has to contain entries `Z` and `propagate_zoom` for the zoom factor and the parameters which are affected by the zoom factor.

kwargs:

The following kwargs are used internally by `OrdinaryDiffEq.solve`. See `OrdinaryDiffEq` documentation for more information.

- `alg`: ODE solving algorithm to use.. 
- `saveat`: Interval or time-points at which to save the solution. Default is 1.
- `reltol`: Relative tolerance of ODE solution, default is 1e-6.
- `model`: ODE system to solve. 
- `callbacks`: A callback set, i.e. pairs of conditions and events. Used to handle discontinuities.  

In addition we have some kwargs that are used to further process inputs and outputs: 

- `statevars_init`: Function with signature `statevars_init(p)` that defines initial state variables as component vector. Components typically match those in the parameter vector. 
- `gen_ind_params: Function that converts species-level parameters to individual-level parameters, for example by replacing parameters which are given as distributions with a random sample from the distribution. The inputs and outputs of this function should contain all components. 
- `returntype`: Indicating of how to return the result. Currently allowed are `dataframe` (complete solution converted to a `DataFrame`) and `odesol` (the ODE solution object as returned by `OrdinaryDiffEq`). Default is `dataframe`.

Run a model as purely ODE-based system: 

```Julia
using EcotoxSystems
sim = DEB.ODE_simulator(EcotoxSystems.defaultparams(p))
```

"""
function ODE_simulator(
    p::ComponentVector; 
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    model = default_ODE!,
    statevars_init = initialize_statevars,
    gen_ind_params = generate_individual_params,
    param_links::Union{Nothing,NamedTuple} = nothing,  
    returntype::ReturnType = dataframe,
    kwargs...
    )

    p_ind = gen_ind_params(p) # converts spc component to ind component
    link_params!(p_ind, param_links) # apply parameter 

    u = statevars_init(p_ind)

    prob = ODEProblem(model, u, (0, p.glb.t_max), p_ind) # define the problem
    sol = solve(prob, alg; saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP

    if returntype == odesol # directly return the ODE solution object
        return sol
    end

    if returntype == dataframe # return solution as dataframe
        return  DataFrame(hcat(sol.t, hcat(sol.u...)'), getcolnames(sol)) # convert solution to dataframe
    end

    error("returntype $returntype not implemented")
end


global_record_to_df(m::AbstractIBM)::DataFrame = DataFrame(hcat(m.global_record...)', getcolnames(m))

function individual_record_to_df(
    m::AbstractIBM; 
    )::DataFrame

    cols = vcat(
        [:t, :id],
        [keys(m.individual_record[1])...]
    )  |> unique

    
    hcat([map(x -> getproperty(x, y), m.individual_record) for y in cols]...) |> 
    x -> DataFrame(x, cols)
end

combine_outputs(outputs::Vector{DataFrame}; idcol = :replicate) = begin

    out = DataFrame()

    for (i,df) in enumerate(outputs)
        df[!,idcol] .= i 
        append!(out, df)
    end

    return out
end

combine_outputs(outputs::Vector{N}; idcol = :replicate) where N <: NamedTuple = begin

    out = [combine_outputs([out[i] for out in outputs], idcol = idcol) for i in 1:length(outputs[1])]

    return NamedTuple(zip(keys(outputs[1]), out))
end

macro replicates(simcall::Expr, nreps::Int)
    quote
        sim = []

        for replicate in 1:$nreps
            sim_i = $(esc(simcall))
            push!(sim, sim_i)
        end
        combine_outputs(Vector{typeof(sim[1])}(sim))
    end
end

function replicates(simulator::Function, defaultparams::ComponentVector, nreps::Int64; kwargs...)
    sim = []

    for replicate in 1:nreps
        sim_i = simulator(defaultparams; kwargs...)
        push!(sim, sim_i)
    end
    
    return combine_outputs(Vector{typeof(sim[1])}(sim))
end

function treplicates(
    simulator::Function, 
    defaultparams::ComponentVector, 
    nreps::Int64; 
    kwargs...)

    sim = Vector{Any}(undef, nreps)

    @threads for replicate in 1:nreps
        sim_i = simulator(defaultparams; kwargs...)
    end
    
    return combine_outputs(Vector{typeof(sim[1])}(sim))
end

set_vecval(C_W::Real) = [C_W]
set_vecval(C_W::AbstractVector) = C_W

function add_idcol(df::DataFrame, col::Symbol, val::Any)::DataFrame
    df[!,col] .= val 
    return df
end

function add_idcol(nt::NT, col::Symbol, val::Any) where NT <: NamedTuple

    for df in nt
        df[!,col] .= val
    end

    return nt
end

reshape_Cmat(Cmat::Matrix{Float64}) = Cmat
reshape_Cmat(Cmat::Vector{Float64}) = hcat(Cmat...)' |> Matrix

function exposure(
    simulator::Function, 
    p::ComponentVector, 
    Cmat::VecOrMat{R}
    ) where R <: Real
    
    Cmat = reshape_Cmat(Cmat) # if a Vector has been provided, we have to reshape

    let original_C_W = deepcopy(p.glb.C_W) # we will modify this value and then reset to the original value
        sim = [] # initialize vector of simulation outputs

        for (i,C_W) in enumerate(eachrow(Cmat)) # iterate over concentrations
            p.glb.C_W = C_W # update concentration
            sim_i = simulator(p) # simulate
            typeof(C_W) <: Number ? sim_i = add_idcol(sim_i, :C_W, C_W) : nothing
            sim_i = add_idcol(sim_i, :treatment_id, i)
            push!(sim, sim_i)
        end
        
        p.glb.C_W = original_C_W # reset C_W, so that the input remains unmodified

        # combine vector of simulation outputs into a single output object
        return combine_outputs(Vector{typeof(sim[1])}(sim); idcol = :treatment_id)
    end
end

function exposure(model, Cmat::VecOrMat{R}; C = :C_W, simfunc = simulate_ODE) where R <: Real
  
    Cmat = reshape_Cmat(Cmat) # if a Vector has been provided, we have to reshape
    p = model.parameters # reference to params 
    
    let original_C = deepcopy(p.glb[C]) # we will modify this value and then reset to the original value
        sim = [] # initialize vector of simulation outputs

        for (i,C_i) in enumerate(eachrow(Cmat)) # iterate over concentrations
            p.glb[C] = C_i # update concentration
            sim_i = simfunc(model) # simulate
            typeof(C_i) <: Number ? sim_i = add_idcol(sim_i, C, C_i) : nothing
            sim_i = add_idcol(sim_i, :treatment_id, i)
            push!(sim, sim_i)
        end
        
        p.glb[C] = original_C # reset C_W, so that the input remains unmodified

        # combine vector of simulation outputs into a single output object
        return combine_outputs(Vector{typeof(sim[1])}(sim); idcol = :treatment_id)
    end

end


"""
Threaded version of `exposure()`.
"""
function texposure(
    simulator::Function, 
    defaultparams::ComponentVector, 
    C_Wvec::Vector{Float64}
    )
    
    let C_W_int = defaultparams.glb.C_W # we will modify this value and then reset to the initial value
        sim = DataFrame()

        @threads for C_W in C_Wvec
            defaultparams.glb.C_W[1] = C_W
            sim_i = simulator(defaultparams)
            append!(sim, sim_i)
        end
        
        defaultparams.glb.C_W = C_W_int 

        return sim
    end
end