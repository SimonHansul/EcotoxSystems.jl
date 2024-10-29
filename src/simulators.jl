

@enum ReturnType odesol dataframe

"""
    simulator(
        p::ComponentVector; 
        alg = Tsit5(),
        saveat = 1,
        reltol = 1e-6,
        model = DEBODE!,
        callbacks = DEBODE_callbacks,
        returntype::ReturnType = dataframe,
        kwargs...
    )::DataFrame

Run the model as ODE system. 


**Example**: 

```Julia
import EnergyBudgetDiffEqs as DEB
p = DEB.params() # loads the default parameters
sim = DEB.ODE_simulator(p)
```

"""
function ODE_simulator(
    p::ComponentVector; 
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    model = DEBODE!,
    statevars_init = initialize_statevars,
    ind_params_init = generate_individual_params,
    callbacks = DEBODE_callbacks,
    returntype::ReturnType = dataframe,
    kwargs...
    )

    individual_params = ind_params_init(p)
    u = statevars_init(individual_params)

    prob = ODEProblem(model, u, (0, p.glb.t_max), individual_params) # define the problem
    sol = solve(prob, alg; callback = callbacks, saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP

    if returntype == odesol # directly return the ODE solution object
        return sol
    end

    if returntype == dataframe # return solution as dataframe
        return  DataFrame(hcat(sol.t, hcat(sol.u...)'), getcolnames(sol)) # convert solution to dataframe
    end

    error("returntype $returntype not implemented")
end


"""
    @replicates(simcall::Expr, nreps::Int64) 

Perform replicated runs of `simcall`, where `simcall` is a call to a simulator function. 

Example:

```Julia
    spc = SpeciesParams(Z = Truncated(Normal(1, 0.1), 0, Inf)) # initialize default parameters with variable zoom factor
    sim = @replicates MechanisticEffectModels.simulator(Params(spc = spc))) 10 # execute replicated runs to simulator
```

In this case, `sim` will contain the output of 10 replicated simulations. For each replicate, the zoom factor is sampled from a truncated Normal distribution. 
`sim` contains an additional column `replicate`.
"""
macro replicates(simcall::Expr, nreps::Int)
    quote
        sim = DataFrame()

        for replicate in 1:$nreps
            sim_i = $(esc(simcall))
            sim_i[!,:replicate] .= replicate
            append!(sim, sim_i)
        end
        sim
    end
end

"""
    IBM_simulator(
        p::ComponentVector; 
        dt = 1/24, 
        saveat = 1,
        showinfo::Number = Inf
        )::DataFrame

Simulate the individual-based version of the default model. 

```
import EnergyBudgetDiffEqs as DEB
p = DEB.params()
sim = DEB.IBMsimulator(p)
```

args

- `p`: The parameter collection with defined global and species parameters.

kwargs

- `dt`: Length of a timestep in the model (unit according to chosen unit of rate parameters)
- `saveat`: Time interval at which to record output
- `showinfo`: Time interval at which to print an update. Nothing will be printed if `showinfo == Inf` (the default).

"""
function IBM_simulator(
    p::ComponentVector; 
    global_ode! = DEBODE_global!,
    individual_ode! = DEBODE_individual!,
    individual_rules! = default_individual_rules!,
    dt = 1/24, 
    saveat = 1,
    record_individuals = true,
    showinfo::Number = Inf
    )::DataFrame

    showinfo < Inf ? @info("Running IndividualBasedModel simulation with t_max=$(p.glb.t_max)") : nothing
    
    m = IndividualBasedModel(
        p; 
        global_ode! = global_ode!, 
        individual_ode! = individual_ode!,
        individual_rules! = individual_rules!,
        dt = dt, 
        saveat = saveat,
        record_individuals = record_individuals
        )

    while !(m.t > m.p.glb.t_max)
        if showinfo < Inf && isapprox(m.t % showinfo, 0, atol = m.dt)
            @info("t=$(m.t)") : nothing
        end

        model_step!(m)
    end

    return individual_record_to_df(m)
end

"""
individual_record_to_df(
    m::AbstractDEBIBM; 
    cols::Union{Symbol,Vector{Symbol}} = :all
    )::DataFrame

Convert individual record to Data Frame.

args 

- `m::AbstractDEBIBM`: Model object containing an `individual_record` field as Vector of Component Vectors

kwargs

- `cols`: A Vector of Symbols indicating the 

"""
function individual_record_to_df(
    m::AbstractDEBIBM; 
    )::DataFrame


    cols = vcat(
        [:t, :id],
        [keys(m.individual_record[1])...]
    )  |> unique

    
    hcat([map(x -> getproperty(x, y), m.individual_record) for y in cols]...) |> 
    x -> DataFrame(x, cols)
end



"""
    replicates(simulator::Function, defaultparams::ComponentVector, nreps::Int64; kwargs...)

Perform replicated runs of `simulator` with parameters `defaultparams` (`simulator(defaultparams)` has to be a valid function call). 
Analogous to `@replicates`, but a bit more flexible.
"""
function replicates(simulator::Function, defaultparams::ComponentVector, nreps::Int64; kwargs...)
    sim = DataFrame()

    for replicate in 1:nreps
        sim_i = simulator(defaultparams; kwargs...)
        sim_i[!,:replicate] .= replicate
        append!(sim, sim_i)
    end
    
    sim
end

"""
    treplicates(
        simulator::Function, 
        defaultparams::ComponentVector, 
        nreps::Int64; 
        kwargs...)

Multi-threaded version of `replicates`. 

Only useful if Julia has been started with multiple threads. 
    
To check the number of threads, run 
```using Base.Threads; Threads.nthreads()```.

In VSCode, you can use the entry "julia.NumThreads" in settings.json to set the default number of threads 
(searching for "julia threads" in the preferences will lead you there). 

Check the [Multi-threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/) 
for more information.

"""
function treplicates(
    simulator::Function, 
    defaultparams::ComponentVector, 
    nreps::Int64; 
    kwargs...)

    sim = Vector{DataFrame}(undef, nreps)

    @threads for replicate in 1:nreps
        sim_i = simulator(defaultparams; kwargs...)
        sim_i[!,:replicate] .= replicate
        sim[replicate] = sim_i
    end
    
    return vcat(sim...)
end

set_vecval(C_W::Real) = [C_W]
set_vecval(C_W::AbstractVector) = C_W


"""
    exposure(simcall::Expr, C_Wvec::Vector{Float64}; kwargs...)

Simulate exposure to a single stressor over a Vector of constant exposure concentrations `C_Wvec`. 
"""
function exposure(
    simulator::Function, 
    p::ComponentVector, 
    C_Wmat::Matrix{R}
    ) where R <: Real
    
    let C_W_int = p.glb.C_W # we will modify this value and then reset to the initial value
        sim = DataFrame()

        for (i,C_W) in enumerate(eachrow(C_Wmat))
            p.glb.C_W = C_W
            sim_i = simulator(p)
            sim_i[!,:treatment_id] .= i
            append!(sim, sim_i)
        end
        
        p.glb.C_W = C_W_int 

        return sim
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


