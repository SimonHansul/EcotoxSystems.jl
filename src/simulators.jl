

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

**args**:

- `p::ComponentVector`: Parameters defined as ComponentVector. Default models assume entries `glb` and `spc` for global and species-specific parameters, respectively. `spc` will be internally converted to `ind` to represent individual-specific parameters. 

**kwargs**:

- `model = EnergyBudgetDiffEqs!`: Definition of the derivatives. A function form `du!(du, u, p, t)::Nothing`. See definition of `EnergyBudgetDiffEqs!` in `derivatives.jl` for an example, or see docs and examples of [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) for more details.
- `alg = Tsit5()`: Algorithm to be used by `solve` function 
- `saveat = 1`: When or how often to save ODE solutions to output
- `reltol = 1e-6`: Relative tolerance of ODE solutions
- `model` : Definition of the derivatives in a form that is compatible with the DifferentialEquations.jl interface
- `callbacks` : Definition of callbacks, e.g. for life-stage transitions
- `returntype` : The return type as enum `ReturnType`. 
- `kwargs...`: Additional keyword argument are passed on to `OrdinaryDiffEq.solve`

**Example**: 

```Julia
sim = simulator(Params())
```

"""
function simulator(
    p::ComponentVector; 
    alg = Tsit5(),
    saveat = 1,
    reltol = 1e-6,
    model = DEBODE!,
    callbacks = DEBODE_callbacks,
    returntype::ReturnType = dataframe,
    kwargs...
    )

    individual_params = generate_individual_params(p)
    u = initialize_statevars(individual_params)

    prob = ODEProblem(model, u, (0, p.glb.t_max), individual_params) # define the problem
    sol = solve(prob, alg; callback = callbacks, saveat = saveat, reltol = reltol, kwargs...) # get solution to the IVP

    if returntype == odesol
        return sol
    end

    if returntype == dataframe
        return  DataFrame(hcat(sol.t, hcat(sol.u...)'), getcolnames(sol)) # convert solution to dataframe
    end

    error("returntype $returntype does not seem to be implemented")
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


