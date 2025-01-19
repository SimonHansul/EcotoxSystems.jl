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
        model = DEBODE!,
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

- `p`: Parameters, given as component vector. Typically at least two components are given: `glb` for global parameters, `spc` for species-level parameters.

kwargs

The following kwargs are used internally by `OrdinaryDiffEq.solve`. See `OrdinaryDiffEq` documentation for more information.

- `alg`: ODE solving algorithm to use.. 
- `saveat`: Interval or time-points at which to save the solution. Default is 1.
- `reltol`: Relative tolerance of ODE solution, default is 1e-6.
- `model`: ODE system to solve. 
- `callbacks`: A callback set, i.e. pairs of conditions and events. Used to handle discontinuities.  

In addition we have some kwargs that are used to further process inputs and outputs: 

- `statevars_init`: Function that defines initial state variables as component vector. Components typically match those in the parameter vector. 
- `gen_ind_params: Function that converts species-level parameters to individual-level parameters, for example by replacing parameters which are given as distributions with a random sample from the distribution. The inputs and outputs of this function should contain all components. 
- `returntype`: Indicating of how to return the result. Currently allowed are `dataframe` (complete solution converted to a `DataFrame`) and `odesol` (the ODE solution object as returned by `OrdinaryDiffEq`). Default is `dataframe`.

**Example**: 

```Julia
import EcotoxSystems as DEB
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
    gen_ind_params = generate_individual_params,
    param_links::Union{Nothing,NamedTuple} = nothing,  
    callbacks = DEBODE_callbacks,
    returntype::ReturnType = dataframe,
    kwargs...
    )

    p_ind = gen_ind_params(p) # converts spc component to ind component
    link_params!(p_ind, param_links) # apply parameter 

    u = statevars_init(p_ind)

    prob = ODEProblem(model, u, (0, p.glb.t_max), p_ind) # define the problem
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
    IBM_simulator(
        p::ComponentVector; 
        dt = 1/24, 
        saveat = 1,
        showinfo::Number = Inf
        )::DataFrame

Simulate the individual-based version of the default model. 

```
import EcotoxSystems as DEB
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
    init_global_statevars = initialize_global_statevars,
    global_ode! = DEBODE_global!,
    global_rules! = default_global_rules!,
    
    individual_ode! = DEBODE_individual!,
    individual_rules! = default_individual_rules!,
    init_individual_statevars = initialize_individual_statevars,
    gen_ind_params = generate_individual_params,
    
    dt = 1/24, 
    saveat = 1,
    record_individuals = true,
    showinfo::Number = Inf
    )

    showinfo < Inf ? @info("Running IndividualBasedModel simulation with t_max=$(p.glb.t_max)") : nothing
    
    global m = IndividualBasedModel(
        p; 

        init_global_statevars = init_global_statevars,
        global_ode! = global_ode!, 
        global_rules! = global_rules!,
        
        individual_ode! = individual_ode!,
        individual_rules! = individual_rules!,
        init_individual_statevars = init_individual_statevars,
        gen_ind_params = gen_ind_params,
        
        dt = dt, 
        saveat = saveat,
        record_individuals = record_individuals
        )

    while !(m.t > m.p.glb.t_max)
        if showinfo < Inf && isapprox(m.t % showinfo, 0, atol = m.dt)
            @info("t=$(m.t), N = $(m.u.glb.N), mem = $(Base.summarysize(m))")
        end

        model_step!(m)
    end


    return (glb = global_record_to_df(m), spc = individual_record_to_df(m))
end


global_record_to_df(m::AbstractDEBIBM)::DataFrame = DataFrame(hcat(m.global_record...)', getcolnames(m))

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
        sim = []

        for replicate in 1:$nreps
            sim_i = $(esc(simcall))
            push!(sim, sim_i)
        end
        combine_outputs(Vector{typeof(sim[1])}(sim))
    end
end


"""
    replicates(simulator::Function, defaultparams::ComponentVector, nreps::Int64; kwargs...)

Perform replicated runs of `simulator` with parameters `defaultparams` (`simulator(defaultparams)` has to be a valid function call). 
Analogous to `@replicates`, but a bit more flexible.
"""
function replicates(simulator::Function, defaultparams::ComponentVector, nreps::Int64; kwargs...)
    sim = []

    for replicate in 1:nreps
        sim_i = simulator(defaultparams; kwargs...)
        push!(sim, sim_i)
    end
    
    return combine_outputs(Vector{typeof(sim[1])}(sim))
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

"""
    exposure(
        simulator::Function, 
        p::ComponentVector, 
        C_Wmat::Matrix{R}
    ) where R <: Real

Simulate exposure to an arbitrary number of stressors over a Matrix of constant exposure concentrations `C_Wmat`. 

The exposure matrix columns are stressors, the rows are treatments. 
That means, to simulate a single-stressor experiment, do 

```Julia 
C_Wmat = [0.; 1.; 2;]
```

, creating a n x 1 matrix with exposure concentrations 0, 1 and 2. 

In contrast, a single treatment with multiple stressors would be defined as 

```Julia
C_Wmat = [0 1 2;]
```

, creating a 1 x n matrix. Here, we would have three stressors with the simultaneous exposure concentrations 0,1,2. <br>

Thus, defining four treatments for two stressors looks like this:

```Julia
C_Wmat = [0.0 0.0; 0.0 0.5; 1.0 0.0; 0.5 1.0]
```

This exposure matrix corresponds to a ray design with constant exposure ratios.
"""
function exposure(
    simulator::Function, 
    p::ComponentVector, 
    C_Wmat::Matrix{R}
    ) where R <: Real
    
    let C_W_int = p.glb.C_W # we will modify this value and then reset to the initial value
        sim = []

        for (i,C_W) in enumerate(eachrow(C_Wmat))
            p.glb.C_W = C_W
            sim_i = simulator(p)
            typeof(C_W) <: Number ? sim_i = add_idcol(sim_i, :C_W, C_W) : nothing
            sim_i = add_idcol(sim_i, :treatment_id, i)
            push!(sim, sim_i)
        end
        
        p.glb.C_W = C_W_int 

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