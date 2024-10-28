function extract_colnames(c::R, k::Symbol) where R <: Real
    return k
end

function extract_colnames(c::AbstractVector, k::Symbol)
    return [Symbol("$(k)_$(i)") for i in 1:length(c)]
end

function extract_colnames(u::ComponentVector)
    colnames = []
    for k in keys(u)
        push!(colnames, extract_colnames(u[k], k))
    end
    return vcat(colnames...)
end

extract_colnames(sol::ODESolution)::Vector{Symbol} = vcat([:t], extract_colnames(sol.u[1]))

function sol_to_mat(sol::O)::Matrix{Float64} where O <: ODESolution
    return hcat(sol.t, hcat(sol.u...)')
end

function sol_to_df(sol::O)::DataFrame where O <: ODESolution
    simout = DataFrame(
        sol_to_mat(sol), # ODE output converted to matrix
        extract_colnames(sol) # column names inferred from component array keys
    )
    return simout
end

robustmin(x) = length(x)>0 ? minimum(x) : Inf

ismin(x::AbstractVector) = x .== minimum(x)


"""
    trim!(spc::ComponentVector)

Trim TKTD parameters to the smallest number of stressors indicated for any PMoA.
For example `k_D_G = [1., 0.], k_D_M = [1.]` will be trimmed to `k_D_G = [1.], k_D_M = [1.]`.

- `spc`: An instance of a species defaultparams type
"""
function trim!(spc::ComponentVector)

    num_stressors = min(
        length(spc.k_D_G),
        length(spc.k_D_M),
        length(spc.k_D_A),
        length(spc.k_D_R),
        length(spc.k_D_h)
    )

    spc.k_D_G = spc.k_D_G[1:num_stressors]
    spc.k_D_M = spc.k_D_M[1:num_stressors]
    spc.k_D_A = spc.k_D_A[1:num_stressors]
    spc.k_D_R = spc.k_D_R[1:num_stressors]
    spc.k_D_h = spc.k_D_h[1:num_stressors]

    spc.drc_functs_G = spc.drc_functs_G[1:num_stressors]
    spc.drc_functs_M = spc.drc_functs_M[1:num_stressors]
    spc.drc_functs_A = spc.drc_functs_A[1:num_stressors]
    spc.drc_functs_R = spc.drc_functs_R[1:num_stressors]
    spc.drc_functs_h = spc.drc_functs_h[1:num_stressors]
    
    spc.drc_params_G = spc.drc_params_G[1:num_stressors]
    spc.drc_params_M = spc.drc_params_M[1:num_stressors]
    spc.drc_params_A = spc.drc_params_A[1:num_stressors]
    spc.drc_params_R = spc.drc_params_R[1:num_stressors]
    spc.drc_params_h = spc.drc_params_h[1:num_stressors]
end



function rr(x::R, x_ref::R)::Float64 where R <: Real
    if x_ref == 0.
        return 1.
    else
        return x / x_ref
    end
end

function rr(x::R, x_ref::Missing)::Missing where R <: Real
    return missing
end

robustmin(x) = length(x)>0 ? minimum(x) : Inf

function robustmean(x)
    xfilt = filter(xi -> isfinite(xi) & !ismissing(xi), x)

    if length(xfilt)==0
        return NaN
    end

    return mean(xfilt)
end

"""
    relative_response(
        sim::D, 
        response_vars::Vector{Symbol},
        treatment_var::Symbol; 
        groupby_vars::Vector{Symbol} = Symbol[],
        identify_control = minimum
        ) where D <: AbstractDataFrame

Calculate relative responses. \n

args:

- `sim::AbstractDataFrame`: results
- `response_vars::Vector{Symbol}`: response variables for which to calculate the relative responses
- `treatment_var::Symbol`: Column indicating the treatment. Column values can be numerical or categorical, but `identify_control` kwarg has to be specified in the latter case

kwargs:

- `groupby_vars::Vector{Symbol}`: relative response will be conditioned on these variables (e.g. time, separate experiments...). Empty by default.
- `identify_control`: function used to identify reference values from `treatment_var`. By default, this is `minimum()` (assuming that column values in `treatment_var` are numerical).

---

## Example 

```Julia@f
using MechanisticEffectModels.EnergyBudgetDiffEqs, MechanisticEffectModels.Utils

simfunct(x) = @replicates EnergyBudgetDiffEqs.simulator(x) 10

sim = exposure(simfunct, Params(), [0., 100., 200.]) |>
x -> relative_response(x, [:S, :H, :R]) # -> data frame will contain columns y_S, y_H, y_R for control-normalized values
```

"""
function relative_response(
    sim::D, 
    response_vars::Vector{Symbol},
    treatment_var::Symbol; 
    groupby_vars::Vector{Symbol} = Symbol[],
    identify_control = minimum
    ) where D <: AbstractDataFrame

    #=
    Calculation of the conditional control mean
    =#
    reference = sim[sim[:,treatment_var] .== identify_control(sim[:,treatment_var]),:] |> # extract the control
    x -> groupby(x, groupby_vars) |> # group by conditioning variables
    x -> combine(x) do df 
        refvals = DataFrame() # (sub-)dataframe of reference values
        for var in response_vars # iterate over response values
            var_ref = Symbol(String(var) * "_ref") # get the reference column name
            refvals[!,var_ref] = [robustmean(df[:,var])] # calculate the conditional control mean
        end
        return refvals
    end

    #=
    Calculation of the relative response
    =#
    sim = leftjoin(sim, reference, on = groupby_vars) # add references as new columns
    for var in response_vars # for every response variable
        y_var = Symbol("y_" * String(var)) # get the relative response column name
        var_ref = Symbol(String(var) * "_ref") # get the reference column name
        sim[!,y_var] = [rr(row[var], row[var_ref]) for row in eachrow(sim)] # calculate the relative response
        select!(sim, Not(var_ref)) # drop the reference column
    end
    return sim
end


function idcol!(sim::AbstractDataFrame, col::Symbol, val)
    sim[!,col] .= val
end


function idcol!(sim::Any, col::Symbol, val)
    for df in sim
        idcol!(df, col, val)
    end
end


