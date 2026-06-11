abstract type AbstractModel end
abstract type AbstractEnergyBudget <: AbstractModel end
abstract type SimplifiedEnergyBudget <: AbstractEnergyBudget end

# TODO: add callback_set to keep this consistent with EnergyBudgetModelZoo
Base.@kwdef mutable struct FullDEBkiss <: AbstractEnergyBudget
    parameters::ComponentVector = defaultparams()
end

simulate(debkiss::FullDEBkiss; kwargs...) = sim_all(debkiss.parameters; kwargs...)

#function simulate_ibm(
#    deb::FullDEBkiss; 
#    dt = 1/24, 
#    saveat = 1,
#    record_individuals = true,
#    showinfo = Inf,
#    kwargs...
#    )
#
#    sim = EcotoxSystems.IBM_simulator(
#        deb.parameters;
#        init_global_statevars = deb.initialize_global_statevars,
#        global_ode! = deb.global_derivatives!,
#        global_rules! = deb.global_rules!,
#
#        individual_ode! = deb.individual_derivatives!,
#        individual_rules! = deb.individual_rules!,
#        init_individual_statevars = deb.initialize_individual_statevars,
#        gen_ind_params = deb.generate_individual_params,
#
#        dt = dt, 
#        saveat = saveat, 
#        record_individuals = record_individuals,
#        showinfo = showinfo,
#        kwargs...
#    )
#
#    return sim
#end


"""
Simulate constant exposure to a single stressor for DEBkiss model.
"""
function simulate_constant_exposure(
    p::ComponentVector, 
    C::Pair{Symbol,Vector{R}}; 
    kwargs...
    ) where R <: Real

    sim = DataFrame()

    C_label = C.first 
    C_values = C.second

    for (i,C_value) in enumerate(C_values) # iterate over concentrations
        setindex!(p.glb, C_value, C_label)
        sim_i = sim_all(p; kwargs...)
        sim_i[!,:treatment_id] .= i
        append!(sim, sim_i)
    end

    return sim
end


"""

```Julia
DEBkiss.simulate_constant_exposure(
    FullDEBkiss().parameters,
    (:C_W1, :C_W2) => [
        0. 0. # control : C_W1 and C_W2 are 0
        0. 1. # exposure conc = 1 for C_W2
        1. 0. # exposure conc = 1 for C_W1
        1. 1. # exposure conc = 1 for both
    ] 
)
```
"""
function simulate_constant_exposure(
    p::ComponentVector, 
    C::Pair;
    kwargs...
    )

    sim = DataFrame()

    for (i,C_values) in enumerate(eachrow(C.second))
        for (C_label,C_value) in zip(C.first,C_values)
            setindex!(p.glb, C_value, C_label)
        end
        
        sim_i = sim_all(p; kwargs...)
        sim_i[!,:treatment_id] .= i
        append!(sim, sim_i)
    end

    return sim
end