abstract type AbstractModel end
abstract type AbstractEnergyBudget <: AbstractModel end
abstract type SimplifiedEnergyBudget <: AbstractEnergyBudget end

# TODO: add callback_set to keep this consistent with EnergyBudgetModelZoo
Base.@kwdef mutable struct FullDEBkiss <: AbstractEnergyBudget
    parameters::ComponentVector = debkiss_defaultparams

    # global component

    #initialize_global_statevars::Union{Function,Nothing} = debkiss_global_statevars
    #global_derivatives!::Union{Function,Nothing} = constant_nutrient_influx!
    #global_rules!::Function = default_global_rules!
    
    # individual-level component
    
    #initialize_individual_statevars::Function = debkiss_individual_statevars
    #individual_derivatives!::Function = debkiss!
    #individual_rules!::Function = default_individual_rules!
    #generate_individual_params::Function = debkiss_individual_params

    # composed model

    #initialize_all_statevars::Union{Function,Nothing} = nothing
    #complete_derivatives!::Union{Function,Nothing} = nothing
    #callback_set::CallbackSet = debkiss_callback_set
end

function instantiate(deb::FullDEBkiss; verbose = false)::FullDEBkiss
    
    #compose_derivatives!(deb)
    #compose_statevars!(deb)
    #compound_parameters!(deb; verbose = verbose)
    
    return deb
end


function compose_statevars!(deb::FullDEBkiss)::Nothing
    let init_global_statevars = isnothing(deb.initialize_global_statevars) ? p -> ComponentVector() : deb.initialize_global_statevars, 
    init_individual_statevars = deb.initialize_individual_statevars

        deb.initialize_all_statevars = function (p::ComponentVector)
            ComponentVector(
                glb = init_global_statevars(p),
                ind = init_individual_statevars(p)
            )
        end

        return nothing
    end
end

simulate(debkiss::FullDEBkiss; kwargs...) = sim_all(debkiss.parameters; kwargs...)


function simulate_ibm(
    deb::FullDEBkiss; 
    dt = 1/24, 
    saveat = 1,
    record_individuals = true,
    showinfo = Inf,
    kwargs...
    )

     sim = EcotoxSystems.IBM_simulator(
        deb.parameters;
        init_global_statevars = deb.initialize_global_statevars,
        global_ode! = deb.global_derivatives!,
        global_rules! = deb.global_rules!,

        individual_ode! = deb.individual_derivatives!,
        individual_rules! = deb.individual_rules!,
        init_individual_statevars = deb.initialize_individual_statevars,
        gen_ind_params = deb.generate_individual_params,

        dt = dt, 
        saveat = saveat, 
        record_individuals = record_individuals,
        showinfo = showinfo,
        kwargs...
        )

    return sim
end

