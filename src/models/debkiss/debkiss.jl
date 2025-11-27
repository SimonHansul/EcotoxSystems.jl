abstract type AbstractEnergyBudget end

Base.@kwdef mutable struct SimplifiedEnergyBudget <: AbstractEnergyBudget
    parameters::ComponentVector = default_debkiss_params

    # global component

    initialize_global_statevars::Union{Function,Nothing} = debkiss_global_statevars
    global_derivatives!::Union{Function,Nothing} = constant_nutrient_influx!
    global_rules!::Function = default_global_rules!
    
    # individual-level component
    
    initialize_individual_statevars::Function = debkiss_individual_statevars
    individual_derivatives!::Function = debkiss!
    individual_rules!::Function = default_individual_rules!
    generate_individual_params::Function = generate_individual_params

    # composed model

    initialize_all_statevars::Union{Function,Nothing} = nothing
    complete_derivatives!::Union{Function,Nothing} = nothing
end

function instantiate(deb::SimplifiedEnergyBudget; verbose = false)::SimplifiedEnergyBudget
    
    compose_derivatives!(deb)
    compose_statevars!(deb)
    #compound_parameters!(deb; verbose = verbose)
    
    return deb
end

function compose_statevars!(deb::AbstractEnergyBudget)::Nothing
    let init_global_statevars

        if isnothing(deb.initialize_global_statevars)
            init_global_statevars = p -> ComponentVector()
        else
            init_global_statevars = deb.initialize_global_statevars
        end

        function initialize_statevars(p::ComponentVector)::ComponentVector
            return ComponentVector(
                glb = init_global_statevars(p), 
                ind = deb.initialize_individual_statevars(p)
            )
        end

        deb.initialize_all_statevars = initialize_statevars

        return nothing
    end
end

function compose_derivatives!(deb::AbstractEnergyBudget)::Nothing
    if !isnothing(deb.global_derivatives!)
        function std_complete_ode!(du, u, p, t)::Nothing
            deb.global_derivatives!(du, u, p, t)
            deb.individual_derivatives!(du, u, p, t)
            return nothing
        end
        deb.complete_derivatives! = std_complete_ode!
        return nothing
    else
        deb.complete_derivatives! = deb.individual_derivatives!
        return nothing
    end
end

function simulate_ode(deb::AbstractEnergyBudget; kwargs...)
    
    sim = EcotoxSystems.ODE_simulator(
            deb.parameters;
            model = deb.complete_derivatives!, 
            statevars_init = deb.initialize_all_statevars,
            gen_ind_params = deb.generate_individual_params, 
            kwargs...
        )

    return sim

end

function simulate_ibm(
    deb::AbstractEnergyBudget; 
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

const simulate_ODE = simulate_ode
const simulate_IBM = simulate_ibm
