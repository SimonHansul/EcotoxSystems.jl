abstract type AbstractEnergyBudget end


Base.@kwdef mutable struct SimplifiedEnergyBudget <: AbstractEnergyBudget
    parameters::ComponentVector = default_debkiss_params
    initialize_global_statevars::Union{Function,Nothing} = debkiss_global_statevars
    initialize_individual_statevars::Function = debkiss_individual_statevars
    initialize_all_statevars::Union{Function,othing} = nothing
    global_derivatives!::Union{Function,Nothing} = constant_nutrient_influx!
    individual_derivatives!::Function = debkiss!
    complete_derivatives!::Union{Function,Nothing} = nothing
    generate_individual_params::Function = generate_individual_params
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

function compose_derivatives!(std::StandardEnergyBudget)::Nothing
    if !isnothing(std.global_derivatives!)
        function std_complete_ode!(du, u, p, t)::Nothing
            std.global_derivatives!(du, u, p, t)
            std.individual_derivatives!(du, u, p, t)
            return nothing
        end
        std.complete_derivatives! = std_complete_ode!
        return nothing
    else
        std.complete_derivatives! = std.individual_derivatives!
        return nothing
    end
end

function simulate(deb::AbstractEnergyBudget; kwargs...)
    sim = EcotoxSystems.ODE_simulator(
            deb.parameters;
            model = deb.complete_derivatives!, 
            statevars_init = deb.initialize_all_statevars,
            gen_ind_params = deb.generate_individual_params, 
            kwargs...
        )

    return sim
end
