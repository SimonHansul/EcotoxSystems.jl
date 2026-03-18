# individuals.jl 
# definition of the default individual and default individual rules
# the mutabe struct for individuals should be generic enough for most applications, 
# the default_individual_rules should be viewed as a starting point for speicifc applications

abstract type AbstractIndividual end

CAUSE_OF_DEATH = Dict(
    0 => "none",
    1 => "age"
)

@inline function determine_S_max_hist(
    S::Float64,
    S_max_hist::Float64
    )::Float64

    return max(S, S_max_hist)

end

@inline function death_by_loss_of_structure(
    S::Float64,
    S_max_hist::Float64,
    W_S_rel_crit::Float64,
    h_S::Float64,
    dt::Float64
    )::Bool

    return ((S/S_max_hist) < W_S_rel_crit) && (rand() > exp(-h_S * dt))

end

@inline function death_by_aging(
    age::Float64,
    a_max::Float64
    )::Bool

    return age >= a_max

end

@inline function check_reproduction_period(time_since_last_repro::Float64, tau_R::Float64)::Bool
    return time_since_last_repro >= tau_R 
end

@inline function calc_num_offspring(R::Float64, X_emb_int::Float64)::Int64
    return trunc(R / X_emb_int)
end


@with_kw mutable struct Individual <: AbstractIndividual

    du::CVOrParamStruct
    u::CVOrParamStruct
    p::CVOrParamStruct
    
    individual_ode!::Function # equation-based portion of the individual step
    individual_rules!::Function # rule-based portion of the individual step
    init_individual_statevars::Function # function to initialize individual state variables
    generate_individual_params::Function # function to initialize individual parameters

    
end

"""
Individual(
        p::CVOrParamStruct, 
        global_statevars::CVOrParamStruct;
        id::Int = 1,
        cohort::Int = 0,
        individual_ode! = default_individual_ODE!,
        individual_rules! = default_individual_rules,
        init_individual_statevars = initialize_individual_statevars,
        )::Individual

Initialization of an individual based on parameters `p` 
and global state `global_statevars`. 

Keyword arguments are used to assure that rules and equations for individual behaviour are inherited correctly.
"""
function Individual(
    p::CVOrParamStruct, 
    global_statevars::CVOrParamStruct;
    id::Int = 1,
    cohort::Int = 0,
    individual_ode! = default_individual_ODE!,
    individual_rules! = default_individual_rules,
    init_individual_statevars = initialize_individual_statevars,
    generate_individual_params = generate_individual_params
    )

    # parameters are instantiated from the individual's method and species-specific parameters
    p_ind = generate_individual_params(p)
    
    # individual stores a reference to global states + copy of own states
    u = ComponentVector(
        glb = global_statevars, # global states
        ind = ComponentVector( # own states
            init_individual_statevars(p_ind, id = id, cohort = cohort);
        )
    )

    # derivatives have the same shape as statevars and start at 0
    du = similar(u)
    du .= 0.


    return Individual(
        du, 
        u,
        p_ind, 
        
        individual_ode!,
        individual_rules!, 
        init_individual_statevars, 
        generate_individual_params

    )
end