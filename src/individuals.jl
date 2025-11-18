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

"""
    default_individual_rules(a::AbstractIndividual, m::AbstractDEBIBM)::Nothing

Defines the default rule-based portion for DEBIndividuals. <br>

The event functions which are used as callbacks during ODE solving are here re-used to apply rules for life stage transitions.

A crude rule for starvation mortality is implemented, applying a constant hazard rate of a certain relative amount of structural mass is lost.

Reproduction is assumed to occur in fixed time intervals, according to `spc.tau_R`.
"""
function default_individual_rules!(a::AbstractIndividual, m::AbstractDEBIBM)::Nothing

    p = a.p
    u = a.u
    u.ind.age += m.dt 

 
    # aging is implemented in a non-mechanistic manner 
    # individuals die when they exceed their maximum age a_max
    # a_max is subject to individual variability
    if death_by_aging(u.ind[:age], p[:ind][:a_max])
        u.ind.cause_of_death = 1.
    end
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    u.ind.S_max_hist = determine_S_max_hist(u.ind[:S], u.ind[:S_max_hist])

    if death_by_loss_of_structure(u.ind[:S], u.ind[:S_max_hist], p[:ind][:W_S_rel_crit], p[:ind][:h_S], m.dt)
        u.ind.cause_of_death = 2.
    end

    # reproduction, assuming a constant reproduction period
    
    # reproduction only occurs if the reproduction period has been exceeded
    if check_reproduction_period(u.ind[:time_since_last_repro], p[:ind][:tau_R]) 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg
        for _ in 1:calc_num_offspring(u.ind[:R], p[:ind][:X_emb_int])
            m.idcount += 1 # increment individual counter
            push!(m.individuals, Individual( # create new individual and push to individuals vector
                m.p, 
                m.u.glb; 
                id = m.idcount, 
                cohort = Int(u.ind.cohort + 1),
                individual_ode! = a.individual_ode!,
                individual_rules! = a.individual_rules!,
                init_individual_statevars = a.init_individual_statevars,
                gen_ind_params = a.generate_individual_params,
                )
            )
            u.ind.R -= p[:ind][:X_emb_int] # decrease reproduction buffer
            u.ind.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
        end
        u.ind.time_since_last_repro = 0. # reset reproduction period
    # if reproduction period has not been exceeded,
    else
        u.ind.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end


@with_kw mutable struct Individual <: AbstractIndividual

    du::CVOrParamStruct
    u::CVOrParamStruct
    p::CVOrParamStruct

    individual_ode!::Function # equation-based portion of the individual step
    individual_rules!::Function # rule-based portion of the individual step
    init_individual_statevars::Function # function to initialize individual state variables
    generate_individual_params::Function # function to initialize individual parameters

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
        gen_ind_params = generate_individual_params
        )
        
        a = new() 

        a.individual_ode! = individual_ode!
        a.individual_rules! = individual_rules!
        a.init_individual_statevars = init_individual_statevars
        a.generate_individual_params = gen_ind_params
        a.p = a.generate_individual_params(p)
        
        # individual stores a reference to global states + copy of own states
        a.u = ComponentVector(
            glb = global_statevars, # global states
            ind = ComponentVector( # own states
                a.init_individual_statevars(a.p, id = id, cohort = cohort);
            )
        )

        a.du = similar(a.u)
        a.du .= 0.
        
        return a
    end
end