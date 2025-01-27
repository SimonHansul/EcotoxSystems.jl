# individuals.jl 
# definition of the default individual and default individual rules
# the mutabe struct for individuals should be generic enough for most applications, 
# the default_individual_rules should be viewed as a starting point for speicifc applications

abstract type AbstractDEBIndividual end

CAUSE_OF_DEATH = Dict(
    0 => "none",
    1 => "age"
)

"""
    default_individual_rules(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing

Defines the default rule-based portion for DEBIndividuals. <br>

The event functions which are used as callbacks during ODE solving are here re-used to apply rules for life stage transitions.

A crude rule for starvation mortality is implemented, applying a constant hazard rate of a certain relative amount of structural mass is lost.

Reproduction is assumed to occur in fixed time intervals, according to `spc.tau_R`.
"""
function default_individual_rules!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing
    
    @unpack glb,ind = a.u
    du = a.du
    p = a.p

    ind.age += m.dt 

    #### life-stage transitions
    # here, we re-use the continuous callback functions defined in EcotoxSystems
    # it happens to be the case that we can treat the individual as an integrator in the callback functions, 
    # since "u" is a field of the individual, just like for integrators

    # check for transition from embryo to juvenile 
    if condition_juvenile(a.u, m.t, a) <= 0
        effect_juvenile!(a)
    end

    # check for transition from juvenile to adult
    if condition_adult(a.u, m.t, a) <= 0
        effect_adult!(a)
    end

    # aging is implemented in a non-mechanistic manner 
    # individuals die when they exceed their maximum age a_max
    # a_max is subject to individual variability
    if ind.age >= a.p.ind.a_max
        ind.cause_of_death = 1.
    end
    
    # for starvation mortality, currently only a limit is set on the amount of mass that can be lost
    # this is basically only a sanity check, and the actual starvation rules should be assessed on a species-by-species basis
    ind.S_max_hist = max(ind.S, ind.S_max_hist)

    if ((ind.S/ind.S_max_hist) < p.ind.S_rel_crit) && (rand() > exp(-p.ind.h_S * m.dt))
    if ((ind.S/ind.S_max_hist) < p.ind.S_rel_crit) && (rand() > exp(-p.ind.h_S * m.dt))
        ind.cause_of_death = 2.
    end

    # reproduction, assuming a constant reproduction period
    
    # reproduction only occurs if the reproduction period has been exceeded
    if ind.time_since_last_repro >= p.ind.tau_R 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg
        let num_offspring = trunc(ind.R / p.ind.X_emb_int)
            for _ in 1:num_offspring
                m.idcount += 1 # increment individual counter
                push!(m.individuals, DEBIndividual( # create new individual and push to individuals vector
                    m.p, 
                    m.u.glb; 
                    id = m.idcount, 
                    cohort = Int(ind.cohort + 1),
                    individual_ode! = a.individual_ode!,
                    individual_rules! = a.individual_rules!,
                    init_individual_statevars = a.init_individual_statevars,
                    gen_ind_params = a.generate_individual_params,
                    )
                )
                ind.R -= p.ind.X_emb_int # decrease reproduction buffer
                ind.cum_repro += 1 # keep track of cumulative reproduction of the mother individual
            end
            ind.time_since_last_repro = 0. # reset reproduction period
        end
    # if reproduction period has not been exceeded,
    else
        ind.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end


@with_kw mutable struct DEBIndividual <: AbstractDEBIndividual

    du::ComponentVector
    u::ComponentVector
    p::ComponentVector

    individual_ode!::Function # equation-based portion of the individual step
    individual_rules!::Function # rule-based portion of the individual step
    init_individual_statevars::Function # function to initialize individual state variables
    generate_individual_params::Function # function to initialize individual parameters

    """
    DEBIndividual(
            p::ComponentVector, 
            global_statevars::ComponentVector;
            id::Int = 1,
            cohort::Int = 0,
            individual_ode! = DEBODE_individual!,
            individual_rules! = default_individual_rules,
            init_individual_statevars = initialize_individual_statevars,
            )::AbstractDEBIndividual

    Initialization of an individual based on parameters `p` 
    and global state `global_statevars`. 

    Keyword arguments are used to assure that rules and equations for individual behaviour are inherited correctly.
    """
    function DEBIndividual(
        p::ComponentVector, 
        global_statevars::ComponentVector;
        id::Int = 1,
        cohort::Int = 0,
        individual_ode! = DEBODE_individual!,
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