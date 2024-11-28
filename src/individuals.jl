abstract type AbstractDEBIndividual end

CAUSE_OF_DEATH = Dict(
    0 => "none",
    1 => "age"
)


"""
    default_individual_rules(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing

Defines the default rule-based portion for DEBIndividuals. <br>

The event functions which are used as callbacks during ODE solving are here re-used to apply rules for life stage transitions.
A crude rule for starvation mortality is implemented, applying a dependency between 
    starvation mortality and the scaled functional response.

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
        ind.cause_of_death = 1
    end
    
    # for starvation mortality, we assume that mortality can occur as soon as the scaled functional response falls below a threshold value f_Xthr
    # below that threshold value, survival probability s_f decreases 
    # s_min is the intercept of s_f vs f_X, 
    # and so related to the survival curve we expect if an individual is depleted from food
    # by calling sig(), we express this if/else statement as a continuous function

    # this is a crude rule, but can serve as a starting point

    #FIXME: this function somehow takes the most computation time of all in the individual step

    ind.S_max_hist = max(ind.S, ind.S_max_hist)

    if ((ind.S/ind.S_max_hist) < p.ind.S_rel_crit) && (rand() <= exp(-p.ind.h_S * m.dt))
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
                    individual_rules! = a.individual_rules!
                    )
                )
                ind.R -= p.ind.X_emb_int # decrease reproduction buffer
                ind.cum_repro += 1
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

    individual_ode!::Function
    individual_rules!::Function
    statevar_init::Function

    function DEBIndividual(
        p::ComponentVector, 
        global_statevars::ComponentVector;
        id::Int = 1,
        cohort::Int = 0,
        individual_ode! = DEBODE_individual!,
        individual_rules! = default_individual_rules,
        statevar_init = initialize_individual_statevars
        )
        
        a = new() 

        a.p = generate_individual_params(p)
        a.individual_ode! = individual_ode!
        a.individual_rules! = individual_rules!
        
        a.u = ComponentVector(
            glb = global_statevars,
            ind = ComponentVector(
                statevar_init(a.p, id = id, cohort = cohort);
            )
        )

        a.du = similar(a.u)
        a.du .= 0.
        
        return a
    end
end