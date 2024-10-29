#schedules.jl

get_recorded_individual_var_indices(m::AbstractDEBIBM) = map(x -> x in m.recorded_individual_vars,  keys(ind)) |> BitVector

"""
    get_global_statevars!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing


"""
function get_global_statevars!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing
    
    a.du.glb = m.du.glb
    a.u.glb = m.du.glb

    return nothing
end


"""
    set_global_statevars!(m::AbstractDEBIBM, a::AbstractDEBIndividual)::Nothing

Updates global state variables and derivatives according to modifications by an individual.
"""
function set_global_statevars!(m::AbstractDEBIBM, a::AbstractDEBIndividual)::Nothing

    m.du.glb = a.du.glb 
    m.u.glb = a.u.glb

    return nothing
end


"""
    individual_step_rulebased!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing

Defines the rule-based portion of an individual step. 

The life stage callbacks defined in the EnergyBudgetDiffEqs are used to apply rules for life stage transitions.

A crude rule for starvation mortality is implemented, applying a constant hazard rqte 
when they lose a given fraction of their structural mass. 

Reproduction is assumed to occur in fixed time intervals, according to `spc.tau_R`.
"""
function individual_step_rulebased!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing
    
    @unpack glb,ind = ind
    du = a.du
    p = a.p

    ind.age += m.dt 

    #### life-stage transitions
    # here, we re-use the continuous callback functions defined in EnergyBudgetDiffEqs
    # it happens to be the case that we can treat the individual as an integrator in the callback functions, 
    # since "u" is a field of the individual, just like for integrators

    # check for transition from embryo to juvenile 
    if condition_juvenile(ind, m.t, a) <= 0
        effect_juvenile!(a)
    end

    # check for transition from juvenile to adult
    if condition_adult(ind, m.t, a) <= 0
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
    # by calling sig(), we express this if/else statement as a continuous function
    
    let s_f = sig(
        ind.f_X, p.ind.f_Xthr, 
        (1 - p.ind.s_min) * ind.f_X / p.ind.f_Xthr + p.ind.s_min,
        1.)^m.dt

        if rand() > s_f
            ind.cause_of_death = 2.
        end
    end

    # reproduction, assuming a constant reproduction period
    
    # reproduction only occurs if the reproduction period has been exceeded
    if ind.time_since_last_repro >= p.ind.tau_R 
        # if that is the case, calculate the number of offspring, 
        # based on the reproduction buffer and the dry mass of an egg
        let num_offspring = trunc(ind.R / ind.X_emb_int)
            for _ in 1:num_offspring
                m.u.idcount += 1 # increment individual counter
                push!(m.individuals, DEBIndividual( # create new individual and push to individuals vector
                    p, 
                    m.u, 
                    m.idcount; 
                    cohort = ind.cohort + 1
                    )
                )
                ind.R -= ind.X_emb_int # decrease reproduction buffer
            end
            ind.time_since_last_repro = 0. # reset reproduction period
        end
    # if reproduction period has not been exceeded,
    else
        ind.time_since_last_repro += m.dt # track reproduction period
    end

    return nothing
end


"""
    individual_step!(a::Agent, m::Model)::Nothing

The individual step follows a generic pattern:

First the ODE-portion of the model is executed, and the corresponding state variables are updated using the Euler scheme. 

Then the rule-based portion of the model is executed. These are all the functions which cannot / should not be expressed as part of an ODE.
At the minimum, this will include life stage transitions, reproduction and death of individuals. 
For a spatially explicit model, movement should also most likely be part of the rule-based portion, 
as well as functions which require direct information exchange between individuals.
"""
function individual_step!(a::AbstractDEBIndividual, m::AbstractDEBIBM)
    
    get_global_statevars!(a, m) # update reference to global states

    a.individual_ode!(a.du, a.u, a.p, m.t) # calculate derivatives of the ODE-portion of the individual model
    Euler!(a.u, a.du, m.dt) # update states using Euler scheme -> this could even be replaced with a more sophisticated scheme
    a.individual_rules!(a, m) # apply the rule-based portion of the individual model

    set_global_statevars!(m, a) # update global states

    return nothing
end


"""
    Euler!(u::ComponentVector, du::ComponentVector, dt::Real, statevar_indices::Vector{Int})::Nothing

Applying the Euler scheme to state variables. 

args

- `u`: State variables
- `du`: Derivatives
- `dt`: Timestep

"""
function Euler!(u::ComponentVector, du::ComponentVector, dt::Real)::Nothing
    u .+= du .* dt
    return nothing
end

"""
    record_individual!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing

Pushes the individual state variables to `m.individual_record` in fixed time intervals, according to `m.saveat`.
"""
function record_individual!(a::AbstractDEBIndividual, m::AbstractDEBIBM)::Nothing

    if m.record_individuals && isapprox(m.t % m.saveat, 0, atol = m.dt)
        push!(
            m.individual_record,
            ComponentVector(a.u.ind; t = m.t)
        )
    end

    return nothing
end


"""
    filter_individuals!(m::AbstractDEBIBM)

Remove dead individuals from a model. 
Agents are flagged to be removed if `cause_of_death>0`.
"""
filter_individuals!(m::AbstractDEBIBM) = m.individuals = filter(x -> x.u.ind.cause_of_death == 0, m.individuals)


"""
    step_all_individuals!(m::AbstractDEBIBM)::Nothing

Executes individual_step! for all individuals in the model. 
Records individual states.
Filters individuals vector to remove dead indiviualds. 
Agents which die will be recorded a last time before they are removed.
"""
function step_all_individuals!(m::AbstractDEBIBM)::Nothing
    
    shuffle!(m.individuals)

    for a in m.individuals
        individual_step!(a, m)
        record_individual!(a, m)
    end
    filter_individuals!(m)

    return nothing
end


"""
    model_step!(m::AbstractDEBIBM)::Nothing

Execute a single model step of a IndividualBasedModel.
"""
function model_step!(m::AbstractDEBIBM)::Nothing
    # calculate global derivatives
    # change in resource abundance, chemical stressor exposure etc.
    
    m.global_ode!(m.du, m.u, m.p, m.t)
    step_all_individuals!(m)
    
    # global statevars are updated after individual derivatives are calculated
    # this is important because individuals affect global states using mutating operators
    Euler!(m.u, m.du, m.dt) 
    m.u.glb.X_p = max(0, m.u.glb.X_p) # HOTFIX : negative resource abundances can cause chaos

    m.t += m.dt

    return nothing
end