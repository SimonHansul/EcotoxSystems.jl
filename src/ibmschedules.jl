# ibmschedules.jl
# generic scheduling for individuals

get_recorded_individual_var_indices(m::AbstractIBM) = map(x -> x in m.recorded_individual_vars,  keys(ind)) |> BitVector

"""
    get_global_statevars!(a::AbstractIndividual, m::AbstractIBM)::Nothing

Retrieve global state variables and derivatives for use in the individual-level model step.
"""
function get_global_statevars!(a::AbstractIndividual, m::AbstractIBM)::Nothing
    
    #a.integrator.du.glb = m.integrator.du.glb
    a.integrator.u.glb = m.integrator.u.glb

    return nothing
end

"""
    set_global_statevars!(m::AbstractIBM, a::AbstractIndividual)::Nothing

Update global state variables and derivatives based on individual-level model step.
"""
function set_global_statevars!(m::AbstractIBM, a::AbstractIndividual)::Nothing

    #m.integrator.du.glb = a.integrator.du.glb 
    m.integrator.u.glb = a.integrator.u.glb

    return nothing
end

"""
    individual_step!(a::Agent, m::Model)::Nothing

The individual-level model step follows a generic pattern:

First the ODE-portion of the model is executed, and the corresponding state variables are updated using the Euler scheme. 

Then the rule-based portion of the model is executed. These are all the functions which cannot / should not be expressed as part of an ODE.
At the minimum, this will include life stage transitions, reproduction and death of individuals. 

For a spatially explicit model, movement should also most likely be part of the rule-based portion, 
as well as functions which require direct information exchange between individuals.
"""
function individual_step!(a::AbstractIndividual, m::AbstractIBM)
    
    get_global_statevars!(a, m) # update reference to global states

    OrdinaryDiffEq.u_modified!(a.integrator, true) # tell the integrator that we might have changed the states
    OrdinaryDiffEq.step!(a.integrator, m.dt, true)
    a.individual_rules!(a, m) # apply the rule-based portion of the individual model

    set_global_statevars!(m, a) # update global states

    return nothing
end

"""
    record_individual!(a::AbstractIndividual, m::AbstractIBM)::Nothing

Store individual-level state variables in `m.individual_record`.
"""
function record_individual!(a::AbstractIndividual, m::AbstractIBM)::Nothing

    if m.record_individuals && isapprox(m.integrator.t % m.saveat, 0, atol = m.dt)
        push!(
            m.individual_record,
            ComponentVector(a.integrator.u.ind; t = m.integrator.t)
        )
    end

    return nothing
end

"""
    record_global!(m::AbstractIBM)::Nothing

Store global state variables in `m.global_record`.
"""
function record_global!(m::AbstractIBM)::Nothing

    if isapprox(m.integrator.t % m.saveat, 0, atol = m.dt)
        push!(m.global_record, ComponentVector(m.integrator.u; t = m.integrator.t))
    end

    return nothing
end


"""
    filter_individuals!(m::AbstractIBM)

Remove individuals which have been flagged to die after the current time-step. 

Individuals for which the condition `u.ind.cause_of_death == 0` applies are retained.
"""
filter_individuals!(m::AbstractIBM) = m.individuals = filter(x -> x.integrator.p.ind.cause_of_death == 0, m.individuals)

function step_all_individuals!(m::AbstractIBM)::Nothing
    
    shuffle!(m.individuals)

    for a in m.individuals
        # before an individual step is executed, the global derivatives are reset to 0. this
        # this is so that individuals can modify the global states, e.g. when ingesting food
        #m.integrator.du.glb .= 0. 
        individual_step!(a, m)
        record_individual!(a, m)
    end

    filter_individuals!(m)

    return nothing
end


"""
    default_global_rules!(m)

Global rule-based portion of the default model. 
"""
function default_global_rules!(m)::Nothing

    i = m.integrator
    i.u.glb.N = length(m.individuals) # tracking population size
    i.u.glb.X = max.(0, m.integrator.u.glb.X) # HOTFIX : negative resource abundances can cause chaos

    return nothing
end


"""
    model_step!(m::AbstractIBM)::Nothing

Generic definition of an individual-based model step.
"""
function model_step!(m::AbstractIBM)::Nothing

    step_all_individuals!(m)
    
    OrdinaryDiffEq.u_modified!(m.integrator, true)
    OrdinaryDiffEq.step!(m.integrator, m.dt, true)
    default_global_rules!(m)
    record_global!(m)

    return nothing
end